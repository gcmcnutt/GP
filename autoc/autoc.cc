
// autoc.cc

/* -------------------------------------------------------------------
From skeleton/skeleton.cc
------------------------------------------------------------------- */

#include <boost/asio.hpp>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <new>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>

#include "gp.h"
#include "gpconfig.h"
#include "autoc.h"
#include "minisim.h"
#include "pathgen.h"
#include "renderer.h"
#include "threadpool.h"
#include "logger.h"

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/support/date_time.hpp>

using namespace std;
namespace logging = boost::log;

Logger logger;

// Define configuration parameters and the neccessary array to
// read/write the configuration to a file.  If you need more
// variables, just add them below and insert an entry in the
// configArray.
GPVariables cfg;
ExtraConfig extraCfg;
struct GPConfigVarInformation configArray[] =
{
  {"PopulationSize", DATAINT, &cfg.PopulationSize},
  {"NumberOfGenerations", DATAINT, &cfg.NumberOfGenerations},
  {"CreationType", DATAINT, &cfg.CreationType},
  {"CrossoverProbability", DATADOUBLE, &cfg.CrossoverProbability},
  {"CreationProbability", DATADOUBLE, &cfg.CreationProbability},
  {"MaximumDepthForCreation", DATAINT, &cfg.MaximumDepthForCreation},
  {"MaximumDepthForCrossover", DATAINT, &cfg.MaximumDepthForCrossover},
  {"SelectionType", DATAINT, &cfg.SelectionType},
  {"TournamentSize", DATAINT, &cfg.TournamentSize},
  {"DemeticGrouping", DATAINT, &cfg.DemeticGrouping},
  {"DemeSize", DATAINT, &cfg.DemeSize},
  {"DemeticMigProbability", DATADOUBLE, &cfg.DemeticMigProbability},
  {"SwapMutationProbability", DATADOUBLE, &cfg.SwapMutationProbability},
  {"ShrinkMutationProbability", DATADOUBLE, &cfg.ShrinkMutationProbability},
  {"AddBestToNewPopulation", DATAINT, &cfg.AddBestToNewPopulation},
  {"SteadyState", DATAINT, &cfg.SteadyState},
  {"SimNumPathsPerGeneration", DATAINT, &extraCfg.simNumPathsPerGen},
  {"EvalThreads", DATAINT, &extraCfg.evalThreads},
  {"MinisimProgram", DATASTRING, &extraCfg.minisimProgram},
  {"MinisimPortOverride", DATAINT, &extraCfg.minisimPortOverride},
  {"", DATAINT, NULL}
};


// Define function and terminal identifiers
enum Operators {
  ADD = 0, SUB, MUL, DIV,
  IF, EQ, GT,
  SIN, COS,
  GETDPHI, GETDTHETA, GETDTARGET, GETDHOME, GETVEL,
  GETPITCH, GETROLL, GETTHROTTLE,
  SETPITCH, SETROLL, SETTHROTTLE,
  PI, ZERO, ONE, TWO, PROGN, _END
};
const int OPERATORS_NR_ITEM = _END;


// Define class identifiers
const int MyGeneID = GPUserID;
const int MyGPID = GPUserID + 1;
const int MyPopulationID = GPUserID + 2;

ThreadPool* threadPool;
std::atomic_bool printEval = false; // verbose (used for rendering best of population)
std::ofstream fout;
Renderer renderer(extraCfg);
std::mutex evalMutex;

class MyGP;

// Inherit the three GP classes GPGene, GP and GPPopulation
class MyGene : public GPGene
{
public:
  // Duplication (mandatory)
  MyGene(const MyGene& gpo) : GPGene(gpo) { }
  virtual GPObject& duplicate() { return *(new MyGene(*this)); }

  // Creation of own class objects (mandatory)
  MyGene(GPNode& gpo) : GPGene(gpo) {}
  virtual GPGene* createChild(GPNode& gpo) {
    return new MyGene(gpo);
  }

  // Tree evaluation (not mandatory, but somehow the trees must be
  // parsed to evaluate the fitness)
  double evaluate(std::vector<Path>& path, MyGP& gp, double arg);

  // Load and save (not mandatory)
  MyGene() {}
  virtual int isA() { return MyGeneID; }
  virtual GPObject* createObject() { return new MyGene; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access children (not mandatory)
  MyGene* NthMyChild(int n) {
    return (MyGene*)GPContainer::Nth(n);
  }
};



class MyGP : public GP
{
public:
  // Duplication (mandatory)
  MyGP(MyGP& gpo) : GP(gpo) { }
  virtual GPObject& duplicate() { return *(new MyGP(*this)); }

  // Creation of own class objects (mandatory)
  MyGP(int genes) : GP(genes) {}
  virtual GPGene* createGene(GPNode& gpo) {
    return new MyGene(gpo);
  }

  // Tree evaluation (mandatory)
  virtual void evaluate();

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Load and save (not mandatory)
  MyGP() {}
  virtual int isA() { return MyGPID; }
  virtual GPObject* createObject() { return new MyGP; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  // Access trees (not mandatory)
  MyGene* NthMyGene(int n) {
    return (MyGene*)GPContainer::Nth(n);
  }

  // async evaluator
  void evalTask(WorkerContext& context);

  AircraftState aircraftState{ 0, Eigen::Quaterniond::Identity(), Eigen::Vector3d(0, 0, 0), 0.0, 0.0, 0.0, 0, false };
  long pathIndex = 0; // current entry on path
};

std::atomic_ulong nanDetector(0);
std::vector<MyGP*> tasks = std::vector<MyGP*>();

class MyPopulation : public GPPopulation
{
public:
  // Constructor (mandatory)
  MyPopulation(GPVariables& GPVar_, GPAdfNodeSet& adfNs_) :
    GPPopulation(GPVar_, adfNs_) {}

  // Duplication (mandatory)
  MyPopulation(MyPopulation& gpo) : GPPopulation(gpo) {}
  virtual GPObject& duplicate() { return *(new MyPopulation(*this)); }

  // Creation of own class objects (mandatory)
  virtual GP* createGP(int numOfTrees) { return new MyGP(numOfTrees); }

  // Load and save (not mandatory)
  MyPopulation() {}
  virtual int isA() { return MyPopulationID; }
  virtual GPObject* createObject() { return new MyPopulation; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  virtual void endOfEvaluation() {
    // dispatch all the GPs now (TODO this may still work inline with evaluate)
    for (auto& task : tasks) {
      threadPool->enqueue([task](WorkerContext& context) {
        task->evalTask(context);
        });
    }

    // wait for all tasks to finish
    threadPool->wait_for_tasks();
    tasks.clear();
  }

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access genetic programs (not mandatory)
  MyGP* NthMyGP(int n) {
    return (MyGP*)GPContainer::Nth(n);
  }
};

int getIndex(std::vector<Path>& path, MyGP& gp, double arg) {
  if (isnan(arg)) {
    return gp.pathIndex;
  }

  // a range of steps to check, can't go lower than the beginning index
  // TODO this checks path next, not the actual simulation steps...
  // XXX for now, this allows forecasting the future path
  int idx = std::clamp((int)arg, -5, 5) + gp.pathIndex;
  idx = std::clamp(idx, 0, (int)path.size() - 1);
  return idx;
}

// This function evaluates the fitness of a genetic tree.  We have the
// freedom to define this function in any way we like.  
double MyGene::evaluate(std::vector<Path>& path, MyGP& run, double arg)
{
  double returnValue = 0.0;

  switch (node->value())
  {
  case ADD: returnValue = NthMyChild(0)->evaluate(path, run, arg) + NthMyChild(1)->evaluate(path, run, arg); break;
  case SUB: returnValue = -NthMyChild(0)->evaluate(path, run, arg) - NthMyChild(1)->evaluate(path, run, arg); break;
  case MUL: returnValue = NthMyChild(0)->evaluate(path, run, arg) * NthMyChild(1)->evaluate(path, run, arg); break;
  case DIV: {
    double dividend = NthMyChild(0)->evaluate(path, run, arg);
    double divisor = NthMyChild(1)->evaluate(path, run, arg);
    returnValue = (divisor == 0) ? 0 : dividend / divisor;
    break;
  }
  case IF: returnValue = NthMyChild(0)->evaluate(path, run, arg) ? NthMyChild(1)->evaluate(path, run, arg) : NthMyChild(2)->evaluate(path, run, arg); break;
  case EQ: returnValue = NthMyChild(0)->evaluate(path, run, arg) == NthMyChild(1)->evaluate(path, run, arg); break;
  case GT: returnValue = NthMyChild(0)->evaluate(path, run, arg) > NthMyChild(1)->evaluate(path, run, arg); break;
  case SETPITCH: returnValue = run.aircraftState.pitchCommand = NthMyChild(0)->evaluate(path, run, arg); break;
  case SETROLL: returnValue = run.aircraftState.rollCommand = NthMyChild(0)->evaluate(path, run, arg); break;
  case SETTHROTTLE: returnValue = run.aircraftState.throttleCommand = NthMyChild(0)->evaluate(path, run, arg); break;
  case GETPITCH: returnValue = run.aircraftState.pitchCommand; break;
  case GETROLL: returnValue = run.aircraftState.rollCommand; break;
  case GETTHROTTLE: returnValue = run.aircraftState.throttleCommand; break;
  case SIN: returnValue = sin(NthMyChild(0)->evaluate(path, run, arg)); break;
  case COS: returnValue = cos(NthMyChild(0)->evaluate(path, run, arg)); break;
  case PI: returnValue = M_PI; break;
  case ZERO: returnValue = 0; break;
  case ONE: returnValue = 1; break;
  case TWO: returnValue = 2; break;
  case PROGN: {
    NthMyChild(0)->evaluate(path, run, arg);
    returnValue = NthMyChild(1)->evaluate(path, run, arg);
    break;
  }
  case GETVEL: returnValue = run.aircraftState.dRelVel; break;

  case GETDPHI: // compute roll goal from current to target
  {
    // Calculate the vector from craft to target in world frame
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    Eigen::Vector3d craftToTarget = path.at(idx).start - run.aircraftState.position;

    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = run.aircraftState.aircraft_orientation.inverse() * craftToTarget;

    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());

    // Calculate the angle between the projected vector and the body Z-axis
    returnValue = std::atan2(projectedVector.y(), -projectedVector.z());
    break;
  }

  case GETDTHETA: // compute pitch goal from current to target
  {
    // Calculate the vector from craft to target in world frame
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    Eigen::Vector3d craftToTarget = path.at(idx).start - run.aircraftState.position;

    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = run.aircraftState.aircraft_orientation.inverse() * craftToTarget;

    // Calculate pitch angle
    returnValue = std::atan2(-target_local.z(), target_local.x());
    break;
  }

  case GETDTARGET: // get distance to the next point
  {
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    returnValue = (path.at(idx).start - run.aircraftState.position).norm();
    break;
  }

  case GETDHOME: // get distance to the home point
  {
    returnValue = (Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE) - run.aircraftState.position).norm();
    break;
  }

  default:
    GPExitSystem("MyGene::evaluate", "Undefined node value");
  }

#define RANGELIMIT 1000000
  if (returnValue < -RANGELIMIT)
    return -RANGELIMIT;
  if (returnValue > RANGELIMIT)
    return RANGELIMIT;
  if (abs(returnValue) < 0.000001)
    return 0;
  return returnValue;
}

void MyGP::evaluate()
{
  tasks.push_back(this);
}

std::string crashReasonToString(CrashReason type) {
  switch (type) {
  case CrashReason::Sim: return "Sim";
  case CrashReason::Eval: return "Eval";
  case CrashReason::None: return "None";
  default: return "*?*";
  }
}

// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evalTask(WorkerContext& context)
{
  stdFitness = 0;

  for (int i = 0; i < renderer.generationPaths.size(); i++) {
    auto& path = renderer.generationPaths.at(i);
    double localFitness = 0;
    std::vector<Eigen::Vector3d> planPath = std::vector<Eigen::Vector3d>();
    std::vector<Eigen::Vector3d> actualPath = std::vector<Eigen::Vector3d>();

    // deal with pre.path on the initial eval..
    if (path.size() == 0) {
      stdFitness = 1000001;
      continue;
    }

    // random initial orientation
    Eigen::Quaterniond aircraft_orientation;
    Eigen::Vector3d initialPosition;
    {
      // Guard random functions
      std::lock_guard<std::mutex> guard(evalMutex);
      aircraft_orientation = Eigen::Quaterniond::UnitRandom();

      initialPosition = Eigen::Vector3d((((double)GPrand() / RAND_MAX) - 0.5) * SIM_INITIAL_LOCATION_DITHER,
        (((double)GPrand() / RAND_MAX) - 0.5) * SIM_INITIAL_LOCATION_DITHER,
        SIM_INITIAL_ALTITUDE - ((double)GPrand() / RAND_MAX) * SIM_INITIAL_LOCATION_DITHER);
    }

    // wait for generic request (we ignore)
    receiveRPC<AircraftState>(*context.socket);

    // send initial state reset
    aircraftState = AircraftState{ SIM_INITIAL_VELOCITY, aircraft_orientation, initialPosition, 0.0, 0.0, SIM_INITIAL_THROTTLE, 0, false };
    sendRPC(*context.socket, MainToSim{ ControlType::AIRCRAFT_STATE, AircraftState{aircraftState} });

    // iterate the simulator
    unsigned long int duration_msec = 0; // how long have we been running
    pathIndex = 0; // where are we on the path?
    bool printHeader = true;

    // as long as we are within the time limit and have not reached the end of the path
    CrashReason crashReason = CrashReason::None;

    // error accumulators
    double distance_error_sum = 0;
    double angle_error_sum = 0;
    double control_smoothness_sum = 0;
    int simulation_steps = 0;

    // initial states
    double roll_prev = aircraftState.rollCommand;
    double pitch_prev = aircraftState.pitchCommand;
    double throttle_prev = aircraftState.throttleCommand;

    while (duration_msec < SIM_TOTAL_TIME_MSEC && pathIndex < path.size() - 2 && crashReason == CrashReason::None) {
      // wait for next state update from sim
      aircraftState = receiveRPC<AircraftState>(*context.socket);

      // did sim crash?
      if (aircraftState.simCrashed) {
        crashReason = CrashReason::Sim;
      }

      // ---------------------- first how did we do?
      // Compute the distance between the aircraft and the goal
      double distanceFromGoal = (path.at(pathIndex).start - aircraftState.position).norm();
      // normalize [100:0]
      distanceFromGoal = distanceFromGoal * 100.0 / SIM_PATH_RADIUS_LIMIT;

      // Compute vector from me to target
      Eigen::Vector3d target_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start);
      Eigen::Vector3d aircraft_to_target = (path.at(pathIndex).start - aircraftState.position);
      double dot_product = target_direction.dot(aircraft_to_target);
      double angle_rad = std::acos(std::clamp(dot_product / (target_direction.norm() * aircraft_to_target.norm()), -1.0, 1.0));
      // normalize [100:0]
      angle_rad = angle_rad * 100.0 / M_PI;

      // control smoothness -- internal vector distance
      double smoothness = pow(roll_prev - aircraftState.rollCommand, 2.0);
      smoothness += pow(pitch_prev - aircraftState.pitchCommand, 2.0);
      smoothness += pow(throttle_prev - aircraftState.throttleCommand, 2.0);
      smoothness = sqrt(smoothness);
      // normalize [100:0]
      smoothness = smoothness * 100.0 / 3.46;

      // ready for next cycle
      roll_prev = aircraftState.rollCommand;
      pitch_prev = aircraftState.pitchCommand;
      throttle_prev = aircraftState.throttleCommand;

      // but eval detected crashed outside the cylinder? (low elevation is detected from sim)
      double distanceFromOrigin = std::sqrt(aircraftState.position[0] * aircraftState.position[0] +
        aircraftState.position[1] * aircraftState.position[1]);
      if (aircraftState.position[2] < (SIM_MIN_ELEVATION - SIM_PATH_RADIUS_LIMIT) || // too high
        aircraftState.position[2] > (SIM_MIN_ELEVATION) || // too low
        distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) { // too far
        crashReason = CrashReason::Eval;
      }

      if (printEval) {
        if (printHeader) {
          fout << "  Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP controlP\n";
          printHeader = false;
        }

        // convert aircraft_orientaton to euler
        Eigen::Vector3d euler = aircraftState.aircraft_orientation.toRotationMatrix().eulerAngles(2, 1, 0);

        char outbuf[1000]; // XXX use c++20
        sprintf(outbuf, "%06ld %3ld % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n",
          aircraftState.simTime, pathIndex,
          path.at(pathIndex).distanceFromStart,
          path.at(pathIndex).start[0],
          path.at(pathIndex).start[1],
          path.at(pathIndex).start[2],
          aircraftState.position[0],
          aircraftState.position[1],
          aircraftState.position[2],
          euler[2],
          euler[1],
          euler[0],
          aircraftState.dRelVel,
          aircraftState.rollCommand,
          aircraftState.pitchCommand,
          aircraftState.throttleCommand,
          distanceFromGoal,
          angle_rad,
          smoothness
        );
        fout << outbuf;

        // now update points for Renderer
        planPath.push_back(path.at(pathIndex).start); // XXX push the full path
        actualPath.push_back(aircraftState.position);
      }

      distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
      angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
      control_smoothness_sum += pow(smoothness, FITNESS_CONTROL_WEIGHT);
      simulation_steps++;

      // TODO ensure time is forward
      unsigned long int duration_msec = aircraftState.simTime;

      // search for location of next timestamp
      double timeDistance = SIM_RABBIT_VELOCITY * duration_msec / 1000.0;
      while (pathIndex < path.size() - 2 && (path.at(pathIndex).distanceFromStart < timeDistance)) {
        pathIndex++;
      }

      // ESTIMATE: pre-compute the estimated roll, pitch, throttle in prep for the evaluator

      // *** ROLL: Calculate the vector from craft to target in world frame
      Eigen::Vector3d craftToTarget = path.at(pathIndex).start - aircraftState.position;

      // Transform the craft-to-target vector to body frame
      Eigen::Vector3d target_local = aircraftState.aircraft_orientation.inverse() * craftToTarget;

      // Project the craft-to-target vector onto the body YZ plane
      Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());

      // Calculate the angle between the projected vector and the body Z-axis 
      double rollEstimate = std::atan2(projectedVector.y(), -projectedVector.z());

      // *** PITCH: Calculate the vector from craft to target in world frame if it did rotate
      Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
      Eigen::Quaterniond virtualOrientation = aircraftState.aircraft_orientation * rollRotation;

      // Transform target vector to new virtual orientation
      Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;

      // Calculate pitch angle
      double pitchEstimate = std::atan2(-newLocalTargetVector.z(), newLocalTargetVector.x());

      // // now try to determine if pitch up or pitch down makes more sense
      // if (std::abs(pitchEstimate) > M_PI / 2) {
      //   pitchEstimate = (pitchEstimate > 0) ? pitchEstimate - M_PI : pitchEstimate + M_PI;
      //   rollEstimate = -rollEstimate;
      // }

      aircraftState.rollCommand = rollEstimate / M_PI;
      aircraftState.pitchCommand = pitchEstimate / M_PI;

      // *** Throttle estimate
      {
        double distance = (path.at(pathIndex).start - aircraftState.position).norm();
        double throttleEstimate = std::clamp((distance - 10) / aircraftState.dRelVel, -1.0, 1.0);
        aircraftState.throttleCommand = throttleEstimate;
      }

      // GP determine control input adjustments
      NthMyGene(0)->evaluate(path, *this, 0);

      // reply with the commands
      MainToSim mainToSim{ ControlType::CONTROL_SIGNAL, ControlSignal{aircraftState.pitchCommand, aircraftState.rollCommand, aircraftState.throttleCommand} };
      sendRPC(*context.socket, mainToSim);
    }

    // tally up the normlized fitness based on steps and progress
    double normalized_distance_error = (distance_error_sum / simulation_steps);
    double normalized_velocity_align = (angle_error_sum / simulation_steps);
    double normalized_control_smoothness = (control_smoothness_sum / simulation_steps);
    localFitness = normalized_distance_error + normalized_velocity_align + normalized_control_smoothness;

    if (isnan(localFitness)) {
      nanDetector++;
    }

    if (crashReason != CrashReason::None) {
      double fractional_distance_remaining = 1.0 - path.at(pathIndex).distanceFromStart / path.back().distanceFromStart;
      localFitness += SIM_CRASH_PENALTY * fractional_distance_remaining;
    }

    if (printEval) {
      // now update actual list of lists
      std::vector<Path> actualElementList;
      for (auto& actualElement : actualPath) {
        actualElementList.push_back({ actualElement, 0 });
      }

      // now update plan list of lists
      std::vector<Path> planElementList;
      for (auto& planElement : planPath) {
        planElementList.push_back({ planElement, 0 });
      }
      renderer.addPathElementList(planElementList, actualElementList);

      planPath.clear();
      actualPath.clear();
    }

    stdFitness += localFitness;

    *logger.debug() << "MyGP: " << this << " path[" << i << "] complete." << endl;
  }

  // normalize
  stdFitness /= renderer.generationPaths.size();

  if (printEval) {
    renderer.update();
  }
}


// Create function and terminal set
void createNodeSet(GPAdfNodeSet& adfNs)
{
  // Reserve space for the node sets
  adfNs.reserveSpace(1);

  // Now define the function and terminal set for each ADF and place
  // function/terminal sets into overall ADF container
  GPNodeSet& ns = *new GPNodeSet(OPERATORS_NR_ITEM);

  adfNs.put(0, ns);

  // Define functions/terminals and place them into the appropriate
  // sets.  Terminals take two arguments, functions three (the third
  // parameter is the number of arguments the function has)
  ns.putNode(*new GPNode(ADD, "ADD", 2));
  ns.putNode(*new GPNode(SUB, "SUB", 2));
  ns.putNode(*new GPNode(MUL, "MUL", 2));
  ns.putNode(*new GPNode(DIV, "DIV", 2));
  ns.putNode(*new GPNode(IF, "IF", 3));
  ns.putNode(*new GPNode(EQ, "EQ", 2));
  ns.putNode(*new GPNode(GT, "GT", 2));
  ns.putNode(*new GPNode(SETPITCH, "SETPITCH", 1));
  ns.putNode(*new GPNode(SETROLL, "SETROLL", 1));
  ns.putNode(*new GPNode(SETTHROTTLE, "SETTHROTTLE", 1));
  ns.putNode(*new GPNode(GETPITCH, "GETPITCH"));
  ns.putNode(*new GPNode(GETROLL, "GETROLL"));
  ns.putNode(*new GPNode(GETTHROTTLE, "GETTHROTTLE"));
  ns.putNode(*new GPNode(SIN, "SIN", 1));
  ns.putNode(*new GPNode(COS, "COS", 1));
  ns.putNode(*new GPNode(PI, "PI"));
  ns.putNode(*new GPNode(ZERO, "0"));
  ns.putNode(*new GPNode(ONE, "1"));
  ns.putNode(*new GPNode(TWO, "2"));
  ns.putNode(*new GPNode(PROGN, "PROGN", 2));
  ns.putNode(*new GPNode(GETDPHI, "GETDPHI", 1));
  ns.putNode(*new GPNode(GETDTHETA, "GETDTHETA", 1));
  ns.putNode(*new GPNode(GETDTARGET, "GETDTARGET", 1));
  ns.putNode(*new GPNode(GETVEL, "GETVEL"));
  ns.putNode(*new GPNode(GETDHOME, "GETDHOME"));
}



void newHandler()
{
  cerr << "\nFatal error: Out of memory." << endl;
  exit(1);
}


int main()
{
  logging::add_console_log(
    std::cout,
    boost::log::keywords::format = (
      boost::log::expressions::stream
      << boost::log::expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
      << ": <" << logging::trivial::severity
      << "> " << boost::log::expressions::smessage
      )
  );
  logging::core::get()->set_filter(
    logging::trivial::severity >= logging::trivial::info
  );

  logging::add_common_attributes();

  logger = Logger();

  // Set up a new-handler, because we might need a lot of memory, and
  // we don't know it's there.
  set_new_handler(newHandler);

  // Init GP system.
  GPInit(1, -1);

  // Read configuration file.
  GPConfiguration config(*logger.info(), "autoc.ini", configArray);
  renderer.extraCfg = extraCfg;

  // initialize workers
  threadPool = new ThreadPool(extraCfg);

  // Print the configuration
  *logger.info() << cfg << endl;
  *logger.info() << "SimNumPathsPerGen: " << extraCfg.simNumPathsPerGen << endl;
  *logger.info() << "EvalThreads: " << extraCfg.evalThreads << endl;
  *logger.info() << "MinisimProgram: " << extraCfg.minisimProgram << endl;
  *logger.info() << "MinisimPortOverride: " << extraCfg.minisimPortOverride << endl << endl;

  // Create the adf function/terminal set and print it out.
  GPAdfNodeSet adfNs;
  createNodeSet(adfNs);
  *logger.info() << adfNs << endl;

  // Open the main output file for the data and statistics file.
  // First set up names for data file.  Remember we should delete the
  // string from the stream, well just a few bytes
  ostringstream strOutFile, strStatFile;
  strOutFile << "data.dat" << ends;
  strStatFile << "data.stc" << ends;
  fout.open(strOutFile.str());
  ofstream bout(strStatFile.str());

  // start rendering screen
  renderer.start();

  // prime the paths?
  renderer.generationPaths = generateSmoothPaths(extraCfg.simNumPathsPerGen, NUM_SEGMENTS_PER_PATH, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);

  // Create a population with this configuration
  *logger.info() << "Creating initial population ..." << endl;
  MyPopulation* pop = new MyPopulation(cfg, adfNs);
  pop->create();
  *logger.info() << "Ok." << endl;
  pop->createGenerationReport(1, 0, fout, bout);

  // This next for statement is the actual genetic programming system
  // which is in essence just repeated reproduction and crossover loop
  // through all the generations ...
  MyPopulation* newPop = NULL;

  for (int gen = 1; gen <= cfg.NumberOfGenerations; gen++)
  {
    // For this generation, build a smooth path goal
    renderer.generationPaths = generateSmoothPaths(extraCfg.simNumPathsPerGen, NUM_SEGMENTS_PER_PATH, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);

    // Create a new generation from the old one by applying the genetic operators
    if (!cfg.SteadyState)
      newPop = new MyPopulation(cfg, adfNs);
    pop->generate(*newPop);

    // TODO fix this pattern to use a dynamic logger
    printEval = true;
    pop->NthMyGP(pop->bestOfPopulation)->evaluate();
    pop->endOfEvaluation();
    printEval = false;

    // Delete the old generation and make the new the old one
    if (!cfg.SteadyState)
    {
      delete pop;
      pop = newPop;
    }

    // Create a report of this generation and how well it is doing
    if (nanDetector > 0) {
      *logger.warn() << "NanDetector count: " << nanDetector << endl;
    }
    pop->createGenerationReport(0, gen, fout, bout);

    // sleep a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }

  // go ahead and dump out the best of the best
  // ofstream bestGP("best.dat");
  // pop->NthMyGP(pop->bestOfPopulation)->save(bestGP);

  // wait for window close
  *logger.info() << "Close window to exit." << endl;
  while (renderer.isRunning()) {
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }
}

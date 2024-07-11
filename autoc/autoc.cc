
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
#include "minisim.h"
#include "pathgen.h"

using namespace std;

// Define configuration parameters and the neccessary array to
// read/write the configuration to a file.  If you need more
// variables, just add them below and insert an entry in the
// configArray.
GPVariables cfg;
ExtraConfig extraCfg;
struct GPConfigVarInformation configArray[]=
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
  {"", DATAINT, NULL}
};


// Define function and terminal identifiers
enum Operators {ADD=0, NEG, MUL, INV,
                IF, EQ, GT, 
                SIN, COS,
                GETDPHI, GETDTHETA, GETDS, GETMX, GETMY, GETMZ, GETOX, GETOY, GETOZ, GETOW, GETVEL,
                PITCH, ROLL, THROTTLE, 
                PI, ZERO, ONE, TWO, PROGN, _END};
const int OPERATORS_NR_ITEM=_END;


// Define class identifiers
const int MyGeneID=GPUserID;
const int MyGPID=GPUserID+1;
const int MyPopulationID=GPUserID+2;

std::atomic_bool printEval = false; // verbose (used for rendering best of population)
std::ofstream fout;
Renderer renderer(extraCfg);

class MyGP;

// Inherit the three GP classes GPGene, GP and GPPopulation
class MyGene : public GPGene
{
public:
  // Duplication (mandatory)
  MyGene (const MyGene& gpo) : GPGene (gpo) { }
  virtual GPObject& duplicate () { return *(new MyGene(*this)); }

  // Creation of own class objects (mandatory)
  MyGene (GPNode& gpo) : GPGene (gpo) {}
  virtual GPGene* createChild (GPNode& gpo) {
    return new MyGene (gpo); }

  // Tree evaluation (not mandatory, but somehow the trees must be
  // parsed to evaluate the fitness)
  double evaluate (std::vector<Path> &path, MyGP &gp, double arg);

  // Load and save (not mandatory)
  MyGene () {}
  virtual int isA () { return MyGeneID; }
  virtual GPObject* createObject() { return new MyGene; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access children (not mandatory)
  MyGene* NthMyChild (int n) {
    return (MyGene*) GPContainer::Nth (n); }
};



class MyGP : public GP 
{
public:
  // Duplication (mandatory)
  MyGP (MyGP& gpo) : GP (gpo) { }
  virtual GPObject& duplicate () { return *(new MyGP(*this)); }

  // Creation of own class objects (mandatory)
  MyGP (int genes) : GP (genes) {}
  virtual GPGene* createGene (GPNode& gpo) {
    return new MyGene (gpo); }

  // Tree evaluation (mandatory)
  virtual void evaluate ();

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Load and save (not mandatory)
  MyGP () {}
  virtual int isA () { return MyGPID; }
  virtual GPObject* createObject() { return new MyGP; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  // Access trees (not mandatory)
  MyGene* NthMyGene (int n) {
    return (MyGene*) GPContainer::Nth (n); }

  // async evaluator
  void evalTask(int gpIndex);

  Aircraft aircraft = Aircraft(0, Eigen::Quaterniond::Identity(), Eigen::Vector3d(0, 0, 0), 0, 0, 0);
  long pathIndex = 0; // current entry on path
};

// Create a Boost.Asio io_context
boost::asio::io_context ioContext;
std::unique_ptr<boost::asio::thread_pool> threadPool = std::make_unique<boost::asio::thread_pool>(4);
std::atomic_ulong taskId(0);
std::atomic_ulong nanDetector(0);
std::vector<MyGP*> tasks = std::vector<MyGP*>();

class MyPopulation : public GPPopulation
{
public:
  // Constructor (mandatory)
  MyPopulation (GPVariables& GPVar_, GPAdfNodeSet& adfNs_) : 
    GPPopulation (GPVar_, adfNs_) {}

  // Duplication (mandatory)
  MyPopulation (MyPopulation& gpo) : GPPopulation(gpo) {}
  virtual GPObject& duplicate () { return *(new MyPopulation(*this)); }

  // Creation of own class objects (mandatory)
  virtual GP* createGP (int numOfTrees) { return new MyGP (numOfTrees); }

  // Load and save (not mandatory)
  MyPopulation () {}
  virtual int isA () { return MyPopulationID; }
  virtual GPObject* createObject() { return new MyPopulation; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);

  virtual void endOfEvaluation () {
    // dispatch all the GPs now (TODO this may still work inline with evaluate)
    for (auto& task : tasks) {
      unsigned long i = taskId++;
      boost::asio::post(*threadPool, [task, i]() {
        task->evalTask(i);
      });
    }

    // wait for all tasks to finish
    threadPool->join();
    tasks.clear();
    
    // TODO how to clean this nasty reload of a thread pool to reactivate after join()
    threadPool  = std::make_unique<boost::asio::thread_pool>(4); // TODO parameter number of cores
  }

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access genetic programs (not mandatory)
  MyGP* NthMyGP (int n) {
    return (MyGP*) GPContainer::Nth (n); }
};

int getIndex(std::vector<Path> &path, MyGP &gp, double arg) {
  if (isnan(arg)) {
    return gp.pathIndex;
  }

  // between now and -10 steps to check, can't go lower than the beginning index
  // TODO this checks path next, not the actual simulation steps...
  // XXX for now, this allows forecasting the future path
  int idx = std::clamp((int) arg, -10, 10) + gp.pathIndex;
  idx = std::clamp(idx, 0, (int) path.size()-1);
  return idx;
}

// This function evaluates the fitness of a genetic tree.  We have the
// freedom to define this function in any way we like.  
double MyGene::evaluate (std::vector<Path> &path, MyGP &run, double arg)
{
  double returnValue = 0.0;

  switch (node->value ())
    {
      case ADD: returnValue=NthMyChild(0)->evaluate (path, run, arg)+NthMyChild(1)->evaluate (path, run, arg); break;
      case NEG: returnValue=-NthMyChild(0)->evaluate (path, run, arg); break;
      case MUL: returnValue=NthMyChild(0)->evaluate (path, run, arg)*NthMyChild(1)->evaluate (path, run, arg); break;
      case INV: {
                double div = NthMyChild(0)->evaluate (path, run, arg);
                returnValue = (div == 0) ? 0 : 1 / div;
                break;
      }
      case IF: returnValue = NthMyChild(0)->evaluate (path, run, arg) ? NthMyChild(1)->evaluate (path, run, arg) : NthMyChild(2)->evaluate (path, run, arg); break;
      case EQ: returnValue = NthMyChild(0)->evaluate (path, run, arg) == NthMyChild(1)->evaluate (path, run, arg); break;
      case GT: returnValue = NthMyChild(0)->evaluate (path, run, arg) > NthMyChild(1)->evaluate (path, run, arg); break;
      case PITCH: returnValue = run.aircraft.setPitchCommand(NthMyChild(0)->evaluate (path, run, arg)); break;
      case ROLL: returnValue = run.aircraft.setRollCommand(NthMyChild(0)->evaluate (path, run, arg)); break;
      case THROTTLE: returnValue = run.aircraft.setThrottleCommand(NthMyChild(0)->evaluate (path, run, arg)); break;
      case SIN: returnValue = sin(NthMyChild(0)->evaluate (path, run, arg)); break;
      case COS: returnValue = cos(NthMyChild(0)->evaluate (path, run, arg)); break;
      case PI: returnValue = M_PI; break;
      case ZERO: returnValue = 0; break;
      case ONE: returnValue = 1; break;
      case TWO: returnValue = 2; break;
      case PROGN: {
                  NthMyChild(0)->evaluate (path, run, arg);
                  returnValue = NthMyChild(1)->evaluate (path, run, arg);
                  break;
      }
      case GETMX: returnValue = run.aircraft.position[0]; break;
      case GETMY: returnValue = run.aircraft.position[1]; break;
      case GETMZ: returnValue = run.aircraft.position[2]; break;
      case GETOX: returnValue = run.aircraft.aircraft_orientation.x(); break;
      case GETOY: returnValue = run.aircraft.aircraft_orientation.y(); break;
      case GETOZ: returnValue = run.aircraft.aircraft_orientation.z(); break;
      case GETOW: returnValue = run.aircraft.aircraft_orientation.w(); break;
      case GETVEL: returnValue = run.aircraft.dRelVel; break;

      case GETDPHI: // compute roll goal from current to target
                  {
                    // Calculate the vector from craft to target in world frame
                    int idx = getIndex(path, run, NthMyChild(0)->evaluate (path, run, arg));
                    Eigen::Vector3d craftToTarget = path.at(idx).start - run.aircraft.position;
                    craftToTarget.normalize();

                    // Transform the craft-to-target vector to body frame
                    Eigen::Vector3d craftToTargetBody = run.aircraft.aircraft_orientation.inverse() * craftToTarget;

                    // Project the craft-to-target vector onto the body YZ plane
                    Eigen::Vector3d projectedVector(0, craftToTargetBody.y(), craftToTargetBody.z());
                    projectedVector.normalize();

                    // Calculate the angle between the projected vector and the body Z-axis
                    returnValue = std::atan2(-projectedVector.y(), projectedVector.z());
                    break;
                  }

      case GETDTHETA: // compute pitch goal from current to target
                  {
                    // Calculate the vector from craft to target in world frame
                    int idx = getIndex(path, run, NthMyChild(0)->evaluate (path, run, arg));
                    Eigen::Vector3d craftToTarget = path.at(idx).start - run.aircraft.position;
                    craftToTarget.normalize();

                    // Transform the craft-to-target vector to body frame
                    Eigen::Vector3d craftToTargetBody = run.aircraft.aircraft_orientation.inverse() * craftToTarget;

                    // Calculate pitch angle
                    returnValue = std::atan2(-craftToTargetBody.x(), 
                                                  std::sqrt(craftToTargetBody.y() * craftToTargetBody.y() + 
                                                            craftToTargetBody.z() * craftToTargetBody.z()));
                    break;
                  }

      case GETDS: // get distance to the next point
                  {
                    int idx = getIndex(path, run, NthMyChild(0)->evaluate (path, run, arg));
                    returnValue = (path.at(idx).start - run.aircraft.position).norm();
                    break;
                  }
      
      default: 
        GPExitSystem ("MyGene::evaluate", "Undefined node value");
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

// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evalTask(int gpIndex)
{
  std::vector<double> fitnesses;

  for (auto &path : renderer.generationPaths) {    
    double localFitness = 0;
    std::vector<Eigen::Vector3d> planPath = std::vector<Eigen::Vector3d>();
    std::vector<Eigen::Vector3d> actualPath = std::vector<Eigen::Vector3d>();

    // deal with pre.path on the initial eval..
    if (path.size() == 0) {
      fitnesses.push_back(1000001);
      continue;
    }
    // north (+x), 5m/s at 10m
    Eigen::Quaterniond aircraft_orientation =
        Eigen::AngleAxisd(SIM_INITIAL_HEADING, Eigen::Vector3d::UnitZ()) *
        Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitY()) *
        Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitX());

    Eigen::Vector3d initialPosition = Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE);
    aircraft = Aircraft(SIM_INITIAL_VELOCITY, aircraft_orientation, initialPosition, 0.0, 0.0, 0.0);

    aircraft.setPitchCommand(0.0);
    aircraft.setRollCommand(0.0);
    aircraft.setThrottleCommand(SIM_INITIAL_THROTTLE);
    
    // iterate the simulator
    double duration = 0.0; // how long have we been running
    pathIndex = 0; // where are we on the path?
    bool printHeader = true;

    // as long as we are within the time limit and have not reached the end of the path
    bool hasCrashed = false;

    // error accumulators
    double distance_error_sum = 0;
    double angle_error_sum = 0;
    double control_smoothness_sum = 0;
    int simulation_steps = 0;

    // initial states
    double roll_prev = aircraft.getRollCommand();
    double pitch_prev = aircraft.getPitchCommand();
    double throttle_prev = aircraft.getThrottleCommand();
    
    while (duration < SIM_TOTAL_TIME && pathIndex < path.size() && !hasCrashed) {

      // walk path looking for next item around TIME_STEP seconds later
      double minDistance = path.at(pathIndex).distanceFromStart + (SIM_TIME_STEP * aircraft.dRelVel);
      int newPathIndex = pathIndex;
      while (newPathIndex < path.size() && (path.at(newPathIndex).distanceFromStart < minDistance)) {
        newPathIndex++;
      }
      // are we off the end?
      if (newPathIndex >= path.size()-1) {
        break;
      }

      // ok, how far is this point from the last point?
      double distance = path.at(newPathIndex).distanceFromStart - path.at(pathIndex).distanceFromStart;
      // so this is the real dT
      double dT = distance / SIM_INITIAL_VELOCITY; // assume a constant speed trajectory

      // advance the simulator
      duration += dT;
      pathIndex = newPathIndex;

      // GP determine control input
      NthMyGene (0)->evaluate (path, *this, 0);

      // and advances the aircract
      aircraft.advanceState(dT);

      // ---------------------- how did we do?
      
      // Compute the distance between the aircraft and the goal
      double distanceFromGoal = (path.at(pathIndex).start - aircraft.position).norm();

      // Compute vector from target to current aircraft position
      Eigen::Vector3d target_direction = (path.at(pathIndex+1).start - path.at(pathIndex).start);
      Eigen::Vector3d aircraft_to_target = (path.at(pathIndex).start - aircraft.position);
      double dot_product = target_direction.dot(aircraft_to_target);
      double angle_rad = std::acos(std::clamp(dot_product / (target_direction.norm() * aircraft_to_target.norm()), -1.0, 1.0));

      // control smoothness -- internal vector distance
      double smoothness = pow(roll_prev - aircraft.getRollCommand(), 2.0);
      smoothness += pow(pitch_prev - aircraft.getPitchCommand(), 2.0);
      smoothness += pow(throttle_prev - aircraft.getThrottleCommand(), 2.0);
      smoothness = sqrt(smoothness);

      // ready for next cycle
      roll_prev = aircraft.getRollCommand();
      pitch_prev = aircraft.getPitchCommand();
      throttle_prev = aircraft.getThrottleCommand();

      // but have we crashed outside the sphere?
      double distanceFromOrigin = (aircraft.position - Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE)).norm();
      if (aircraft.position[2] > SIM_MIN_ELEVATION || distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) {
        hasCrashed = true;
      }

      if (printEval) {
        if (printHeader) {
          fout << "    Time Idx       dT  totDist   pathX    pathY    pathZ        X        Y        Z       dW       dX       dY       dZ   relVel     roll    pitch    power    distP   angleP controlP\n";
          printHeader = false;
        }

        char outbuf[1000]; // XXX use c++20
        sprintf(outbuf, "% 8.2f %3ld % 8.2f % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n", 
          duration, pathIndex, dT, 
          path.at(pathIndex).distanceFromStart,
          path.at(pathIndex).start[0],
          path.at(pathIndex).start[1],
          path.at(pathIndex).start[2],
          aircraft.position[0],
          aircraft.position[1],
          aircraft.position[2],
          aircraft.aircraft_orientation.w(),
          aircraft.aircraft_orientation.x(),
          aircraft.aircraft_orientation.y(),
          aircraft.aircraft_orientation.z(),
          aircraft.dRelVel,
          aircraft.getRollCommand(),
          aircraft.getPitchCommand(),
          aircraft.getThrottleCommand(),
          distanceFromGoal,
          angle_rad,
          smoothness
        );
          fout << outbuf;

        // now update points for Renderer
        planPath.push_back(path.at(pathIndex).start); // XXX push the full path
        actualPath.push_back(aircraft.position);
      }

      distance_error_sum += distanceFromGoal;
      angle_error_sum += angle_rad;
      control_smoothness_sum += smoothness;
      simulation_steps++;

    }

    // tally up the normlized fitness
    double normalized_distance_error = (distance_error_sum / simulation_steps) * FITNESS_DISTANCE_WEIGHT;
    double normalized_velocity_align = (angle_error_sum / simulation_steps) * FITNESS_ALIGNMENT_WEIGHT;
    double normalized_control_smoothness = (control_smoothness_sum / simulation_steps) * FITNESS_CONTROL_WEIGHT;

    localFitness = normalized_distance_error + normalized_velocity_align + normalized_control_smoothness;

    if (isnan(localFitness)) {
      nanDetector++;
    }

    if (hasCrashed) {
      localFitness += 1000; // TODO
    }
          
    if (printEval) {
      // now update actual list of lists
      std::vector<Path> actualElementList;
      for (auto& actualElement : actualPath) {
        actualElementList.push_back({actualElement, 0});
      }

      // now update plan list of lists
      std::vector<Path> planElementList;
      for (auto& planElement : planPath) {
        planElementList.push_back({planElement, 0});
      }
      renderer.addPathElementList(planElementList, actualElementList);

      planPath.clear();
      actualPath.clear();
    }

    fitnesses.push_back(localFitness);
  }

  // median value is the fitness for this bunch
  std::sort(fitnesses.begin(), fitnesses.end());
  stdFitness = fitnesses.at(fitnesses.size() / 2);

  if (printEval) {
    renderer.update();
  }
}


// Create function and terminal set
void createNodeSet (GPAdfNodeSet& adfNs)
{
  // Reserve space for the node sets
  adfNs.reserveSpace (1);
  
  // Now define the function and terminal set for each ADF and place
  // function/terminal sets into overall ADF container
  GPNodeSet& ns=*new GPNodeSet (OPERATORS_NR_ITEM);

  adfNs.put (0, ns);
  
  // Define functions/terminals and place them into the appropriate
  // sets.  Terminals take two arguments, functions three (the third
  // parameter is the number of arguments the function has)
  ns.putNode (*new GPNode (ADD, "ADD", 2));
  ns.putNode (*new GPNode (NEG, "NEG", 1));
  ns.putNode (*new GPNode (MUL, "MUL", 2));
  ns.putNode (*new GPNode (INV, "INV", 1));
  ns.putNode (*new GPNode (IF, "IF", 3));
  ns.putNode (*new GPNode (EQ, "EQ", 2));
  ns.putNode (*new GPNode (GT, "GT", 2));
  ns.putNode (*new GPNode (PITCH, "PITCH", 1));
  ns.putNode (*new GPNode (ROLL, "ROLL", 1));
  ns.putNode (*new GPNode (THROTTLE, "THROTTLE", 1));
  ns.putNode (*new GPNode (SIN, "SIN", 1));
  ns.putNode (*new GPNode (COS, "COS", 1));
  ns.putNode (*new GPNode (PI, "PI"));
  ns.putNode (*new GPNode (ZERO, "0"));
  ns.putNode (*new GPNode (ONE, "1"));
  ns.putNode (*new GPNode (TWO, "2"));
  ns.putNode (*new GPNode (PROGN, "PROGN", 2));
  ns.putNode (*new GPNode (GETDPHI, "GETDPHI", 1));
  ns.putNode (*new GPNode (GETDTHETA, "GETDTHETA", 1));
  ns.putNode (*new GPNode (GETDS, "GETDS", 1));
  ns.putNode (*new GPNode (GETMX, "GETMX"));
  ns.putNode (*new GPNode (GETMY, "GETMY"));
  ns.putNode (*new GPNode (GETMZ, "GETMZ"));
  ns.putNode (*new GPNode (GETOX, "GETOX"));
  ns.putNode (*new GPNode (GETOY, "GETOY"));
  ns.putNode (*new GPNode (GETOZ, "GETOZ"));
  ns.putNode (*new GPNode (GETOW, "GETOW"));
  ns.putNode (*new GPNode (GETVEL, "GETVEL"));
}



void newHandler ()
{
  cerr << "\nFatal error: Out of memory." << endl;
  exit (1);
}



int main ()
{
  // Set up a new-handler, because we might need a lot of memory, and
  // we don't know it's there.
  set_new_handler (newHandler);

  // Init GP system.
  GPInit (1, -1);
  
  // Read configuration file.
  GPConfiguration config (cout, "autoc.ini", configArray);
  renderer.extraCfg = extraCfg;
  
  // Print the configuration
  cout << cfg << endl;
  cout << "SimNumPathsPerGen: " << extraCfg.simNumPathsPerGen << endl;
  
  // Create the adf function/terminal set and print it out.
  GPAdfNodeSet adfNs;
  createNodeSet (adfNs);
  cout << adfNs << endl; 
  
  // Open the main output file for the data and statistics file.
  // First set up names for data file.  Remember we should delete the
  // string from the stream, well just a few bytes
  ostringstream strOutFile, strStatFile;
  strOutFile  << "data.dat" << ends;
  strStatFile << "data.stc" << ends;
  fout.open(strOutFile.str());
  ofstream bout (strStatFile.str());

  // start rendering screen
  renderer.start();
  
  // prime the paths?
  renderer.generationPaths = generateSmoothPaths(extraCfg.simNumPathsPerGen, NUM_SEGMENTS_PER_PATH, SIM_PATH_BOUNDS);
  
  // Create a population with this configuration
  cout << "Creating initial population ..." << endl;
  MyPopulation* pop=new MyPopulation (cfg, adfNs);
  pop->create ();
  cout << "Ok." << endl;
  pop->createGenerationReport (1, 0, fout, bout);
  
  // This next for statement is the actual genetic programming system
  // which is in essence just repeated reproduction and crossover loop
  // through all the generations ...
  MyPopulation* newPop=NULL;
  
  for (int gen=1; gen<=cfg.NumberOfGenerations; gen++)
  {
    // For this generation, build a smooth path goal
    renderer.generationPaths = generateSmoothPaths(extraCfg.simNumPathsPerGen, NUM_SEGMENTS_PER_PATH, SIM_PATH_BOUNDS);
    
    // Create a new generation from the old one by applying the genetic operators
    if (!cfg.SteadyState)
      newPop=new MyPopulation (cfg, adfNs);
    pop->generate (*newPop);
    
    // TODO fix this pattern to use a dynamic logger
    printEval = true;
    pop->NthMyGP(pop->bestOfPopulation)->evaluate();
    pop->endOfEvaluation();
    printEval = false;

    // Delete the old generation and make the new the old one
    if (!cfg.SteadyState)
    {
      delete pop;
      pop=newPop;
    }

    // Create a report of this generation and how well it is doing
    if (nanDetector > 0) {
      cout << "NanDetector count: " << nanDetector << endl;
    }
    pop->createGenerationReport (0, gen, fout, bout);

    // sleep a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }

  // wait for window close
  cout << "Close window to exit." << endl;
  while (renderer.isRunning()) {
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }
}

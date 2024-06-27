
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
  {"", DATAINT, NULL}
};


// Define function and terminal identifiers
enum Operators {ADD=0, NEG, MUL, INV,
                IF, EQ, GT, 
                SIN, COS,
                GETDPHI, GETDTHETA, GETDS, GETMX, GETMY, GETMZ, GETOX, GETOY, GETOZ, GETOW,
                PITCH, ROLL, THROTTLE, 
                PI, ZERO, ONE, TWO, PROGN, _END};
const int OPERATORS_NR_ITEM=_END;


// Define class identifiers
const int MyGeneID=GPUserID;
const int MyGPID=GPUserID+1;
const int MyPopulationID=GPUserID+2;

// Create a Boost.Asio io_context
boost::asio::io_context ioContext;
std::unique_ptr<boost::asio::thread_pool> threadPool = std::make_unique<boost::asio::thread_pool>(4);
std::atomic_ulong taskId(0);

// the path(s) a population will attempt
std::vector<Path> path = std::vector<Path>();
bool printEval = false; // verbose (used for rendering best of population)
std::ofstream fout;
Renderer renderer = Renderer();

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
  double evaluate (SimRun *run, double arg);

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
  void evalTask(unsigned long taskIndex);
};



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
    threadPool->join();
    
    // TODO how to clean this nasty reload of a thread pool to reactivate after join()
    threadPool  = std::make_unique<boost::asio::thread_pool>(4); // TODO parameter number of cores
  }

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access genetic programs (not mandatory)
  MyGP* NthMyGP (int n) {
    return (MyGP*) GPContainer::Nth (n); }
};

int getIndex(SimRun *run, double arg) {
  if (isnan(arg)) {
    return run->pathIndex;
  }

  int idx = min(max((int) arg, -1), 5) + run->pathIndex;

  return idx < path.size()-1 ? idx : path.size()-1;
}

// This function evaluates the fitness of a genetic tree.  We have the
// freedom to define this function in any way we like.  
double MyGene::evaluate (SimRun *run, double arg)
{
  double returnValue = 0.0;

  switch (node->value ())
    {
      case ADD: returnValue=NthMyChild(0)->evaluate (run, arg)+NthMyChild(1)->evaluate (run, arg); break;
      case NEG: returnValue=-NthMyChild(0)->evaluate (run, arg); break;
      case MUL: returnValue=NthMyChild(0)->evaluate (run, arg)*NthMyChild(1)->evaluate (run, arg); break;
      case INV: {
                double div = NthMyChild(0)->evaluate (run, arg);
                returnValue = (div == 0) ? 0 : 1 / div;
                break;
      }
      case IF: returnValue = NthMyChild(0)->evaluate (run, arg) ? NthMyChild(1)->evaluate (run, arg) : NthMyChild(2)->evaluate (run, arg); break;
      case EQ: returnValue = NthMyChild(0)->evaluate (run, arg) == NthMyChild(1)->evaluate (run, arg); break;
      case GT: returnValue = NthMyChild(0)->evaluate (run, arg) > NthMyChild(1)->evaluate (run, arg); break;
      case PITCH: returnValue = run->aircraft->setPitchCommand(NthMyChild(0)->evaluate (run, arg)); break;
      case ROLL: returnValue = run->aircraft->setRollCommand(NthMyChild(0)->evaluate (run, arg)); break;
      case THROTTLE: returnValue = run->aircraft->setThrottleCommand(NthMyChild(0)->evaluate (run, arg)); break;
      case SIN: returnValue = sin(NthMyChild(0)->evaluate (run, arg)); break;
      case COS: returnValue = cos(NthMyChild(0)->evaluate (run, arg)); break;
      case PI: returnValue = M_PI; break;
      case ZERO: returnValue = 0; break;
      case ONE: returnValue = 1; break;
      case TWO: returnValue = 2; break;
      case PROGN: {
                  NthMyChild(0)->evaluate (run, arg);
                  returnValue = NthMyChild(1)->evaluate (run, arg);
                  break;
      }
      case GETMX: returnValue = run->aircraft->getState()->position[0]; break;
      case GETMY: returnValue = run->aircraft->getState()->position[1]; break;
      case GETMZ: returnValue = run->aircraft->getState()->position[2]; break;
      case GETOX: returnValue = run->aircraft->getState()->aircraft_orientation.x(); break;
      case GETOY: returnValue = run->aircraft->getState()->aircraft_orientation.y(); break;
      case GETOZ: returnValue = run->aircraft->getState()->aircraft_orientation.z(); break;
      case GETOW: returnValue = run->aircraft->getState()->aircraft_orientation.w(); break;

      case GETDPHI: // compute roll goal from current to target
                  {
                    // Aircraft's initial position and orientation (in world frame)
                    Eigen::Vector3d aircraft_position = run->aircraft->getState()->position;

                    // Target's position (in world frame)
                    int idx = getIndex(run, NthMyChild(0)->evaluate (run, arg));
                    Eigen::Vector3d target_position(path.at(idx).start[0], path.at(idx).start[1], path.at(idx).start[2]);

                    // Compute the vector from the aircraft to the target
                    Eigen::Vector3d target_vector = target_position - aircraft_position;

                    // Normalize the target vector to get the direction
                    Eigen::Vector3d target_direction = target_vector.normalized();

                    // Aircraft's current z-axis vector in world frame
                    Eigen::Vector3d aircraft_z_axis(0, 0, 1); // Assuming unit vector along the z-axis
                    aircraft_z_axis = run->aircraft->getState()->aircraft_orientation * aircraft_z_axis;

                    // Compute the quaternion representing the rotation from aircraft z-axis to target direction
                    Eigen::Quaterniond quat_current_to_target = Eigen::Quaterniond::FromTwoVectors(aircraft_z_axis, target_direction);

                    // Extract the roll delta from the quaternion
                    Eigen::Vector3d euler_angles_delta = quat_current_to_target.toRotationMatrix().eulerAngles(2, 1, 0); // ZYX order: yaw, pitch, roll
                    double returnValue = euler_angles_delta(2);
                    break;
                  }

      case GETDTHETA: // compute pitch goal from current to target
                  {
                    // Aircraft's initial position and orientation (in world frame)
                    Eigen::Vector3d aircraft_position = run->aircraft->getState()->position;

                    // Target's position (in world frame)
                    int idx = getIndex(run, NthMyChild(0)->evaluate (run, arg));
                    Eigen::Vector3d target_position(path.at(idx).start[0], path.at(idx).start[1], path.at(idx).start[2]);

                    // Compute the vector from the aircraft to the target
                    Eigen::Vector3d target_vector = target_position - aircraft_position;

                    // Normalize the target vector to get the direction
                    Eigen::Vector3d target_direction = target_vector.normalized();

                    // Aircraft's current x-axis velocity vector in world frame
                    Eigen::Vector3d aircraft_velocity(1, 0, 0); // Assuming unit vector along the x-axis
                    aircraft_velocity = run->aircraft->getState()->aircraft_orientation * aircraft_velocity;

                    // Compute the quaternion representing the rotation from aircraft velocity to target direction
                    Eigen::Quaterniond quat_current_to_target = Eigen::Quaterniond::FromTwoVectors(aircraft_velocity, target_direction);

                    // Extract the pitch delta from the quaternion
                    Eigen::Vector3d euler_angles_delta = quat_current_to_target.toRotationMatrix().eulerAngles(2, 1, 0); // ZYX order: yaw, pitch, roll
                    double returnValue = euler_angles_delta(1);
                    break;
                  }

      case GETDS: // get distance to the next point
                  {
                    int idx = getIndex(run, NthMyChild(0)->evaluate (run, arg));
                    returnValue = (path.at(idx).start - run->aircraft->getState()->position).norm();
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
  unsigned long i = taskId++;

  boost::asio::post(*threadPool, [this, i]() {
    this->evalTask(i);
  });
}

// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evalTask(unsigned long taskIndex)
{
  // deal with pre.path on the initial eval..
  if (path.size() == 0) {
    stdFitness = 1000001;
    return;
  } else {
    stdFitness = 0.0;
  }

  SimRun *run = new SimRun();

  // north (+x), 5m/s at 10m
  Eigen::Quaterniond aircraft_orientation =
      Eigen::AngleAxisd(SIM_INITIAL_HEADING, Eigen::Vector3d::UnitZ()) *
      Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitY()) *
      Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitX());

  Eigen::Vector3d initialPosition = Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE);
  AircraftState state = {SIM_INITIAL_VELOCITY, aircraft_orientation, initialPosition, 0.0, 0.0, 0.0};
  run->aircraft->setState(&state);
  run->aircraft->setPitchCommand(0.0);
  run->aircraft->setRollCommand(0.0);
  run->aircraft->setThrottleCommand(SIM_INITIAL_THROTTLE);
  
  // char output[200];
  // aircraft->toString(output);
  // cout << output << endl;

  // iterate the simulator
  double duration = 0.0; // how long have we been running
  run->pathIndex = 0; // where are we on the path?
  bool printHeader = true;

  std::vector<Eigen::Vector3d> planPath = std::vector<Eigen::Vector3d>();
  std::vector<Eigen::Vector3d> actualPath = std::vector<Eigen::Vector3d>();

  planPath.push_back(path.at(run->pathIndex).start); // XXX push the full path
  actualPath.push_back(run->aircraft->getState()->position);

  // as long as we are within the time limit and have not reached the end of the path
  bool hasCrashed = false; 
  while (duration < SIM_TOTAL_TIME && run->pathIndex < path.size() && !hasCrashed) {

    // walk path looking for next item around TIME_STEP seconds later
    double minDistance = path.at(run->pathIndex).distanceFromStart + (SIM_TIME_STEP * SIM_INITIAL_VELOCITY);
    int newPathIndex = run->pathIndex;
    while (newPathIndex < path.size() && (path.at(newPathIndex).distanceFromStart < minDistance)) {
      newPathIndex++;
    }
    // are we off the end?
    if (newPathIndex >= path.size()) {
      break;
    }

    // ok, how far is this point from the last point?
    double distance = path.at(newPathIndex).distanceFromStart - path.at(run->pathIndex).distanceFromStart;
    // so this is the real dT
    double dT = distance / SIM_INITIAL_VELOCITY;

    // advance the simulator
    duration += dT;
    run->pathIndex = newPathIndex;

    // GP determine control input
    NthMyGene (0)->evaluate (run, 0);

    // and advances the aircraft
    run->aircraft->advanceState(dT);

    // how did we do?
    // compute the delta vector from the aircraft vector to the goal vector based on
    // alignment (vector direction diff to velocity vector) and distance

    // Aircraft's current x-axis normalized velocity vector in world frame
    Eigen::Vector3d aircraft_velocity(1, 0, 0); // Assuming unit vector along the x-axis
    aircraft_velocity = run->aircraft->getState()->aircraft_orientation * aircraft_velocity;
    Eigen::Vector3d aircraft_velocity_normalized = aircraft_velocity.normalized();

    // Target's direction normalized vector (in world frame)
    Eigen::Vector3d target_direction(path.at(run->pathIndex).start - path.at(run->pathIndex-1).start);
    target_direction.normalize();

    // Compute the dot product of the two normalized vectors
    double dot_product = target_direction.dot(aircraft_velocity_normalized);

    // Clamp the dot product to avoid numerical issues with acos
    dot_product = std::clamp(dot_product, -1.0, 1.0);
    double angle_rad = std::acos(dot_product);

    // Compute the distance between the aircraft and the goal
    double distanceFromGoal = (path.at(run->pathIndex).start - run->aircraft->getState()->position).norm();

    // Compute the fitness value (distance * some power function) + (angle penalty * some power function)
    double anglePenalty = pow(SIM_ANGLE_SCALE_FACTOR * fabs(angle_rad), SIM_ANGLE_PENALTY_FACTOR);
    double distancePenalty = pow(distanceFromGoal, SIM_DISTANCE_PENALTY_FACTOR);
    double fitness = distancePenalty + anglePenalty;

    // add in distance component
    if (!isnan(fitness)) {
      stdFitness += fitness;
    } else {
      stdFitness += 1000000;
    }

    // but have we crashed outside the sphere?
    double distanceFromOrigin = (run->aircraft->getState()->position - Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE)).norm();
    if (run->aircraft->getState()->position[2] > SIM_MIN_ELEVATION || distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) {
      // ok we are outside the bounds -- penalize at whatever rate of error we had so far
      double timeRemaining = max(0.0, SIM_TOTAL_TIME - duration);
      double preCrashFitness = stdFitness / duration;
      double projectedFitness = preCrashFitness * timeRemaining;
      stdFitness += pow(projectedFitness, SIM_CRASH_FITNESS_PENALTY_FACTOR);
      hasCrashed = true;
    }

    if (printEval) {
      if (printHeader) {
        fout << "    Time Idx       dT  totDist   pathX    pathY    pathZ        X        Y        Z       dW       dX       dY       dZ   relVel       dG     roll    pitch    power  fitness   angleP    distP\n";
        printHeader = false;
      }

      char outbuf[1000]; // XXX use c++20
      sprintf(outbuf, "% 8.2f %3ld % 8.2f % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n", 
        duration, run->pathIndex, dT, 
        path.at(run->pathIndex).distanceFromStart,
        path.at(run->pathIndex).start[0],
        path.at(run->pathIndex).start[1],
        path.at(run->pathIndex).start[2],
        run->aircraft->getState()->position[0],
        run->aircraft->getState()->position[1],
        run->aircraft->getState()->position[2],
        run->aircraft->getState()->aircraft_orientation.w(),
        run->aircraft->getState()->aircraft_orientation.x(),
        run->aircraft->getState()->aircraft_orientation.y(),
        run->aircraft->getState()->aircraft_orientation.z(),
        run->aircraft->getState()->dRelVel,
        distanceFromGoal,
        run->aircraft->getRollCommand(),
        run->aircraft->getPitchCommand(),
        run->aircraft->getThrottleCommand(),
        stdFitness,
        anglePenalty,
        distancePenalty
        );
        fout << outbuf;

      // now update points for Renderer
      planPath.push_back(path.at(run->pathIndex).start); // XXX push the full path
      actualPath.push_back(run->aircraft->getState()->position);
    }
  }

  if (printEval) {
    renderer.update(planPath, actualPath);
    planPath.clear();
    actualPath.clear();
  }

  delete run;
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
  
  // Print the configuration
  cout << cfg << endl;
  
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
    path = generateSmoothPath(18, SIM_PATH_BOUNDS); // TODO parameterize points
    
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
    pop->createGenerationReport (0, gen, fout, bout);

    // clean out prior fitness case
    path.clear();

    // sleep a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }

  // wait for window close
  cout << "Close window to exit." << endl;
  while (renderer.isRunning()) {
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
  }
}

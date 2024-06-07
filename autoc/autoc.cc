
// autoc.cc

/* -------------------------------------------------------------------
From skeleton/skeleton.cc
------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include <new>
#include <fstream>
#include <sstream>
#include <random>

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
enum Operators {ADD=0, NEG, MUL, INV, IF, EQ, GT, SIN, COS, TAN, GETDR, GETDS, GETMX, GETMY, GETPSI, PITCH, ROLL, THROTTLE, PI, ZERO, ONE, TWO, PROGN, _END};
const int OPERATORS_NR_ITEM=_END;


// Define class identifiers
const int MyGeneID=GPUserID;
const int MyGPID=GPUserID+1;
const int MyPopulationID=GPUserID+2;

std::uniform_real_distribution<double> dist(-M_PI, M_PI);  //(min, max)
//Mersenne Twister: Good quality random number generator
std::mt19937 rng; 

Aircraft *aircraft = new Aircraft(new AircraftState());
std::vector<Path> path = std::vector<Path>();
unsigned long pathIndex = 0; // current entry on path
bool printEval = false;
std::ofstream fout;

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
  double evaluate (double arg);

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

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access genetic programs (not mandatory)
  MyGP* NthMyGP (int n) {
    return (MyGP*) GPContainer::Nth (n); }
};



// This function evaluates the fitness of a genetic tree.  We have the
// freedom to define this function in any way we like.  
double MyGene::evaluate (double arg)
{
  double returnValue = 0.0;

  switch (node->value ())
    {
      case ADD: returnValue=NthMyChild(0)->evaluate (arg)+NthMyChild(1)->evaluate (arg); break;
      case NEG: returnValue=-NthMyChild(0)->evaluate (arg); break;
      case MUL: returnValue=NthMyChild(0)->evaluate (arg)*NthMyChild(1)->evaluate (arg); break;
      case INV: {
                double div = NthMyChild(0)->evaluate (arg);
                returnValue = (div == 0) ? 0 : 1 / div;
                break;
      }
      case IF: returnValue = NthMyChild(0)->evaluate (arg) ? NthMyChild(1)->evaluate (arg) : NthMyChild(2)->evaluate (arg); break;
      case EQ: returnValue = NthMyChild(0)->evaluate (arg) == NthMyChild(1)->evaluate (arg); break;
      case GT: returnValue = NthMyChild(0)->evaluate (arg) > NthMyChild(1)->evaluate (arg); break;
      case PITCH: returnValue = aircraft->setPitchCommand(NthMyChild(0)->evaluate (arg)); break;
      case ROLL: returnValue = aircraft->setRollCommand(NthMyChild(0)->evaluate (arg)); break;
      case THROTTLE: returnValue = aircraft->setThrottleCommand(NthMyChild(0)->evaluate (arg)); break;
      case SIN: returnValue = sin(NthMyChild(0)->evaluate (arg)); break;
      case COS: returnValue = cos(NthMyChild(0)->evaluate (arg)); break;
      case TAN: returnValue = tan(NthMyChild(0)->evaluate (arg)); break;
      case PI: returnValue = M_PI; break;
      case ZERO: returnValue = 0; break;
      case ONE: returnValue = 1; break;
      case TWO: returnValue = 2; break;
      case PROGN: {
                  NthMyChild(0)->evaluate (arg);
                  returnValue = NthMyChild(1)->evaluate (arg);
                  break;
      }
      case GETMX: returnValue = aircraft->getState()->X; break;
      case GETMY: returnValue = aircraft->getState()->Y; break;
      case GETPSI: returnValue = aircraft->getState()->dPsi; break;

      case GETDR: // get left/right angle from my path to the next point
                  {
                    unsigned long idx = max(min(NthMyChild(0)->evaluate (arg) + pathIndex, (double) min(path.size()-1, (pathIndex+5))), (double) pathIndex);
                    double dx = path.at(idx).start.x - aircraft->getState()->X;
                    double dy = path.at(idx).start.y - aircraft->getState()->Y;
                    double angle = atan2(dy, dx);
                    double dPsi = aircraft->getState()->dPsi;
                    returnValue = remainder(angle - dPsi, M_PI * 2);
                    break;
                  }
      case GETDS: // get distance to the next point
                  {
                    unsigned long idx = max(min(NthMyChild(0)->evaluate (arg) + pathIndex, (double) min(path.size()-1, (pathIndex+5))), (double) pathIndex);
                    double dx = path.at(idx).start.x - aircraft->getState()->X;
                    double dy = path.at(idx).start.y - aircraft->getState()->Y;
                    double dz = path.at(idx).start.z - aircraft->getState()->Z;
                    returnValue = sqrt(dx*dx + dy*dy + dz*dz);
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



// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evaluate ()
{
  // deal with pre.path on the initial eval..
  if (path.size() == 0) {
    stdFitness = 1000001;
    return;
  } else {
    stdFitness = 0.0;
  }

  // north, 5m/s at 10m
  AircraftState state = {SIM_INITIAL_VELOCITY, 0.0, 0.0, 0.0, 0.0, 0.0, SIM_INITIAL_ALTITUDE, 0.0, 0.0, 0.0};
  aircraft->setState(&state);
  aircraft->setPitchCommand(0.0);
  aircraft->setRollCommand(0.0);
  aircraft->setThrottleCommand(SIM_INITIAL_THROTTLE);
  // char output[200];
  // aircraft->toString(output);
  // cout << output << endl;

  // iterate the simulator
  double duration = 0.0; // how long have we been running
  pathIndex = 0; // where are we on the path?
  bool printHeader = true;

  // as long as we are within the time limit and have not reached the end of the path
  while (duration < SIM_TOTAL_TIME && pathIndex < path.size()) {

    // walk path looking for next item around 0.5 seconds later
    double minDistance = path.at(pathIndex).distanceFromStart + (0.5 * SIM_INITIAL_VELOCITY);
    int newPathIndex = pathIndex;
    while (newPathIndex < path.size() && (path.at(newPathIndex).distanceFromStart < minDistance)) {
      newPathIndex++;
    }
    // are we off the end?
    if (newPathIndex == path.size()) {
      break;
    }

    // ok, how far is this point from the last point?
    double distance = path.at(newPathIndex).distanceFromStart - path.at(pathIndex).distanceFromStart;
    // so this is the real dT
    double dT = distance / SIM_INITIAL_VELOCITY;

    // advance the simulator
    duration += dT;
    pathIndex = newPathIndex;

    // GP determine control input
    NthMyGene (0)->evaluate (0);

    // and advances the aircraft
    aircraft->advanceState(dT);

    // how did we do?
    // for now how close are we to the goal point?
    double distanceFromGoal = std::sqrt(
          std::pow(aircraft->getState()->X - path.at(pathIndex).start.x, 2) +
          std::pow(aircraft->getState()->Y - path.at(pathIndex).start.y, 2) +
          std::pow(aircraft->getState()->Z - path.at(pathIndex).start.z, 2));

    // add in distance component
    stdFitness += distanceFromGoal;

    // but have we crashed outside the sphere?
    double x = aircraft->getState()->X;
    double y = aircraft->getState()->Y;
    double z = aircraft->getState()->Z;
    if (aircraft->getState()->Z < 0 || std::sqrt(x*x + y*y + z*z) > SIM_PATH_RADIUS_LIMIT) { // XXX parameterize
      // ok we are outside the bounds -- penalize but 
      double timeRemaining = max(0.0, SIM_TOTAL_TIME - duration); // XXX parameterize
      stdFitness += SIM_CRASH_FITNESS_PENALTY + timeRemaining * SIM_INITIAL_VELOCITY; // big-xx our velocity
      break;
    }

    if (printEval) {
      if (printHeader) {
        fout << "    Time Idx       dT    pathX    pathY    pathZ        X        Y        Z     dPsi     dPhi   dTheta   relVel       dG     roll    pitch    power  fitness\n";
        printHeader = false;
      }

      char outbuf[1000]; // XXX use c++20
      sprintf(outbuf, "% 8.2f %3ld % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n", 
        duration, pathIndex, dT, 
        path.at(pathIndex).start.x,
        path.at(pathIndex).start.y,
        path.at(pathIndex).start.z,
        x, y, z,
        aircraft->getState()->dPsi,
        aircraft->getState()->dPhi,
        aircraft->getState()->dTheta,
        aircraft->getState()->dRelVel,
        distanceFromGoal,
        aircraft->getRollCommand(),
        aircraft->getPitchCommand(),
        aircraft->getThrottleCommand(),
        stdFitness);
        fout << outbuf;
    }
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
  ns.putNode (*new GPNode (TAN, "TAN", 1));
  ns.putNode (*new GPNode (PI, "PI"));
  ns.putNode (*new GPNode (ZERO, "0"));
  ns.putNode (*new GPNode (ONE, "1"));
  ns.putNode (*new GPNode (TWO, "2"));
  ns.putNode (*new GPNode (PROGN, "PROGN", 2));
  ns.putNode (*new GPNode (GETDR, "GETDR", 1));
  ns.putNode (*new GPNode (GETDS, "GETDS", 1));
  ns.putNode (*new GPNode (GETMX, "GETMX"));
  ns.putNode (*new GPNode (GETMY, "GETMY"));
  ns.putNode (*new GPNode (GETPSI, "GETPSI"));
  // ns.putNode (*new GPNode (GETMZ, "GETMZ"));
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

  //Initialize with non-deterministic seeds
  rng.seed(std::random_device{}()); 
  
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
  
  // THIS IS FOR ALL GENERATIONS...
  // path = generateSmoothPath(20, SIM_PATH_BOUNDS); // TODO parameterize points

  for (int gen=1; gen<=cfg.NumberOfGenerations; gen++)
    {
      // For this generation, build a smooth path goal
      path = generateSmoothPath(20, SIM_PATH_BOUNDS); // TODO parameterize points
      // char *output = new char[200];
      // for (int i = 0; i < path.size(); i++) {
      //   path.at(i).toString(output);
      //   cout << i << ":" << output << endl;
      // }
      // printf("Path: %lu points, total length %f\n", path.size(), path.at(path.size() - 1).distanceFromStart);

      // Create a new generation from the old one by applying the
      // genetic operators
      if (!cfg.SteadyState)
	      newPop=new MyPopulation (cfg, adfNs);
      pop->generate (*newPop);
      
            // XXX fix this pattern to use a dynamic logger
      printEval = true;
      pop->NthMyGP(pop->bestOfPopulation)->evaluate();
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
    }
}

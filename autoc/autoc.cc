
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
enum Operators {ADD=0, SUB, MUL, DIV, REM, ABS, NEG, MAX, MIN, IF, GT, LT, EQ, ZERO, ONE, LEFT, RIGHT, UP, DOWN, FAST, SLOW, GET};
const int OPERATORS_NR_ITEM=22;


// Define class identifiers
const int MyGeneID=GPUserID;
const int MyGPID=GPUserID+1;
const int MyPopulationID=GPUserID+2;

// temporary buffers
#define INPUT_SIZE 1
class Payload {
public:
  double inputVar[INPUT_SIZE];
  double roll;
  double pitch;
  double throttle;
};

Payload *payload = new Payload();

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
  double evaluate ();

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
double MyGene::evaluate ()
{
  double returnValue = 0.0;

  switch (node->value ())
    {
      case ADD: returnValue=NthMyChild(0)->evaluate ()+NthMyChild(1)->evaluate (); break;
      case SUB: returnValue=NthMyChild(0)->evaluate ()-NthMyChild(1)->evaluate (); break;
      case MUL: returnValue=NthMyChild(0)->evaluate ()*NthMyChild(1)->evaluate (); break;
      case DIV: { double divisor = NthMyChild(1)->evaluate ();
                  divisor == 0 ? 0 : returnValue=NthMyChild(0)->evaluate ()/divisor; 
                  break;
      }
      case REM: { double divisor = NthMyChild(1)->evaluate ();
                  divisor == 0 ? 0 : returnValue=remainder(NthMyChild(0)->evaluate (), divisor); 
                  break;
      }
      case ABS: returnValue=abs (NthMyChild(0)->evaluate ()); break;
      case NEG: returnValue=-NthMyChild(0)->evaluate (); break;
      case MAX: returnValue=max (NthMyChild(0)->evaluate (), NthMyChild(1)->evaluate ()); break;
      case MIN: returnValue=min (NthMyChild(0)->evaluate (), NthMyChild(1)->evaluate ()); break;
      case IF: returnValue=NthMyChild(0)->evaluate () ? NthMyChild(1)->evaluate () : NthMyChild(2)->evaluate (); break;
      case GT: returnValue=NthMyChild(0)->evaluate () > NthMyChild(1)->evaluate (); break;
      case LT: returnValue=NthMyChild(0)->evaluate () < NthMyChild(1)->evaluate (); break;
      case EQ: returnValue=NthMyChild(0)->evaluate () == NthMyChild(1)->evaluate (); break;
      case ZERO: returnValue=0; break;
      case ONE: returnValue=1; break;
      case LEFT: payload->roll = max(-100.0, payload->roll - 1); break;
      case RIGHT: payload->roll = min(100.0, payload->roll + 1); break;
      case UP: payload->pitch = min(100.0, payload->pitch + 1); break; // TODO: polarity?
      case DOWN: payload->pitch = max(-100.0, payload->pitch - 1); break;
      case FAST: payload->throttle = min(100.0, payload->throttle + 1); break;
      case SLOW: payload->throttle = max(-100.0, payload->throttle - 1); break;
      case GET: { double index = (fmod(fabs(NthMyChild(0)->evaluate ()), (double)(INPUT_SIZE)));
                  returnValue =  (isnan(index) || isinf(index)) ? 0 : payload->inputVar[(int) index];
                  break;
      }

    default: 
      GPExitSystem ("MyGene::evaluate", "Undefined node value");
    }
  
  return returnValue;
}



// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evaluate ()
{
  // Evaluate main tree
  std::uniform_real_distribution<double> dist(-1, 1);  //(min, max)
  //Mersenne Twister: Good quality random number generator
  std::mt19937 rng; 
  //Initialize with non-deterministic seeds
  rng.seed(std::random_device{}()); 

  stdFitness = 0;
  for (int i = 0; i < 1000; i++) {
    payload->inputVar[0] = dist(rng); // TODO Thread instance

    // what we should get
    double actual = sin(payload->inputVar[0] * M_PI);

    // eval
    payload->roll = 0;
    NthMyGene (0)->evaluate ();
    double found = payload->roll / 100.0;

    // how did we do?
    double delta = abs(actual - found);
    if (!isnan(delta) && !isinf(delta)) {
      stdFitness += (delta);
    } else {
      stdFitness += 1000000;
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
  ns.putNode (*new GPNode (SUB, "SUB", 2));
  ns.putNode (*new GPNode (MUL, "MUL", 2));
  ns.putNode (*new GPNode (DIV, "DIV", 2));
  ns.putNode (*new GPNode (REM, "REM", 2));
  ns.putNode (*new GPNode (ABS, "ABS", 1));
  ns.putNode (*new GPNode (NEG, "NEG", 1));
  ns.putNode (*new GPNode (MAX, "MAX", 2));
  ns.putNode (*new GPNode (MIN, "MIN", 2));
  ns.putNode (*new GPNode (IF, "IF", 3));
  ns.putNode (*new GPNode (GT, "GT", 2));
  ns.putNode (*new GPNode (LT, "LT", 2));
  ns.putNode (*new GPNode (EQ, "EQ", 2));
  ns.putNode (*new GPNode (ZERO, "ZERO"));
  ns.putNode (*new GPNode (ONE, "ONE"));
  ns.putNode (*new GPNode (LEFT, "LEFT", 0));
  ns.putNode (*new GPNode (RIGHT, "RIGHT", 0));
  ns.putNode (*new GPNode (UP, "UP", 0));
  ns.putNode (*new GPNode (DOWN, "DOWN", 0));
  ns.putNode (*new GPNode (FAST, "FAST", 0));
  ns.putNode (*new GPNode (SLOW, "SLOW", 0)); 
  ns.putNode (*new GPNode (GET, "GET", 1));
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
  ofstream fout (strOutFile.str());
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
  for (int gen=1; gen<=cfg.NumberOfGenerations; gen++)
    {
      // Create a new generation from the old one by applying the
      // genetic operators
      if (!cfg.SteadyState)
	      newPop=new MyPopulation (cfg, adfNs);
      pop->generate (*newPop);
      
      // Delete the old generation and make the new the old one
      if (!cfg.SteadyState)
	    {
	      delete pop;
	      pop=newPop;
	    }

      // Create a report of this generation and how well it is doing
      pop->createGenerationReport (0, gen, fout, bout);
    }

  return 0;
}

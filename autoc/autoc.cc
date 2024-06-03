
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
enum Operators {ADD=0, SUB, MUL, DIV, ABS, IF, GT, EQ, ROLL, PITCH, THROTTLE, GET, ADF0, ADF1, ADF2};
const int OPERATORS_NR_ITEM=15;
enum ADF0 {A_ADD=OPERATORS_NR_ITEM, A_MUL, A_NEG, A_INV, A_SIN, A_PI, A_UNITY, A_CROSS};
const int ADF0_NR_ITEM=8;


// Define class identifiers
const int MyGeneID=GPUserID;
const int MyGPID=GPUserID+1;
const int MyPopulationID=GPUserID+2;

Aircraft *aircraft = new Aircraft(new AircraftState());

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
                  returnValue = (divisor == 0) ? 0 : NthMyChild(0)->evaluate ()/divisor; 
                  break;
      }
      case ABS: returnValue=abs (NthMyChild(0)->evaluate ()); break;
      case IF: returnValue=NthMyChild(0)->evaluate () ? NthMyChild(1)->evaluate () : NthMyChild(2)->evaluate (); break;
      case GT: returnValue=NthMyChild(0)->evaluate () > NthMyChild(1)->evaluate (); break;
      case EQ: returnValue=NthMyChild(0)->evaluate () == NthMyChild(1)->evaluate (); break;
      case ROLL: returnValue = aircraft->setRollCommand(NthMyChild(0)->evaluate ()); break;
      case PITCH: returnValue = aircraft->setRollCommand(NthMyChild(0)->evaluate ()); break;
      case THROTTLE: returnValue = aircraft->setThrottleCommand(NthMyChild(0)->evaluate ()); break;
      case GET: returnValue = aircraft->getState()->dRelVel; break;
      case ADF0: returnValue = NthMyChild(0)->evaluate (); break;
      case ADF1: returnValue = NthMyChild(0)->evaluate (); break;
      case ADF2: returnValue = NthMyChild(0)->evaluate (); break;
      
      case A_ADD: returnValue=NthMyChild(0)->evaluate ()+NthMyChild(1)->evaluate (); break;
      case A_MUL: returnValue=NthMyChild(0)->evaluate ()*NthMyChild(1)->evaluate (); break;
      case A_INV: {
              double value = NthMyChild(0)->evaluate ();
              returnValue = value != 0 ? 1.0/value : 0;
              break;
      }
      case A_NEG: returnValue=-NthMyChild(0)->evaluate (); break;
      case A_SIN: returnValue=sin (NthMyChild(0)->evaluate ()); break;
      case A_PI: returnValue=M_PI; break;
      case A_UNITY: returnValue=1.0; break;
      case A_CROSS: returnValue=NthMyChild(0)->evaluate ()*NthMyChild(1)->evaluate (); break;

    default: 
      GPExitSystem ("MyGene::evaluate", "Undefined node value");
    }
  
  #define RANGELIMIT 1000000
  if (returnValue < -RANGELIMIT)
    return -RANGELIMIT;
  if (returnValue > RANGELIMIT)
    return RANGELIMIT;
  return returnValue;
}



// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evaluate ()
{
  // Evaluate main tree
  std::uniform_real_distribution<double> dist(-M_PI, M_PI);  //(min, max)
  //Mersenne Twister: Good quality random number generator
  std::mt19937 rng; 
  //Initialize with non-deterministic seeds
  rng.seed(std::random_device{}()); 

  stdFitness = 0;
  for (int i = 0; i < 1000; i++) {

    // AircraftState state = AircraftState(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
    // Aircraft aircraft = Aircraft(state);
    // aircraft.setPitchCommand(20);
    // aircraft.setRollCommand(30);
    // aircraft.setThrottleCommand(40);
    // char output[200];
    // aircraft.toString(output);
    // cout << output << endl;


    double input = dist(rng); // TODO Thread instance

    // what we should get
    double actual = sin(input);

    // eval
    aircraft->getState()->dRelVel = input;
    NthMyGene (0)->evaluate ();
    double found = aircraft->getRollCommand();

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
  adfNs.reserveSpace (4);
  
  // Now define the function and terminal set for each ADF and place
  // function/terminal sets into overall ADF container
  GPNodeSet& ns=*new GPNodeSet (OPERATORS_NR_ITEM);
  GPNodeSet& adf0=*new GPNodeSet (ADF0_NR_ITEM);
  GPNodeSet& adf1=*new GPNodeSet (ADF0_NR_ITEM);
  GPNodeSet& adf2=*new GPNodeSet (ADF0_NR_ITEM);

  adfNs.put (0, ns);
  adfNs.put (1, adf0);
  adfNs.put (2, adf1);
  adfNs.put (3, adf2);
  
  // Define functions/terminals and place them into the appropriate
  // sets.  Terminals take two arguments, functions three (the third
  // parameter is the number of arguments the function has)
  ns.putNode (*new GPNode (ADD, "ADD", 2));
  ns.putNode (*new GPNode (SUB, "SUB", 2));
  ns.putNode (*new GPNode (MUL, "MUL", 2));
  ns.putNode (*new GPNode (DIV, "DIV", 2));
  ns.putNode (*new GPNode (ABS, "ABS", 1));
  ns.putNode (*new GPNode (IF, "IF", 3));
  ns.putNode (*new GPNode (GT, "GT", 2));
  ns.putNode (*new GPNode (EQ, "EQ", 2));
  ns.putNode (*new GPNode (ROLL, "ROLL", 1));
  ns.putNode (*new GPNode (PITCH, "PITCH", 1));
  ns.putNode (*new GPNode (THROTTLE, "THROTTLE", 1));
  ns.putNode (*new GPNode (GET, "GET"));
  ns.putNode (*new GPNode (ADF0, "ADF0", 1));
  ns.putNode (*new GPNode (ADF1, "ADF1", 1));
  ns.putNode (*new GPNode (ADF2, "ADF2", 1));

  adf0.putNode (*new GPNode (A_ADD, "A_ADD", 2));
  adf0.putNode (*new GPNode (A_MUL, "A_MUL", 2));
  adf0.putNode (*new GPNode (A_INV, "A_INV", 1));
  adf0.putNode (*new GPNode (A_NEG, "A_NEG", 1));
  adf0.putNode (*new GPNode (A_SIN, "A_SIN", 1));
  adf0.putNode (*new GPNode (A_PI, "A_PI"));
  adf0.putNode (*new GPNode (A_UNITY, "A_UNITY"));
  adf0.putNode (*new GPNode (A_CROSS, "A_CROSS", 2));

  adf1.putNode (*new GPNode (A_ADD, "A_ADD", 2));
  adf1.putNode (*new GPNode (A_MUL, "A_MUL", 2));
  adf1.putNode (*new GPNode (A_INV, "A_INV", 1));
  adf1.putNode (*new GPNode (A_NEG, "A_NEG", 1));
  adf1.putNode (*new GPNode (A_SIN, "A_SIN", 1));
  adf1.putNode (*new GPNode (A_PI, "A_PI"));
  adf1.putNode (*new GPNode (A_UNITY, "A_UNITY"));
  adf1.putNode (*new GPNode (A_CROSS, "A_CROSS", 2));

  adf2.putNode (*new GPNode (A_ADD, "A_ADD", 2));
  adf2.putNode (*new GPNode (A_MUL, "A_MUL", 2));
  adf2.putNode (*new GPNode (A_INV, "A_INV", 1));
  adf2.putNode (*new GPNode (A_NEG, "A_NEG", 1));
  adf2.putNode (*new GPNode (A_SIN, "A_SIN", 1));
  adf2.putNode (*new GPNode (A_PI, "A_PI"));
  adf2.putNode (*new GPNode (A_UNITY, "A_UNITY"));
  adf2.putNode (*new GPNode (A_CROSS, "A_CROSS", 2));
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

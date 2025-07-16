#ifndef AUTOC_H
#define AUTOC_H

#define FITNESS_DISTANCE_WEIGHT 1.5
#define FITNESS_ALIGNMENT_WEIGHT 1.3
#define FITNESS_CONTROL_WEIGHT 1.0

class ExtraConfig {
public:
  int simNumPathsPerGen = 1;
  char* generatorMethod = "classic";
  int evalThreads = 1;
  char* minisimProgram = "../build/minisim";
  unsigned short minisimPortOverride = 0;
  char* s3Bucket = "autoc-storage";
  char* s3Profile = "default";
  int evaluateMode = 0;

  // // Custom implementation of the << operator for the extraCfg type
  // std::ostream& operator << (std::ostream& os) {
  //   os << "simNumPathsPerGen: " + simNumPathsPerGen;
  //   return os;
  // }
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

class MyGP;

extern GPAdfNodeSet adfNs;
extern void createNodeSet(GPAdfNodeSet& adfNs);
extern AircraftState aircraftState;
extern void initializeSimGP();

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
  virtual void evalTask(WorkerContext& context);
};



#endif
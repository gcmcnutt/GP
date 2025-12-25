#ifndef AUTOC_H
#define AUTOC_H

#include <vector>

#include "gp.h"

// Forward declarations
struct WorkerContext;

#define FITNESS_DISTANCE_WEIGHT 1.5f
#define FITNESS_ALIGNMENT_WEIGHT 1.3f
#define FITNESS_CONTROL_WEIGHT 1.0f

// Movement efficiency penalty constants - reduced for better exploration
#define CONTROL_SATURATION_THRESHOLD 0.9f
#define CONTROL_SATURATION_PENALTY 5.0f
#define CONTROL_RATE_PENALTY 2.0f
#define MOVEMENT_EFFICIENCY_WEIGHT 1.5f
#define MOVEMENT_EFFICIENCY_MAX_PENALTY 8.0f

// Cross-track and attitude penalties
#define CROSS_TRACK_WEIGHT 1.2f
#define CROSS_TRACK_SIGN_THRESHOLD 1.0f
#define CROSS_TRACK_OSC_WEIGHT 1.4f
#define ORIENTATION_TANGENT_WEIGHT 1.1f
#define ORIENTATION_UP_WEIGHT 1.3f
#define PATH_UP_INVERSION_PENALTY 75.0f

struct WindScenarioConfig {
  unsigned int windSeed = 0;
  int windVariantIndex = 0;
};

struct ScenarioDescriptor {
  std::vector<std::vector<Path>> pathList;
  std::vector<WindScenarioConfig> windScenarios;
  unsigned int windSeed = 0;
  int pathVariantIndex = 0;
  int windVariantIndex = 0;
};

class ExtraConfig {
public:
  int simNumPathsPerGen = 1;
  char* generatorMethod = "classic";
  int evalThreads = 1;
  char* minisimProgram = "../build/minisim";
  unsigned short minisimPortOverride = 0;
  char* s3Bucket = "autoc-storage";
  char* s3Profile = "default";
  int evaluateMode = 0;  // 0=normal GP evolution, 1=bytecode verification
  char* bytecodeFile = "gp_program.dat";  // Bytecode file for verification mode
  int windScenarioCount = 1;
  int windSeedBase = 1337;
  int windSeedStride = 1;
  unsigned int randomPathSeedB = 67890;  // Seed for SeededRandomB path generation

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
  GETALPHA, GETBETA, GETVELX, GETVELY, GETVELZ,
  GETROLL_RAD, GETPITCH_RAD,
  CLAMP, ATAN2, ABS, SQRT, MIN, MAX,
  OP_PI, ZERO, ONE, TWO, PROGN, _END
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
  gp_scalar evaluate(std::vector<Path>& path, MyGP& gp, gp_scalar arg);

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

  void setScenarioIndex(int idx) { scenarioIndex = idx; }
  int getScenarioIndex() const { return scenarioIndex; }
  void setBakeoffMode(bool mode) { bakeoffMode = mode; }
  bool isBakeoffMode() const { return bakeoffMode; }
  void setFitness(gp_scalar fitness) { stdFitness = fitness; fitnessValid = 1; }

  // async evaluator
  virtual void evalTask(WorkerContext& context);

private:
  int scenarioIndex = 0;
  bool bakeoffMode = false;
};



#endif

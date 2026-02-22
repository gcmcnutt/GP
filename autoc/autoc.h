#ifndef AUTOC_H
#define AUTOC_H

#include <vector>

#include "gp.h"
#include "gp_types.h"

// Forward declarations
struct WorkerContext;

// ========================================================================
// SIMPLIFIED FITNESS FUNCTION (v3 - two objectives, path-relative scaling)
// Goal: Stay close to rabbit, fly smooth. Strategy belongs in higher layers.
// See specs/FITNESS_SIMPLIFY_20260221.md for rationale.
// ========================================================================

// Distance penalty: nonlinear to keep small excursions expensive
#define DISTANCE_POWER 1.2

// Attitude scaling: computed per-path from path geometry (no manual tuning)
// attitude_scale = path_distance / path_turn_rad (meters per radian)

// Crash penalty: soft lexicographic multiplier
// Completion dominates (1e6 scale), quality provides gradient within similar completion levels
#define CRASH_COMPLETION_WEIGHT 1e6         // Multiplier for (1 - fraction_completed)

struct WindScenarioConfig {
  unsigned int windSeed = 0;
  int windVariantIndex = 0;
};

struct ScenarioDescriptor {
  std::vector<std::vector<Path>> pathList;
  std::vector<WindScenarioConfig> windScenarios;
  unsigned int windSeed = 0;
  int pathVariantIndex = -1;   // -1 = unset/aggregated
  int windVariantIndex = -1;   // -1 = unset/aggregated
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
  // windSeedBase and windSeedStride removed - wind seeds now derived from GPrand() (single PRNG architecture)
  int randomPathSeedB = 67890;  // Seed for path generation (-1 = derive from GPrand())
  int gpSeed = -1;  // Seed for GP initialization (-1 = use time-based seed)
  char* trainingNodes = "";  // Comma-separated list of node names for training (empty = all nodes)

  // VARIATIONS1: Entry and wind direction variations (see specs/VARIATIONS1.md)
  int enableEntryVariations = 0;  // 0=disabled, 1=enabled (requires crrcsim support)
  int enableWindVariations = 0;   // 0=disabled, 1=enabled (requires crrcsim support)
  // Sigma values in degrees (for documentation/logging, actual values in variation_generator.h)
  double entryHeadingSigma = 45.0;
  double entryRollSigma = 22.5;
  double entryPitchSigma = 7.5;
  double entrySpeedSigma = 0.1;
  double windDirectionSigma = 45.0;

  // Variable rabbit speed (see specs/VARIABLE_RABBIT.md)
  // Set rabbitSpeedSigma=0 for constant speed at rabbitSpeedNominal
  double rabbitSpeedNominal = 16.0;   // m/s - center of distribution (default matches SIM_RABBIT_VELOCITY)
  double rabbitSpeedSigma = 0.0;      // m/s - 1σ deviation (0 = constant speed)
  double rabbitSpeedMin = 8.0;        // m/s - hard floor
  double rabbitSpeedMax = 25.0;       // m/s - hard ceiling
  double rabbitSpeedCycleMin = 0.5;   // seconds - min variation cycle duration
  double rabbitSpeedCycleMax = 5.0;   // seconds - max variation cycle duration

  // // Custom implementation of the << operator for the extraCfg type
  // std::ostream& operator << (std::ostream& os) {
  //   os << "simNumPathsPerGen: " + simNumPathsPerGen;
  //   return os;
  // }
};

// Define function and terminal identifiers
// NOTE: New operators must be added at the END (before _END) to preserve
// backward compatibility with serialized bytecode programs.
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
  OP_PI, ZERO, ONE, TWO, PROGN,
  // Temporal state terminals (added 2026-02) - see specs/TEMPORAL_STATE.md
  GETDPHI_PREV, GETDTHETA_PREV,   // History lookback (1 arg: tick index)
  GETDPHI_RATE, GETDTHETA_RATE,   // Error derivatives (nullary)
  _END
};
const int OPERATORS_NR_ITEM = _END;

// Define class identifiers
const int MyGeneID = GPUserID;
const int MyGPID = GPUserID + 1;
const int MyPopulationID = GPUserID + 2;

class MyGP;

extern GPAdfNodeSet adfNs;
extern void createNodeSet(GPAdfNodeSet& adfNs);
extern void setTrainingNodesMask(const char* mask);  // Call before createNodeSet() to filter nodes
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
  MyGP(MyGP& gpo) : GP(gpo),
                     scenarioIndex(gpo.scenarioIndex),
                     bakeoffMode(gpo.bakeoffMode) { }
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
  void setFitness(gp_fitness fitness) { stdFitness = fitness; fitnessValid = 1; }

  // async evaluator
  virtual void evalTask(WorkerContext& context);

private:
  int scenarioIndex = 0;
  bool bakeoffMode = false;
};



#endif

#ifndef AUTOC_H
#define AUTOC_H

#include <vector>

#include "gp.h"
#include "gp_types.h"

// Forward declarations
struct WorkerContext;

// ========================================================================
// SIMPLIFIED FITNESS FUNCTION (v2 - raw units, no percentage scaling)
// Goal: Smooth control to intercept and track path with minimum energy
// All metrics in raw units (meters, J/kg) with power-law weighting
// Step fitness normalized by time (per-second basis at 10Hz sampling)
// ========================================================================

// Time normalization: scale step fitness to per-second basis
// At 10Hz (100ms steps), each step contributes 0.1 seconds worth of error
#define STEP_TIME_WEIGHT (SIM_TIME_STEP_MSEC / 1000.0)  // 0.1 at 10Hz

// Primary objective: Reach waypoints on time (meters)
// Asymmetric based on position relative to rabbit's velocity vector
// Ahead (in front of rabbit's path tangent): harsher penalty (overshooting)
// Behind (trailing rabbit): lighter penalty (following correctly)
#define WAYPOINT_AHEAD_POWER 1.8            // Power when ahead of rabbit (worse)
#define WAYPOINT_BEHIND_POWER 1.2           // Power when behind rabbit (acceptable)

// Secondary objective: Move efficiently along path
#define MOVEMENT_DIRECTION_WEIGHT 1.3       // Power for direction alignment
#define DIRECTION_SCALE 20.0                // Scale 0-2 direction error to ~meter equivalent

// Tertiary objective: Match rabbit's energy state (kinetic + potential)
// Energy deviation in J/kg (= m^2/s^2), asymmetric by altitude
// NED coords: altitude = -z, so craft_alt < rabbit_alt means craft is BELOW
#define ALTITUDE_LOW_POWER 1.5              // Power when below target (harder to recover)
#define ALTITUDE_HIGH_POWER 1.0             // Power when above target (easy to correct)

// Path normal alignment: aircraft orientation error relative to coordinated flight
// For coordinated turns, aircraft "up" should point toward inside of turn circle
// For straight flight, aircraft should be wings-level or transitioning smoothly
// Error measured as angle (radians) between desired and actual aircraft normal
#define NORMAL_ERROR_POWER 1.5              // Power for normal alignment error
#define NORMAL_ERROR_SCALE 10.0             // Scale factor: normalize ~0.5 rad error to ~meter-equivalent

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
  int windSeedBase = 1337;
  int windSeedStride = 1;
  int randomPathSeedB = 67890;  // Seed for SeededRandomB path generation (-1 = use time-based seed)
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

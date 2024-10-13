#ifndef AUTOC_H
#define AUTOC_H

#define FITNESS_DISTANCE_WEIGHT 1.5
#define FITNESS_ALIGNMENT_WEIGHT 1.3
#define FITNESS_CONTROL_WEIGHT 1.0

class ExtraConfig {
public:
  int simNumPathsPerGen = 1;
  int evalThreads = 1;
  char* minisimProgram = "../build/minisim";
  unsigned short minisimPortOverride = 0;
  char* s3Bucket = "autoc-storage";

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

extern std::atomic_bool printEval;
extern std::vector<MyGP*> tasks;
extern std::vector<std::vector<Path>> generationPaths;
extern std::ofstream fout;
extern std::atomic_ulong nanDetector;
extern void createNodeSet(GPAdfNodeSet& adfNs);
extern EvalResults evalResults;

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

  AircraftState aircraftState{ 0, 0, Eigen::Quaterniond::Identity(), Eigen::Vector3d(0, 0, 0), 0.0, 0.0, 0.0, 0, false };
  long pathIndex = 0; // current entry on path

  // TODO this is unrelated to MyGP
  void printAircraftState(std::ostream& os, EvalResults &evalResults) {
    for (int i = 0; i < evalResults.aircraftStateList.size(); i++) {
      os << "  Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP controlP\n";

      for (int j = 0; j < evalResults.aircraftStateList.at(i).size(); j++) {
        AircraftState stepAircraftState = evalResults.aircraftStateList.at(i).at(j);
        Path path = evalResults.pathList.at(i).at(stepAircraftState.getThisPathIndex());
        Eigen::Vector3d euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);

        char outbuf[1000]; // XXX use c++20
        sprintf(outbuf, "%06ld %3ld % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n",
          stepAircraftState.getSimTime(), pathIndex,
          path.distanceFromStart,
          path.start[0],
          path.start[1],
          path.start[2],
          stepAircraftState.getPosition()[0],
          stepAircraftState.getPosition()[1],
          stepAircraftState.getPosition()[2],
          euler[2],
          euler[1],
          euler[0],
          stepAircraftState.getRelVel(),
          stepAircraftState.getRollCommand(),
          stepAircraftState.getPitchCommand(),
          stepAircraftState.getThrottleCommand()
        );
        os << outbuf;
      }
    }
  }
};



#endif
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
  char* sqsUrl = "https://sqs.us-west-2.amazonaws.com/499918285206/autoc-tasks";
  char* s3Bucket = "autoc-storage";

  // // Custom implementation of the << operator for the extraCfg type
  // std::ostream& operator << (std::ostream& os) {
  //   os << "simNumPathsPerGen: " + simNumPathsPerGen;
  //   return os;
  // }
};

enum class CrashReason {
  None,
  Sim,
  Eval,
};

#endif
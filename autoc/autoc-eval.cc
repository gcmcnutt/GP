
// autoc-eval.cc

#include <vector>
#include <math.h>
#include <set>
#include <string>
#include <sstream>

#include "gp.h"
#include "minisim.h"
#include "autoc.h"
#include "gp_evaluator_portable.h"

GPAdfNodeSet adfNs;

// Track which nodes are enabled for training (for logging purposes)
std::vector<std::string> trainingNodeNames;

// Training node mask - set by autoc.cc before calling createNodeSet()
// Empty string means use all nodes (default for minisim, crrcsim, xiao-gp)
static const char* gTrainingNodesMask = "";

void setTrainingNodesMask(const char* mask) {
  gTrainingNodesMask = mask ? mask : "";
}

// TODO really get rid of this global!
AircraftState aircraftState;

std::string crashReasonToString(CrashReason type) {
  switch (type) {
  case CrashReason::None: return "None";
  case CrashReason::Boot: return "Boot";
  case CrashReason::Sim: return "Sim";
  case CrashReason::Eval: return "Eval";
  case CrashReason::Time: return "Time";
  case CrashReason::Distance: return "Distance";
  default: return "*?*";
  }
};

// initialize the GP engine
void initializeSimGP() {

  // Create the adf function/terminal set and print it out.
  createNodeSet(adfNs);

  // partial registration ( TODO may be easier just to load all in clients )
  GPRegisterClass(new GPContainer());
  GPRegisterClass(new GPNode());
  GPRegisterClass(new GPNodeSet());
  GPRegisterClass(new GPAdfNodeSet());
  GPRegisterClass(new GPGene());
  GPRegisterClass(new GP());

  // manually add our classes for load operation
  GPRegisterClass(new MyGene());
  GPRegisterClass(new MyGP());
}

// Parse comma-separated node names into a set
static std::set<std::string> parseTrainingMask(const char* nodeList) {
  std::set<std::string> mask;
  if (!nodeList || nodeList[0] == '\0') {
    return mask;  // Empty = no filtering (all nodes enabled)
  }
  std::istringstream ss(nodeList);
  std::string node;
  while (std::getline(ss, node, ',')) {
    // Trim whitespace
    size_t start = node.find_first_not_of(" \t");
    size_t end = node.find_last_not_of(" \t");
    if (start != std::string::npos && end != std::string::npos) {
      mask.insert(node.substr(start, end - start + 1));
    }
  }
  return mask;
}

// Node definition table for consistent iteration
struct NodeDef {
  int opcode;
  const char* name;
  int args;  // 0 = terminal
};

static const NodeDef allNodes[] = {
  {ADD, "ADD", 2},
  {SUB, "SUB", 2},
  {MUL, "MUL", 2},
  {DIV, "DIV", 2},
  {IF, "IF", 3},
  {EQ, "EQ", 2},
  {GT, "GT", 2},
  {SETPITCH, "SETPITCH", 1},
  {SETROLL, "SETROLL", 1},
  {SETTHROTTLE, "SETTHROTTLE", 1},
  {GETPITCH, "GETPITCH", 0},
  {GETROLL, "GETROLL", 0},
  {GETTHROTTLE, "GETTHROTTLE", 0},
  {SIN, "SIN", 1},
  {COS, "COS", 1},
  {OP_PI, "OP_PI", 0},
  {ZERO, "0", 0},
  {ONE, "1", 0},
  {TWO, "2", 0},
  {PROGN, "PROGN", 2},
  {GETDPHI, "GETDPHI", 1},
  {GETDTHETA, "GETDTHETA", 1},
  {GETDTARGET, "GETDTARGET", 1},
  {GETVEL, "GETVEL", 0},
  {GETDHOME, "GETDHOME", 0},
  {GETALPHA, "GETALPHA", 0},
  {GETBETA, "GETBETA", 0},
  {GETVELX, "GETVELX", 0},
  {GETVELY, "GETVELY", 0},
  {GETVELZ, "GETVELZ", 0},
  {GETROLL_RAD, "GETROLL_RAD", 0},
  {GETPITCH_RAD, "GETPITCH_RAD", 0},
  {CLAMP, "CLAMP", 3},
  {ATAN2, "ATAN2", 2},
  {ABS, "ABS", 1},
  {SQRT, "SQRT", 1},
  {MIN, "MIN", 2},
  {MAX, "MAX", 2},
};
static const int allNodesCount = sizeof(allNodes) / sizeof(allNodes[0]);

// Create function and terminal set
void createNodeSet(GPAdfNodeSet& adfNs)
{
  // Get training mask (set by setTrainingNodesMask() before this call)
  std::set<std::string> mask = parseTrainingMask(gTrainingNodesMask);
  bool filterEnabled = !mask.empty();

  // Clear the training node names list
  trainingNodeNames.clear();

  // First pass: count how many nodes will be added
  int nodeCount = 0;
  for (int i = 0; i < allNodesCount; i++) {
    if (!filterEnabled || mask.count(allNodes[i].name) > 0) {
      nodeCount++;
    }
  }

  // Reserve space for the node sets
  adfNs.reserveSpace(1);

  // Create GPNodeSet with exact size needed
  GPNodeSet& ns = *new GPNodeSet(nodeCount);

  adfNs.put(0, ns);

  // Second pass: add the nodes
  for (int i = 0; i < allNodesCount; i++) {
    const NodeDef& def = allNodes[i];
    if (!filterEnabled || mask.count(def.name) > 0) {
      if (def.args > 0) {
        ns.putNode(*new GPNode(def.opcode, const_cast<char*>(def.name), def.args));
      } else {
        ns.putNode(*new GPNode(def.opcode, const_cast<char*>(def.name)));
      }
      trainingNodeNames.push_back(def.name);
    }
  }
}



// Legacy getIndex function - now delegates to portable implementation
int getIndex(std::vector<Path>& path, MyGP& gp, gp_scalar arg) {
  VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
  return getPathIndex(pathProvider, aircraftState, arg);
}

// This function evaluates the fitness of a genetic tree using the portable evaluator.
// This maintains compatibility while delegating to the portable implementation.
gp_scalar MyGene::evaluate(std::vector<Path>& path, MyGP& run, gp_scalar arg)
{
  // Create path provider for this evaluation
  VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
  
  // For leaf nodes (terminals), evaluate directly
  if (containerSize() == 0) {
    return evaluateGPOperator(node->value(), pathProvider, aircraftState, nullptr, 0, static_cast<gp_scalar>(arg));
  }
  
  // For function nodes, evaluate children first
  gp_scalar childArgs[3] = {0.0f, 0.0f, 0.0f}; // Max 3 args (IF, CLAMP)
  int numArgs = containerSize();
  
  for (int i = 0; i < numArgs && i < 3; i++) {
    childArgs[i] = static_cast<gp_scalar>(NthMyChild(i)->evaluate(path, run, arg));
  }
  
  return evaluateGPOperator(node->value(), pathProvider, aircraftState, childArgs, numArgs, static_cast<gp_scalar>(arg));
}

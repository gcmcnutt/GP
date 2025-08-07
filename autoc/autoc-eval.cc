
// autoc-eval.cc

#include <vector>
#include <math.h>

#include "gp.h"
#include "minisim.h"
#include "autoc.h"
#include "gp_evaluator_portable.h"

GPAdfNodeSet adfNs;

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

// Create function and terminal set
void createNodeSet(GPAdfNodeSet& adfNs)
{
  // Reserve space for the node sets
  adfNs.reserveSpace(1);

  // Now define the function and terminal set for each ADF and place
  // function/terminal sets into overall ADF container
  GPNodeSet& ns = *new GPNodeSet(OPERATORS_NR_ITEM);

  adfNs.put(0, ns);

  // Define functions/terminals and place them into the appropriate
  // sets.  Terminals take two arguments, functions three (the third
  // parameter is the number of arguments the function has)
  ns.putNode(*new GPNode(ADD, "ADD", 2));
  ns.putNode(*new GPNode(SUB, "SUB", 2));
  ns.putNode(*new GPNode(MUL, "MUL", 2));
  ns.putNode(*new GPNode(DIV, "DIV", 2));
  ns.putNode(*new GPNode(IF, "IF", 3));
  ns.putNode(*new GPNode(EQ, "EQ", 2));
  ns.putNode(*new GPNode(GT, "GT", 2));
  ns.putNode(*new GPNode(SETPITCH, "SETPITCH", 1));
  ns.putNode(*new GPNode(SETROLL, "SETROLL", 1));
  ns.putNode(*new GPNode(SETTHROTTLE, "SETTHROTTLE", 1));
  ns.putNode(*new GPNode(GETPITCH, "GETPITCH"));
  ns.putNode(*new GPNode(GETROLL, "GETROLL"));
  ns.putNode(*new GPNode(GETTHROTTLE, "GETTHROTTLE"));
  ns.putNode(*new GPNode(SIN, "SIN", 1));
  ns.putNode(*new GPNode(COS, "COS", 1));
  ns.putNode(*new GPNode(OP_PI, "OP_PI"));
  ns.putNode(*new GPNode(ZERO, "0"));
  ns.putNode(*new GPNode(ONE, "1"));
  ns.putNode(*new GPNode(TWO, "2"));
  ns.putNode(*new GPNode(PROGN, "PROGN", 2));
  ns.putNode(*new GPNode(GETDPHI, "GETDPHI", 1));
  ns.putNode(*new GPNode(GETDTHETA, "GETDTHETA", 1));
  ns.putNode(*new GPNode(GETDTARGET, "GETDTARGET", 1));
  ns.putNode(*new GPNode(GETVEL, "GETVEL"));
  ns.putNode(*new GPNode(GETDHOME, "GETDHOME"));
  ns.putNode(*new GPNode(GETALPHA, "GETALPHA"));
  ns.putNode(*new GPNode(GETBETA, "GETBETA"));
  ns.putNode(*new GPNode(GETVELX, "GETVELX"));
  ns.putNode(*new GPNode(GETVELY, "GETVELY"));
  ns.putNode(*new GPNode(GETVELZ, "GETVELZ"));
  ns.putNode(*new GPNode(GETROLL_RAD, "GETROLL_RAD"));
  ns.putNode(*new GPNode(GETPITCH_RAD, "GETPITCH_RAD"));
  ns.putNode(*new GPNode(CLAMP, "CLAMP", 3));
  ns.putNode(*new GPNode(ATAN2, "ATAN2", 2));
  ns.putNode(*new GPNode(ABS, "ABS", 1));
  ns.putNode(*new GPNode(SQRT, "SQRT", 1));
  ns.putNode(*new GPNode(MIN, "MIN", 2));
  ns.putNode(*new GPNode(MAX, "MAX", 2));
}



// Legacy getIndex function - now delegates to portable implementation
int getIndex(std::vector<Path>& path, MyGP& gp, double arg) {
  VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
  return getPathIndex(pathProvider, aircraftState, arg);
}

// This function evaluates the fitness of a genetic tree using the portable evaluator.
// This maintains compatibility while delegating to the portable implementation.
double MyGene::evaluate(std::vector<Path>& path, MyGP& run, double arg)
{
  // Create path provider for this evaluation
  VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
  
  // For leaf nodes (terminals), evaluate directly
  if (containerSize() == 0) {
    return evaluateGPOperator(node->value(), pathProvider, aircraftState, nullptr, 0, arg);
  }
  
  // For function nodes, evaluate children first
  double childArgs[3] = {0.0, 0.0, 0.0}; // Max 3 args (IF, CLAMP)
  int numArgs = containerSize();
  
  for (int i = 0; i < numArgs && i < 3; i++) {
    childArgs[i] = NthMyChild(i)->evaluate(path, run, arg);
  }
  
  return evaluateGPOperator(node->value(), pathProvider, aircraftState, childArgs, numArgs, arg);
}

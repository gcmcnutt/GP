
// autoc-eval.cc

#include <vector>
#include <math.h>

#include "gp.h"
#include "minisim.h"
#include "autoc.h"

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
  ns.putNode(*new GPNode(PI, "PI"));
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



int getIndex(std::vector<Path>& path, MyGP& gp, double arg) {
  if (isnan(arg)) {
    return aircraftState.getThisPathIndex();
  }

  // how many rabbit steps are we estimating?
  int steps = std::clamp((int)arg, -5, 5);

  // determine the current distance and then the distance to the selected point as a function of rabbit speed
  double distanceSoFar = (path.at(aircraftState.getThisPathIndex()).distanceFromStart);
  double distanceGoal = distanceSoFar + steps * SIM_RABBIT_VELOCITY * (SIM_TIME_STEP_MSEC / 1000.0);

  // ok, we have a goal and our current position, walk up or down the paths to find the closest point
  int currentStep = aircraftState.getThisPathIndex();
  if (steps > 0) {
    while (path.at(currentStep).distanceFromStart < distanceGoal && currentStep < path.size() - 1) {
      currentStep++;
    }
  }
  else {
    while (path.at(currentStep).distanceFromStart > distanceGoal && currentStep > 0) {
      currentStep--;
    }
  }
  return currentStep;
}

// This function evaluates the fitness of a genetic tree.  We have the
// freedom to define this function in any way we like.
double MyGene::evaluate(std::vector<Path>& path, MyGP& run, double arg)
{
  double returnValue = 0.0;

  switch (node->value())
  {
  case ADD: returnValue = NthMyChild(0)->evaluate(path, run, arg) + NthMyChild(1)->evaluate(path, run, arg); break;
  case SUB: returnValue = -NthMyChild(0)->evaluate(path, run, arg) - NthMyChild(1)->evaluate(path, run, arg); break;
  case MUL: returnValue = NthMyChild(0)->evaluate(path, run, arg) * NthMyChild(1)->evaluate(path, run, arg); break;
  case DIV: {
    double dividend = NthMyChild(0)->evaluate(path, run, arg);
    double divisor = NthMyChild(1)->evaluate(path, run, arg);
    returnValue = (divisor == 0) ? 0 : dividend / divisor;
    break;
  }
  case IF: returnValue = NthMyChild(0)->evaluate(path, run, arg) ? NthMyChild(1)->evaluate(path, run, arg) : NthMyChild(2)->evaluate(path, run, arg); break;
  case EQ: returnValue = NthMyChild(0)->evaluate(path, run, arg) == NthMyChild(1)->evaluate(path, run, arg); break;
  case GT: returnValue = NthMyChild(0)->evaluate(path, run, arg) > NthMyChild(1)->evaluate(path, run, arg); break;
  case SETPITCH: returnValue = aircraftState.setPitchCommand(NthMyChild(0)->evaluate(path, run, arg)); break;
  case SETROLL: returnValue = aircraftState.setRollCommand(NthMyChild(0)->evaluate(path, run, arg)); break;
  case SETTHROTTLE: returnValue = aircraftState.setThrottleCommand(NthMyChild(0)->evaluate(path, run, arg)); break;
  case GETPITCH: returnValue = aircraftState.getPitchCommand(); break;
  case GETROLL: returnValue = aircraftState.getRollCommand(); break;
  case GETTHROTTLE: returnValue = aircraftState.getThrottleCommand(); break;
  case SIN: returnValue = sin(NthMyChild(0)->evaluate(path, run, arg)); break;
  case COS: returnValue = cos(NthMyChild(0)->evaluate(path, run, arg)); break;
  case PI: returnValue = M_PI; break;
  case ZERO: returnValue = 0; break;
  case ONE: returnValue = 1; break;
  case TWO: returnValue = 2; break;
  case PROGN: {
    NthMyChild(0)->evaluate(path, run, arg);
    returnValue = NthMyChild(1)->evaluate(path, run, arg);
    break;
  }
  case GETVEL: returnValue = aircraftState.getRelVel(); break;

  case GETDPHI: // compute roll goal from current to target
  {
    // Calculate the vector from craft to target in world frame
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    Eigen::Vector3d craftToTarget = path.at(idx).start - aircraftState.getPosition();

    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;

    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());

    // Calculate the angle between the projected vector and the body Z-axis
    returnValue = std::atan2(projectedVector.y(), -projectedVector.z());
    break;
  }

  case GETDTHETA: // compute pitch goal from current to target assuming we did a good roll
  {
    // Calculate the vector from craft to target in world frame
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    Eigen::Vector3d craftToTarget = path.at(idx).start - aircraftState.getPosition();

    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;

    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());

    // Calculate the angle between the projected vector and the body Z-axis
    double rollEstimate = std::atan2(projectedVector.y(), -projectedVector.z());

    // *** PITCH: Calculate the vector from craft to target in world frame if it did rotate
    Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
    Eigen::Quaterniond virtualOrientation = aircraftState.getOrientation() * rollRotation;

    // Transform target vector to new virtual orientation
    Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;

    // Calculate pitch angle
    returnValue = std::atan2(-newLocalTargetVector.z(), newLocalTargetVector.x());
    break;
  }

  case GETDTARGET: // get distance to the next point
  {
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    double distance = (path.at(idx).start - aircraftState.getPosition()).norm();
    returnValue = std::clamp((distance - 10) / aircraftState.getRelVel(), -1.0, 1.0);
    break;
  }

  case GETDHOME: // get distance to the home point
  {
    returnValue = (Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE) - aircraftState.getPosition()).norm();
    break;
  }

  case GETALPHA: // get angle of attack
  {
    // Transform actual velocity vector to body frame
    Eigen::Vector3d velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
    
    // Angle of attack is angle between velocity and body X axis (forward)
    // atan2(vertical_component, forward_component)
    returnValue = std::atan2(-velocity_body.z(), velocity_body.x());
    break;
  }

  case GETBETA: // get sideslip angle
  {
    // Transform actual velocity vector to body frame
    Eigen::Vector3d velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
    
    // Sideslip is angle between velocity and body XZ plane
    // atan2(side_component, forward_component)
    returnValue = std::atan2(velocity_body.y(), velocity_body.x());
    break;
  }

  case GETVELX: // get velocity X component (North)
  {
    returnValue = aircraftState.getVelocity().x();
    break;
  }

  case GETVELY: // get velocity Y component (East)
  {
    returnValue = aircraftState.getVelocity().y();
    break;
  }

  case GETVELZ: // get velocity Z component (Down/climb rate)
  {
    returnValue = aircraftState.getVelocity().z();
    break;
  }

  case GETROLL_RAD: // get roll angle in radians
  {
    // Convert quaternion to Euler angles using standard aerospace convention
    Eigen::Vector3d euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
    returnValue = euler[2]; // Roll angle (rotation around X-axis)
    break;
  }

  case GETPITCH_RAD: // get pitch angle in radians  
  {
    // Convert quaternion to Euler angles using standard aerospace convention
    Eigen::Vector3d euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
    returnValue = euler[1]; // Pitch angle (rotation around Y-axis)
    break;
  }

  case CLAMP: // clamp value between min and max
  {
    double value = NthMyChild(0)->evaluate(path, run, arg);
    double minVal = NthMyChild(1)->evaluate(path, run, arg);
    double maxVal = NthMyChild(2)->evaluate(path, run, arg);
    returnValue = CLAMP_DEF(value, minVal, maxVal);
    break;
  }

  case ATAN2: // atan2(y, x)
  {
    double y = NthMyChild(0)->evaluate(path, run, arg);
    double x = NthMyChild(1)->evaluate(path, run, arg);
    returnValue = ATAN2_DEF(y, x);
    break;
  }

  case ABS: // absolute value
  {
    double value = NthMyChild(0)->evaluate(path, run, arg);
    returnValue = ABS_DEF(value);
    break;
  }

  case SQRT: // square root
  {
    double value = NthMyChild(0)->evaluate(path, run, arg);
    returnValue = (value >= 0) ? SQRT_DEF(value) : 0.0;
    break;
  }

  case MIN: // minimum of two values
  {
    double a = NthMyChild(0)->evaluate(path, run, arg);
    double b = NthMyChild(1)->evaluate(path, run, arg);
    returnValue = MIN_DEF(a, b);
    break;
  }

  case MAX: // maximum of two values
  {
    double a = NthMyChild(0)->evaluate(path, run, arg);
    double b = NthMyChild(1)->evaluate(path, run, arg);
    returnValue = MAX_DEF(a, b);
    break;
  }

  default:
    GPExitSystem("MyGene::evaluate", "Undefined node value");
  }

#define RANGELIMIT 1000000
  if (returnValue < -RANGELIMIT)
    return -RANGELIMIT;
  if (returnValue > RANGELIMIT)
    return RANGELIMIT;
  if (abs(returnValue) < 0.000001)
    return 0;
  return returnValue;
}

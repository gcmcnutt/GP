
// autoc-eval.cc

#include <boost/asio.hpp>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <new>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>

#include "gp.h"
#include "minisim.h"
#include "autoc.h"

// TODO really get rid of these globals!
std::atomic_bool printEval = false; // verbose (used for rendering best of population)
std::vector<MyGP*> tasks = std::vector<MyGP*>();
std::vector<std::vector<Path>> generationPaths;
std::atomic_ulong nanDetector = 0;
EvalResults evalResults;

std::string crashReasonToString(CrashReason type) {
  switch (type) {
  case CrashReason::None: return "None";
  case CrashReason::Boot: return "Boot";
  case CrashReason::Sim: return "Sim";
  case CrashReason::Eval: return "Eval";
  default: return "*?*";
  }
};

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
}



int getIndex(std::vector<Path>& path, MyGP& gp, double arg) {
  if (isnan(arg)) {
    return gp.pathIndex;
  }

  // a range of steps to check, can't go lower than the beginning index
  // TODO this checks path next, not the actual simulation steps...
  // XXX for now, this allows forecasting the future path
  int idx = std::clamp((int)arg, -5, 5) + gp.pathIndex;
  idx = std::clamp(idx, 0, (int)path.size() - 1);
  return idx;
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
  case SETPITCH: returnValue = /*run.aircraftState.setPitchCommand(*/NthMyChild(0)->evaluate(path, run, arg); break;
  case SETROLL: returnValue = /*run.aircraftState.setRollCommand(*/NthMyChild(0)->evaluate(path, run, arg); break;
  case SETTHROTTLE: returnValue = /*run.aircraftState.setThrottleCommand(*/NthMyChild(0)->evaluate(path, run, arg); break;
  case GETPITCH: returnValue = run.aircraftState.getPitchCommand(); break;
  case GETROLL: returnValue = run.aircraftState.getRollCommand(); break;
  case GETTHROTTLE: returnValue = run.aircraftState.getThrottleCommand(); break;
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
  case GETVEL: returnValue = run.aircraftState.getRelVel(); break;

  case GETDPHI: // compute roll goal from current to target
  {
    // Calculate the vector from craft to target in world frame
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    Eigen::Vector3d craftToTarget = path.at(idx).start - run.aircraftState.getPosition();

    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = run.aircraftState.getOrientation().inverse() * craftToTarget;

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
    Eigen::Vector3d craftToTarget = path.at(idx).start - run.aircraftState.getPosition();

    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = run.aircraftState.getOrientation().inverse() * craftToTarget;

    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());

    // Calculate the angle between the projected vector and the body Z-axis
    double rollEstimate = std::atan2(projectedVector.y(), -projectedVector.z());

    // *** PITCH: Calculate the vector from craft to target in world frame if it did rotate
    Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
    Eigen::Quaterniond virtualOrientation = run.aircraftState.getOrientation() * rollRotation;

    // Transform target vector to new virtual orientation
    Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;

    // Calculate pitch angle
    returnValue = std::atan2(-newLocalTargetVector.z(), newLocalTargetVector.x());
    break;
  }

  case GETDTARGET: // get distance to the next point
  {
    int idx = getIndex(path, run, NthMyChild(0)->evaluate(path, run, arg));
    double distance = (path.at(idx).start - run.aircraftState.getPosition()).norm();
    returnValue = std::clamp((distance - 10) / run.aircraftState.getRelVel(), -1.0, 1.0);
    break;
  }

  case GETDHOME: // get distance to the home point
  {
    returnValue = (Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE) - run.aircraftState.getPosition()).norm();
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

void MyGP::evaluate()
{
  tasks.push_back(this);
}

// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evalTask(WorkerContext& context)
{
  stdFitness = 0;

  // TODO this save is probably cheaper when GP knows about boost archives...
  std::vector<char> buffer;
  boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> outStream(buffer);
  save(outStream);
  outStream.flush();
  EvalData evalData = { buffer, generationPaths };
  sendRPC(*context.socket, evalData);

  // How did it go?
  evalResults = receiveRPC<EvalResults>(*context.socket);

  // Compute the fitness results for each path
  for (int i = 0; i < evalResults.pathList.size(); i++) {
    bool printHeader = true;

    // get path, actual and aircraft state
    auto& path = evalResults.pathList.at(i);
    auto& aircraftState = evalResults.aircraftStateList.at(i);
    auto& crashReason = evalResults.crashReasonList.at(i);

    // compute this path fitness
    double localFitness = 0;
    int stepIndex = 0; // where are we on the flight path?

    // error accumulators
    double distance_error_sum = 0;
    double angle_error_sum = 0;
    double control_smoothness_sum = 0;
    int simulation_steps = 0;

    // initial states
    auto& stepAircraftState = aircraftState.at(stepIndex++);
    double roll_prev = stepAircraftState.getRollCommand();
    double pitch_prev = stepAircraftState.getPitchCommand();
    double throttle_prev = stepAircraftState.getThrottleCommand();

    std::vector<Eigen::Vector3d> planPath;
    std::vector<Eigen::Vector3d> actualPath;
    std::vector<Eigen::Vector3d> actualOrientation;

    // now walk next steps of actual path
    while (stepIndex++ < aircraftState.size()) {
      auto& stepAircraftState = aircraftState.at(stepIndex++);
      int pathIndex = stepAircraftState.getThisPathIndex();

      // Compute the distance between the aircraft and the goal
      double distanceFromGoal = (path.at(pathIndex).start - stepAircraftState.getPosition()).norm();
      // normalize [100:0]
      distanceFromGoal = distanceFromGoal * 100.0 / SIM_PATH_RADIUS_LIMIT;

      // Compute vector from me to target
      Eigen::Vector3d target_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start);
      Eigen::Vector3d aircraft_to_target = (path.at(pathIndex).start - stepAircraftState.getPosition());
      double dot_product = target_direction.dot(aircraft_to_target);
      double angle_rad = std::acos(std::clamp(dot_product / (target_direction.norm() * aircraft_to_target.norm()), -1.0, 1.0));
      // normalize [100:0]
      angle_rad = angle_rad * 100.0 / M_PI;

      // control smoothness -- internal vector distance
      double smoothness = pow(roll_prev - stepAircraftState.getRollCommand(), 2.0);
      smoothness += pow(pitch_prev - stepAircraftState.getPitchCommand(), 2.0);
      smoothness += pow(throttle_prev - stepAircraftState.getThrottleCommand(), 2.0);
      smoothness = sqrt(smoothness);
      // normalize [100:0]
      smoothness = smoothness * 100.0 / 3.46;

      // ready for next cycle
      roll_prev = stepAircraftState.getRollCommand();
      pitch_prev = stepAircraftState.getPitchCommand();
      throttle_prev = stepAircraftState.getThrottleCommand();

      // accumulate the error
      distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
      angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
      control_smoothness_sum += pow(smoothness, FITNESS_CONTROL_WEIGHT);
      simulation_steps++;
    }

    // tally up the normlized fitness based on steps and progress
    double normalized_distance_error = (distance_error_sum / simulation_steps);
    double normalized_velocity_align = (angle_error_sum / simulation_steps);
    double normalized_control_smoothness = (control_smoothness_sum / simulation_steps);
    localFitness = normalized_distance_error + normalized_velocity_align + normalized_control_smoothness;

    if (isnan(localFitness)) {
      nanDetector++;
    }

    if (crashReason != CrashReason::None) {
      double fractional_distance_remaining = 1.0 - path.at(pathIndex).distanceFromStart / path.back().distanceFromStart;
      localFitness += SIM_CRASH_PENALTY * fractional_distance_remaining;
    }

    if (printEval) {
      // now update actual list of lists
      std::vector<Path> actualElementList;
      for (int i = 0; i < actualPath.size(); i++) {
        Eigen::Vector3d path = actualPath.at(i);
        Eigen::Vector3d orientation = actualOrientation.at(i);
        actualElementList.push_back({ path, orientation, 0, 0 });
      }

      // now update plan list of lists
      std::vector<Path> planElementList;
      for (auto& planElement : planPath) {
        planElementList.push_back({ planElement, Eigen::Vector3d::UnitX(), 0, 0 });
      }

      planPath.clear();
      actualPath.clear();
      actualOrientation.clear();
    }

    // accumulate the local fitness
    stdFitness += localFitness;
  }

  // normalize
  stdFitness /= evalResults.pathList.size();
}

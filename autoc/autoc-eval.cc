
// autoc-eval.cc

/* -------------------------------------------------------------------
From skeleton/skeleton.cc
------------------------------------------------------------------- */

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
#include "gpconfig.h"
#include "minisim.h"
#include "autoc.h"
#include "pathgen.h"
#include "threadpool.h"
#include "logger.h"

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/PutObjectRequest.h>

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/support/date_time.hpp>

std::mutex evalMutex;

std::atomic_ulong nanDetector(0);

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
  case SETPITCH: returnValue = run.aircraftState.setPitchCommand(NthMyChild(0)->evaluate(path, run, arg)); break;
  case SETROLL: returnValue = run.aircraftState.setRollCommand(NthMyChild(0)->evaluate(path, run, arg)); break;
  case SETTHROTTLE: returnValue = run.aircraftState.setThrottleCommand(NthMyChild(0)->evaluate(path, run, arg)); break;
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


// TODO move this outside of MyGP
void MyGP::evaluate()
{
  tasks.push_back(this);
}



// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evalTask(WorkerContext& context)
{
  stdFitness = 0;

  for (int i = 0; i < generationPaths.size(); i++) {
    auto& path = generationPaths.at(i);
    double localFitness = 0;
    std::vector<Eigen::Vector3d> planPath = std::vector<Eigen::Vector3d>();
    std::vector<Eigen::Vector3d> actualPath = std::vector<Eigen::Vector3d>();
    std::vector<Eigen::Vector3d> actualOrientation = std::vector<Eigen::Vector3d>();

    // deal with pre.path on the initial eval..
    if (path.size() == 0) {
      stdFitness = 1000001;
      continue;
    }

    // random initial orientation
    Eigen::Quaterniond aircraft_orientation;
    Eigen::Vector3d initialPosition;
    {
      // Guard random functions
      std::lock_guard<std::mutex> guard(evalMutex);
      aircraft_orientation = Eigen::Quaterniond::UnitRandom();

      initialPosition = Eigen::Vector3d((((double)GPrand() / RAND_MAX) - 0.5) * SIM_INITIAL_LOCATION_DITHER,
        (((double)GPrand() / RAND_MAX) - 0.5) * SIM_INITIAL_LOCATION_DITHER,
        SIM_INITIAL_ALTITUDE - ((double)GPrand() / RAND_MAX) * SIM_INITIAL_LOCATION_DITHER);

      // override orientation of aircraft to be upright and level
      aircraft_orientation = Eigen::Quaterniond(Eigen::AngleAxisd(M_PI, Eigen::Vector3d::UnitZ())) *
        Eigen::Quaterniond(Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY())) *
        Eigen::Quaterniond(Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX()));

      // override position of aircraft to be at the origin
      initialPosition = Eigen::Vector3d(-2.19, 5.49, -36.93);
    }

    // wait for generic request (we ignore)
    receiveRPC<AircraftState>(*context.socket);

    // send initial state reset
    aircraftState = AircraftState{ SIM_INITIAL_VELOCITY, aircraft_orientation, initialPosition, 0.0, 0.0, SIM_INITIAL_THROTTLE, 0, false };
    sendRPC(*context.socket, MainToSim{ ControlType::AIRCRAFT_STATE, aircraftState });

    // iterate the simulator
    unsigned long int duration_msec = 0; // how long have we been running
    pathIndex = 0; // where are we on the path?
    bool printHeader = true;

    // as long as we are within the time limit and have not reached the end of the path
    CrashReason crashReason = CrashReason::None;

    // error accumulators
    double distance_error_sum = 0;
    double angle_error_sum = 0;
    double control_smoothness_sum = 0;
    int simulation_steps = 0;

    // initial states
    double roll_prev = aircraftState.getRollCommand();
    double pitch_prev = aircraftState.getPitchCommand();
    double throttle_prev = aircraftState.getThrottleCommand();

    while (duration_msec < SIM_TOTAL_TIME_MSEC && pathIndex < path.size() - 2 && crashReason == CrashReason::None) {
      // wait for next state update from sim
      aircraftState = receiveRPC<AircraftState>(*context.socket);

      // did sim crash?
      if (aircraftState.isSimCrashed()) {
        crashReason = CrashReason::Sim;
      }

      // ---------------------- first how did we do?
      // Compute the distance between the aircraft and the goal
      double distanceFromGoal = (path.at(pathIndex).start - aircraftState.getPosition()).norm();
      // normalize [100:0]
      distanceFromGoal = distanceFromGoal * 100.0 / SIM_PATH_RADIUS_LIMIT;

      // Compute vector from me to target
      Eigen::Vector3d target_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start);
      Eigen::Vector3d aircraft_to_target = (path.at(pathIndex).start - aircraftState.getPosition());
      double dot_product = target_direction.dot(aircraft_to_target);
      double angle_rad = std::acos(std::clamp(dot_product / (target_direction.norm() * aircraft_to_target.norm()), -1.0, 1.0));
      // normalize [100:0]
      angle_rad = angle_rad * 100.0 / M_PI;

      // control smoothness -- internal vector distance
      double smoothness = pow(roll_prev - aircraftState.getRollCommand(), 2.0);
      smoothness += pow(pitch_prev - aircraftState.getPitchCommand(), 2.0);
      smoothness += pow(throttle_prev - aircraftState.getThrottleCommand(), 2.0);
      smoothness = sqrt(smoothness);
      // normalize [100:0]
      smoothness = smoothness * 100.0 / 3.46;

      // ready for next cycle
      roll_prev = aircraftState.getRollCommand();
      pitch_prev = aircraftState.getPitchCommand();
      throttle_prev = aircraftState.getThrottleCommand();

      // but eval detected crashed outside the cylinder? (low elevation is detected from sim)
      double distanceFromOrigin = std::sqrt(aircraftState.getPosition()[0] * aircraftState.getPosition()[0] +
        aircraftState.getPosition()[1] * aircraftState.getPosition()[1]);
      if (aircraftState.getPosition()[2] < (SIM_MIN_ELEVATION - SIM_PATH_RADIUS_LIMIT) || // too high
        aircraftState.getPosition()[2] > (SIM_MIN_ELEVATION) || // too low
        distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) { // too far
        crashReason = CrashReason::Eval;
      }

      if (printEval) {
        if (printHeader) {
          fout << "  Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP controlP\n";
          printHeader = false;
        }

        // convert aircraft_orientaton to euler
        Eigen::Vector3d euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);

        char outbuf[1000]; // XXX use c++20
        sprintf(outbuf, "%06ld %3ld % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n",
          aircraftState.getSimTime(), pathIndex,
          path.at(pathIndex).distanceFromStart,
          path.at(pathIndex).start[0],
          path.at(pathIndex).start[1],
          path.at(pathIndex).start[2],
          aircraftState.getPosition()[0],
          aircraftState.getPosition()[1],
          aircraftState.getPosition()[2],
          euler[2],
          euler[1],
          euler[0],
          aircraftState.getRelVel(),
          aircraftState.getRollCommand(),
          aircraftState.getPitchCommand(),
          aircraftState.getThrottleCommand(),
          distanceFromGoal,
          angle_rad,
          smoothness
        );
        fout << outbuf;

        // now update points for Renderer
        planPath.push_back(path.at(pathIndex).start); // XXX push the full path
        actualPath.push_back(aircraftState.getPosition());

        // negative here as world up is negative Z
        actualOrientation.push_back(aircraftState.getOrientation() * -Eigen::Vector3d::UnitZ());
      }

      distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
      angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
      control_smoothness_sum += pow(smoothness, FITNESS_CONTROL_WEIGHT);
      simulation_steps++;

      // TODO ensure time is forward
      unsigned long int duration_msec = aircraftState.getSimTime();

      // search for location of next timestamp
      double timeDistance = SIM_RABBIT_VELOCITY * duration_msec / 1000.0;
      while (pathIndex < path.size() - 2 && (path.at(pathIndex).distanceFromStart < timeDistance)) {
        pathIndex++;
      }

#if 0
      // ESTIMATE: pre-compute the estimated roll, pitch, throttle in prep for the evaluator

      // *** ROLL: Calculate the vector from craft to target in world frame
      Eigen::Vector3d craftToTarget = path.at(pathIndex).start - aircraftState.getPosition();

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
      double pitchEstimate = std::atan2(-newLocalTargetVector.z(), newLocalTargetVector.x());

      // // now try to determine if pitch up or pitch down makes more sense
      // if (std::abs(pitchEstimate) > M_PI / 2) {
      //   pitchEstimate = (pitchEstimate > 0) ? pitchEstimate - M_PI : pitchEstimate + M_PI;
      //   rollEstimate = -rollEstimate;
      // }

      // range is -1:1
      aircraftState.setRollCommand(rollEstimate / M_PI);
      aircraftState.setPitchCommand(pitchEstimate / M_PI);

      // Throttle estimate range is -1:1
      {
        double distance = (path.at(pathIndex).start - aircraftState.getPosition()).norm();
        double throttleEstimate = std::clamp((distance - 10) / aircraftState.getRelVel(), -1.0, 1.0);
        aircraftState.setThrottleCommand(throttleEstimate);
      }

      {
        char outbuf[1000];
        Eigen::Matrix3d r = aircraftState.getOrientation().toRotationMatrix();
        sprintf(outbuf, "Estimate airXaxis:%f %f %f  airZaxis:%f %f %f toTarget:%f %f %f rolledTarget:%f %f %f  r:%f p:%f t:%f\n",
          r.col(0).x(), r.col(0).y(), r.col(0).z(),
          r.col(2).x(), r.col(2).y(), r.col(2).z(),
          projectedVector.x(), projectedVector.y(), projectedVector.z(),
          newLocalTargetVector.x(), newLocalTargetVector.y(), newLocalTargetVector.z(),
          aircraftState.getPitchCommand(), aircraftState.getRollCommand(), aircraftState.getThrottleCommand());

        *logger.debug() << outbuf;
      }
#endif

      // GP determine control input adjustments
      NthMyGene(0)->evaluate(path, *this, 0);

      // reply with the commands
      MainToSim mainToSim{ ControlType::CONTROL_SIGNAL, ControlSignal{aircraftState.getPitchCommand(), aircraftState.getRollCommand(), aircraftState.getThrottleCommand()} };
      sendRPC(*context.socket, mainToSim);
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
      evalResults.actualList.push_back(actualElementList);
      evalResults.pathList.push_back(planElementList);

      planPath.clear();
      actualPath.clear();
      actualOrientation.clear();
    }

    stdFitness += localFitness;

    *logger.debug() << "MyGP: " << this << " path[" << i << "] complete." << endl;
  }

  // normalize
  stdFitness /= generationPaths.size();
}


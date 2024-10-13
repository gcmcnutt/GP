/* test sim for aircraft */
#include <boost/asio.hpp>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <new>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <iostream>
#include <cmath>
#include <array>

#include "gp.h"
#include "minisim.h"
#include "autoc.h"

using namespace std;
using boost::asio::ip::tcp;

boost::iostreams::stream<boost::iostreams::array_source> charArrayToIstream(const std::vector<char>& charArray) {
  return boost::iostreams::stream<boost::iostreams::array_source>(
    boost::iostreams::array_source(charArray.data(), charArray.size())
  );
}

GPAdfNodeSet adfNs;

class Aircraft {
public:
  // TODO should be private
  double dRelVel; // reltive forward velocity on +x airplane axis m/s

  // world frame for now
  Eigen::Quaterniond aircraft_orientation;

  // NED convention for location x+ north, y+ east, z+ down
  Eigen::Vector3d position;

  // not used yet
  double R_X;     // rotationX
  double R_Y;     // rotationY
  double R_Z;     // rotationZ

  Aircraft(double dRelVel, Eigen::Quaterniond aircraft_orientation, Eigen::Vector3d position, double R_X, double R_Y, double R_Z);

  double setPitchCommand(double pitchCommand);
  double getPitchCommand();
  double setRollCommand(double rollCommand);
  double getRollCommand();
  double setThrottleCommand(double throttleCommand);
  double getThrottleCommand();
  void advanceState(double dt);
  void toString(char* output);

private:

  // aircraft command values
  double pitchCommand;  // -1:1
  double rollCommand;   // -1:1
  double throttleCommand; // -1:1
};

Aircraft::Aircraft(double dRelVel, Eigen::Quaterniond aircraft_orientation, Eigen::Vector3d position, double R_X, double R_Y, double R_Z) {
  this->dRelVel = dRelVel;
  this->aircraft_orientation = aircraft_orientation;
  this->position = position;
  this->R_X = R_X;
  this->R_Y = R_Y;
  this->R_Z = R_Z;

  this->pitchCommand = 0;
  this->rollCommand = 0;
  this->throttleCommand = 0;
}

double Aircraft::setPitchCommand(double pitchCommand) {
  this->pitchCommand = pitchCommand;
  return pitchCommand;
}

double Aircraft::getPitchCommand() {
  return std::clamp(pitchCommand, -1.0, 1.0);
}

double Aircraft::setRollCommand(double rollCommand) {
  this->rollCommand = rollCommand;
  return rollCommand;
}

double Aircraft::getRollCommand() {
  return std::clamp(rollCommand, -1.0, 1.0);
}

double Aircraft::setThrottleCommand(double throttleCommand) {
  this->throttleCommand = throttleCommand;
  return throttleCommand;
}

double Aircraft::getThrottleCommand() {
  return std::clamp(throttleCommand, -1.0, +1.0);
}

void Aircraft::advanceState(double dt) {
  // get current roll state, compute left/right force (positive roll is right)
  double delta_roll = remainder(getRollCommand() * dt * SIM_MAX_ROLL_RATE_RADSEC, M_PI);

  // get current pitch state, compute up/down force (positive pitch is up)
  double delta_pitch = remainder(getPitchCommand() * dt * SIM_MAX_PITCH_RATE_RADSEC, M_PI);

  // adjust velocity as a function of throttle (-1:1)
  dRelVel = SIM_INITIAL_VELOCITY + (getThrottleCommand() * SIM_THROTTLE_SCALE);

  // Convert pitch and roll updates to quaternions (in the body frame)
  Eigen::Quaterniond delta_roll_quat(Eigen::AngleAxisd(delta_roll, Eigen::Vector3d::UnitX()));
  Eigen::Quaterniond delta_pitch_quat(Eigen::AngleAxisd(delta_pitch, Eigen::Vector3d::UnitY()));

  // Apply the roll and pitch adjustments to the aircraft's orientation
  aircraft_orientation = aircraft_orientation * delta_roll_quat;
  aircraft_orientation = aircraft_orientation * delta_pitch_quat;

  // Normalize the resulting quaternion
  aircraft_orientation.normalize();

  // Define the initial velocity vector in the body frame
  Eigen::Vector3d velocity_body(dRelVel * dt, 0, 0);

  // Rotate the velocity vector using the updated quaternion
  Eigen::Vector3d velocity_world = aircraft_orientation * velocity_body;

  // adjust position
  position += velocity_world;
}

void Aircraft::toString(char* output) {
  sprintf(output, "AircraftState: %f %f %f %f %f %f %f %f %f %f %f  Command: %f %f %f\n", dRelVel,
    aircraft_orientation.w(), aircraft_orientation.x(), aircraft_orientation.y(), aircraft_orientation.z(),
    position[0], position[1], position[2], R_X, R_Y, R_Z,
    pitchCommand, rollCommand, throttleCommand);
}


class SimProcess {
public:
  Aircraft aircraft = Aircraft(0.0, Eigen::Quaterniond(1, 0, 0, 0), Eigen::Vector3d(0, 0, 0), 0.0, 0.0, 0.0);
  unsigned long int simTime = 0;
  bool simCrashed = false;
  SimProcess(boost::asio::io_context& io_context, unsigned short port) : socket_(io_context) {
    tcp::resolver resolver(io_context);
    auto endpoints = resolver.resolve("localhost", std::to_string(port));
    boost::asio::connect(socket_, endpoints);

    // send checkin call
    sendRPC(socket_, EvalResults{});
  }

  void run() {
    while (true) { // TODO have reply have a controlled loop exit
      EvalResults evalResults;

      // ok what does main say to do
      EvalData evalData = receiveRPC<EvalData>(socket_);

      // for each path, evaluate
      for (int i = 0; i < evalData.pathList.size(); i++) {
        std::vector<Path> path = evalData.pathList.at(i);
        MyGP gp;
        boost::iostreams::stream<boost::iostreams::array_source> is = charArrayToIstream(evalData.gp);
        gp.load(is);
        gp.resolveNodeValues(adfNs);

        // accumulators
        std::vector<Path> actualPath;
        std::vector<AircraftState> aircraftState;

        // random initial orientation
        Eigen::Quaterniond aircraft_orientation;
        Eigen::Vector3d initialPosition;
        {
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

        // reset sim state
        gp.aircraftState = AircraftState{ 0, SIM_INITIAL_VELOCITY, aircraft_orientation, initialPosition, 0.0, 0.0, SIM_INITIAL_THROTTLE, 0, false };
        gp.pathIndex = 0;

        // iterate the simulator
        unsigned long int duration_msec = 0; // how long have we been running
        CrashReason crashReason = CrashReason::None;

        // as long as we are within the time limit and have not reached the end of the path
        while (duration_msec < SIM_TOTAL_TIME_MSEC && gp.pathIndex < path.size() - 2 && crashReason == CrashReason::None) {
          gp.aircraftState.setThisPathIndex(gp.pathIndex);

          // approximate pitch/roll/throttle to achieve goal

          // *** ROLL: Calculate the vector from craft to target in world frame
          Eigen::Vector3d craftToTarget = path.at(gp.pathIndex).start - gp.aircraftState.getPosition();

          // Transform the craft-to-target vector to body frame
          Eigen::Vector3d target_local = gp.aircraftState.getOrientation().inverse() * craftToTarget;

          // Project the craft-to-target vector onto the body YZ plane
          Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());

          // Calculate the angle between the projected vector and the body Z-axis 
          double rollEstimate = std::atan2(projectedVector.y(), -projectedVector.z());

          // *** PITCH: Calculate the vector from craft to target in world frame if it did rotate
          Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
          Eigen::Quaterniond virtualOrientation = gp.aircraftState.getOrientation() * rollRotation;

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
          gp.aircraftState.setRollCommand(rollEstimate / M_PI);
          gp.aircraftState.setPitchCommand(pitchEstimate / M_PI);

          // Throttle estimate range is -1:1
          {
            double distance = (path.at(gp.pathIndex).start - gp.aircraftState.getPosition()).norm();
            double throttleEstimate = std::clamp((distance - 10) / gp.aircraftState.getRelVel(), -1.0, 1.0);
            gp.aircraftState.setThrottleCommand(throttleEstimate);
          }

          {
            char outbuf[1000];
            Eigen::Matrix3d r = gp.aircraftState.getOrientation().toRotationMatrix();
            sprintf(outbuf, "Estimate airXaxis:%f %f %f  airZaxis:%f %f %f toTarget:%f %f %f rolledTarget:%f %f %f  r:%f p:%f t:%f\n",
              r.col(0).x(), r.col(0).y(), r.col(0).z(),
              r.col(2).x(), r.col(2).y(), r.col(2).z(),
              projectedVector.x(), projectedVector.y(), projectedVector.z(),
              newLocalTargetVector.x(), newLocalTargetVector.y(), newLocalTargetVector.z(),
              gp.aircraftState.getPitchCommand(), gp.aircraftState.getRollCommand(), gp.aircraftState.getThrottleCommand());

            std::cout << outbuf;
          }

          // run the GP controller
          gp.NthMyGene(0)->evaluate(path, gp, 0);

          // advance the aircraft
          aircraft.setPitchCommand(gp.aircraftState.getPitchCommand());
          aircraft.setRollCommand(gp.aircraftState.getRollCommand());
          aircraft.setThrottleCommand(gp.aircraftState.getThrottleCommand());
          aircraft.advanceState(SIM_TIME_STEP_MSEC / 1000.0);

          // send our updated state to evaluator
          gp.aircraftState.setRelVel(aircraft.dRelVel);
          gp.aircraftState.setOrientation(aircraft.aircraft_orientation);
          gp.aircraftState.setPosition(aircraft.position);
          gp.aircraftState.setSimTime(gp.aircraftState.getSimTime() + SIM_TIME_STEP_MSEC);

          // did we crash?
          double distanceFromOrigin = std::sqrt(gp.aircraftState.getPosition()[0] * gp.aircraftState.getPosition()[0] +
            gp.aircraftState.getPosition()[1] * gp.aircraftState.getPosition()[1]);
          if (gp.aircraftState.getPosition()[2] < (SIM_MIN_ELEVATION - SIM_PATH_RADIUS_LIMIT) || // too high
            gp.aircraftState.getPosition()[2] > (SIM_MIN_ELEVATION) || // too low
            distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) { // too far
            crashReason = CrashReason::Eval;
          }

          // TODO ensure time is forward
          unsigned long int duration_msec = gp.aircraftState.getSimTime();

          // search for location of next timestamp
          double timeDistance = SIM_RABBIT_VELOCITY * duration_msec / 1000.0;
          while (gp.pathIndex < path.size() - 2 && (path.at(gp.pathIndex).distanceFromStart < timeDistance)) {
            gp.pathIndex++;
          }

          // record progress
          aircraftState.push_back(gp.aircraftState);
        }

        // save results of this run
        evalResults.pathList.push_back(path);
        evalResults.aircraftStateList.push_back(aircraftState);
        evalResults.crashReasonList.push_back(crashReason);
      }

      // always send our state -- covers initial state
      sendRPC(socket_, evalResults);

      // prepare results for next cycle
      evalResults.pathList.clear();
      evalResults.aircraftStateList.clear();
      evalResults.crashReasonList.clear();
    }
  }

private:
  tcp::socket socket_;
};

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: minisim <port>" << std::endl;
    return 1;
  }

  // Create the adf function/terminal set and print it out.
  createNodeSet(adfNs);

  // partial registration ( TODO may be easier just to load all in clients )
  GPRegisterClass(new GPContainer());
  GPRegisterClass(new GPNode());
  GPRegisterClass(new GPNodeSet());
  GPRegisterClass(new GPAdfNodeSet());
  // GPRegisterClass (new GPVariables());
  GPRegisterClass(new GPGene());
  GPRegisterClass(new GP());
  // GPRegisterClass (new GPPopulation());

  // manually add our classes for load operation
  GPRegisterClass(new MyGene());
  GPRegisterClass(new MyGP());

  unsigned short port = std::atoi(argv[1]);
  boost::asio::io_context io_context;
  SimProcess sim_process(io_context, port);
  sim_process.run();

  return 0;
}

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

class SimProcess {
public:
  SimProcess(boost::asio::io_context& io_context, unsigned short port) : socket_(io_context) {
    tcp::resolver resolver(io_context);
    auto endpoints = resolver.resolve("localhost", std::to_string(port));
    boost::asio::connect(socket_, endpoints);
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

        // accumulate steps
        std::vector<AircraftState> aircraftStateSteps;

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
        aircraftState = AircraftState{ 0, SIM_INITIAL_VELOCITY, aircraft_orientation, initialPosition, 0.0, 0.0, SIM_INITIAL_THROTTLE, 0 };

        // iterate the simulator
        unsigned long int duration_msec = 0; // how long have we been running
        CrashReason crashReason = CrashReason::None;

        while (crashReason == CrashReason::None) {

          // approximate pitch/roll/throttle to achieve goal

          // *** ROLL: Calculate the vector from craft to target in world frame
          Eigen::Vector3d craftToTarget = path.at(aircraftState.getThisPathIndex()).start - aircraftState.getPosition();

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
            double distance = (path.at(aircraftState.getThisPathIndex()).start - aircraftState.getPosition()).norm();
            double throttleEstimate = std::clamp((distance - 10) / aircraftState.getRelVel(), -1.0, 1.0);
            aircraftState.setThrottleCommand(throttleEstimate);
          }

#if 0
          {
            char outbuf[1000];
            Eigen::Matrix3d r = aircraftState.getOrientation().toRotationMatrix();
            sprintf(outbuf, "Estimate airXaxis:%f %f %f  airZaxis:%f %f %f toTarget:%f %f %f rolledTarget:%f %f %f  r:%f p:%f t:%f\n",
              r.col(0).x(), r.col(0).y(), r.col(0).z(),
              r.col(2).x(), r.col(2).y(), r.col(2).z(),
              projectedVector.x(), projectedVector.y(), projectedVector.z(),
              newLocalTargetVector.x(), newLocalTargetVector.y(), newLocalTargetVector.z(),
              aircraftState.getPitchCommand(), aircraftState.getRollCommand(), aircraftState.getThrottleCommand());

            std::cout << outbuf;
          }
#endif

          // run the GP controller
          gp.NthMyGene(0)->evaluate(path, gp, 0);

          // advance the aircraft
          aircraftState.minisimAdvanceState(SIM_TIME_STEP_MSEC);
          duration_msec += SIM_TIME_STEP_MSEC;
          aircraftState.setSimTimeMsec(duration_msec);

          // did we crash?
          double distanceFromOrigin = std::sqrt(aircraftState.getPosition()[0] * aircraftState.getPosition()[0] +
            aircraftState.getPosition()[1] * aircraftState.getPosition()[1]);
          if (aircraftState.getPosition()[2] < (SIM_MIN_ELEVATION - SIM_PATH_RADIUS_LIMIT) || // too high
            aircraftState.getPosition()[2] > (SIM_MIN_ELEVATION) || // too low
            distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) { // too far
            crashReason = CrashReason::Eval;
          }

          // search for location of next timestamp
          double timeDistance = SIM_RABBIT_VELOCITY * duration_msec / 1000.0;
          while (aircraftState.getThisPathIndex() < path.size() - 2 && (path.at(aircraftState.getThisPathIndex()).distanceFromStart < timeDistance)) {
            aircraftState.setThisPathIndex(aircraftState.getThisPathIndex() + 1);
          }

          // record progress
          aircraftStateSteps.push_back(aircraftState);

          // as long as we are within the time limit and have not reached the end of the path
          if (duration_msec >= SIM_TOTAL_TIME_MSEC) {
            crashReason = CrashReason::Time;
          }

          if (aircraftState.getThisPathIndex() >= path.size() - 2) {
            crashReason = CrashReason::Distance;
          }
        }

        // save results of this run
        evalResults.pathList.push_back(path);
        evalResults.aircraftStateList.push_back(aircraftStateSteps);
        evalResults.crashReasonList.push_back(crashReason);
      }

      // always send our state -- covers initial state
      // evalResults.dump(std::cout);
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
  if (argc != 4) {
    std::cerr << "Usage: minisim <base> <id> <port>" << std::endl;
    return 1;
  }
std:string id = argv[2];

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

  unsigned short port = std::atoi(argv[3]);
  boost::asio::io_context io_context;
  SimProcess sim_process(io_context, port);
  sim_process.run();

  return 0;
}

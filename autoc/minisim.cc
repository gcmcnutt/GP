/* test sim for aircraft */
#include <boost/asio.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
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
#include "gp_bytecode.h"
#include "fastmath/metrics.h"

using namespace std;
using boost::asio::ip::tcp;

namespace {

void ensureScenarioMetadata(EvalData& evalData) {
  if (evalData.scenarioList.size() == evalData.pathList.size()) {
    return;
  }
  evalData.scenarioList.assign(evalData.pathList.size(), evalData.scenario);
  for (size_t idx = 0; idx < evalData.scenarioList.size(); ++idx) {
    if (evalData.scenarioList[idx].pathVariantIndex < 0) {
      evalData.scenarioList[idx].pathVariantIndex = static_cast<int>(idx);
    }
  }
}

ScenarioMetadata scenarioForPathIndex(const EvalData& evalData, size_t idx) {
  if (idx < evalData.scenarioList.size()) {
    return evalData.scenarioList.at(idx);
  }
  ScenarioMetadata meta = evalData.scenario;
  if (meta.pathVariantIndex < 0) {
    meta.pathVariantIndex = static_cast<int>(idx);
  }
  return meta;
}

} // namespace

boost::iostreams::stream<boost::iostreams::array_source> charArrayToIstream(const std::vector<char>& charArray) {
  return boost::iostreams::stream<boost::iostreams::array_source>(
    boost::iostreams::array_source(charArray.data(), charArray.size())
  );
}

void MyGP::evaluate() {}
void MyGP::evalTask(WorkerContext& context) {}

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

      ensureScenarioMetadata(evalData);

      // flip this back for return trip
      evalResults.gp = evalData.gp;
      if (!evalData.scenarioList.empty()) {
        evalResults.scenario = evalData.scenarioList.front();
      } else {
        evalResults.scenario = evalData.scenario;
      }
      evalResults.scenarioList.clear();
      evalResults.scenarioList.reserve(evalData.scenarioList.size());

      // Detect if we have GP tree data or bytecode data
      bool isGPTreeData = false;
      bool isBytecodeData = false;
      MyGP gp;
      GPBytecodeInterpreter interpreter;
      
      // Check if the data starts with binary archive header (binary archives have specific magic bytes)
      // Binary archives start with specific byte patterns, while GP tree data is typically text-based
      bool looksLikeBinaryArchive = false;
      if (evalData.gp.size() >= 4) {
        // Check for binary archive magic bytes (boost binary archive typically starts with specific patterns)
        unsigned char firstBytes[4] = {
          (unsigned char)evalData.gp[0], 
          (unsigned char)evalData.gp[1], 
          (unsigned char)evalData.gp[2], 
          (unsigned char)evalData.gp[3]
        };
        // Binary archives often start with version info in binary format
        looksLikeBinaryArchive = (firstBytes[0] >= 0x16 && firstBytes[0] <= 0x20) && 
                                (firstBytes[1] == 0x00 || firstBytes[2] == 0x00);
      }
      
      if (looksLikeBinaryArchive) {
        // This looks like Boost binary serialized bytecode data
        try {
          boost::iostreams::stream<boost::iostreams::array_source> bytecodeStream = charArrayToIstream(evalData.gp);
          boost::archive::binary_iarchive archive(bytecodeStream);
          archive >> interpreter;
          isBytecodeData = true;
        } catch (const std::exception& e) {
          std::cerr << "Error loading bytecode data: " << e.what() << std::endl;
          continue;
        }
      } else {
        // This looks like GP tree data
        try {
          boost::iostreams::stream<boost::iostreams::array_source> gpStream = charArrayToIstream(evalData.gp);
          gp.load(gpStream);
          gp.resolveNodeValues(adfNs);
          isGPTreeData = true;
        } catch (const std::exception& e) {
          std::cerr << "Error loading GP tree data: " << e.what() << std::endl;
          continue;
        }
      }
      
      // for each path, evaluate
      for (int i = 0; i < evalData.pathList.size(); i++) {
        std::vector<Path> path = evalData.pathList.at(i);

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

        // compute initial velocity vector based on aircraft orientation (keeping minisim simple)
        Eigen::Vector3d initial_velocity = aircraft_orientation * Eigen::Vector3d(SIM_INITIAL_VELOCITY, 0, 0);
        
        // reset sim state
        aircraftState = AircraftState{ 0, SIM_INITIAL_VELOCITY, initial_velocity, aircraft_orientation, initialPosition, 0.0, 0.0, SIM_INITIAL_THROTTLE, 0 };

        // Record initial aircraft state at time 0 to match path start
        aircraftStateSteps.push_back(aircraftState);

        // iterate the simulator
        unsigned long int duration_msec = 0; // how long have we been running
        CrashReason crashReason = CrashReason::None;

        while (crashReason == CrashReason::None) {

          // BASELINE CONTROLLER DISABLED FOR GP-ONLY LEARNING (consistent with crrcsim)
          // Control commands start at 0 for each path (consistent with crrcsim reset behavior)
          // Within each path evaluation, commands persist between timesteps for incremental GP adjustments

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

          // Store control values before evaluation for comparison
          double pre_roll = aircraftState.getRollCommand();
          double pre_pitch = aircraftState.getPitchCommand();
          double pre_throttle = aircraftState.getThrottleCommand();
          
          // run the controller (GP tree or bytecode interpreter) - BOTH NOW HAVE SAME BASELINE
          double evaluation_result = 0.0;
          if (isGPTreeData) {
            evaluation_result = gp.NthMyGene(0)->evaluate(path, gp, 0);
          } else if (isBytecodeData) {
            // Use bytecode interpreter to evaluate and set control commands
            // NOTE: Bytecode interpreter now starts with proper baseline estimates set above
            evaluation_result = interpreter.evaluate(aircraftState, path, 0.0);
            // Note: The bytecode interpreter can modify control commands via SETPITCH, SETROLL, SETTHROTTLE opcodes
          }
          

          // advance the aircraft
          aircraftState.minisimAdvanceState(SIM_TIME_STEP_MSEC);
          duration_msec += SIM_TIME_STEP_MSEC;
          aircraftState.setSimTimeMsec(duration_msec);

          // did we crash?
          double distanceFromOrigin = std::sqrt(aircraftState.getPosition()[0] * aircraftState.getPosition()[0] +
            aircraftState.getPosition()[1] * aircraftState.getPosition()[1]);
          if (aircraftState.getPosition()[2] < (SIM_MAX_ELEVATION) || // too high
            aircraftState.getPosition()[2] > (SIM_MIN_ELEVATION) || // too low
            distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) { // too far
            crashReason = CrashReason::Eval;
          }

          // search for location of next timestamp using time-based targeting
          while (aircraftState.getThisPathIndex() < path.size() - 2 && (path.at(aircraftState.getThisPathIndex()).simTimeMsec < duration_msec)) {
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
        {
          ScenarioMetadata meta = scenarioForPathIndex(evalData, static_cast<size_t>(i));
          evalResults.scenarioList.push_back(meta);
        }
        evalResults.pathList.push_back(path);
        evalResults.aircraftStateList.push_back(aircraftStateSteps);
        evalResults.crashReasonList.push_back(crashReason);
      }

      // always send our state -- covers initial state
      // evalResults.dump(std::cout);
      sendRPC(socket_, evalResults);
#ifdef GP_FASTMATH_TRACE
      double evalCount = static_cast<double>(evalResults.pathList.size());
      if (evalCount <= 0.0) evalCount = 1.0;
      fastmath::logFastMathMetrics("minisim", std::cout, evalCount);
#endif

      // prepare results for next cycle
      evalResults.pathList.clear();
      evalResults.aircraftStateList.clear();
      evalResults.crashReasonList.clear();
      evalResults.scenarioList.clear();
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

  initializeSimGP();

  unsigned short port = std::atoi(argv[3]);
  boost::asio::io_context io_context;
  SimProcess sim_process(io_context, port);
  sim_process.run();

  return 0;
}

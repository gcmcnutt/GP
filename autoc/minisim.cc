/* test sim for aircraft */
#include <boost/asio.hpp>
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
#include <unistd.h>

#include "gp.h"
#include "minisim.h"
#include "autoc.h"
#include "gp_bytecode.h"

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
    workerPid = static_cast<int>(getpid());
    workerId = 0; // single-worker context uses index 0
  }

  void run() {
    while (true) { // TODO have reply have a controlled loop exit
      EvalResults evalResults;
      evalResults.workerPid = workerPid;
      evalResults.workerId = workerId;
      evalResults.workerEvalCounter = ++evalCounter;

      // ok what does main say to do
      EvalData evalData = receiveRPC<EvalData>(socket_);

      const uint64_t localGpHash = hashByteVector(evalData.gp);
      if (evalData.gpHash == 0) {
        evalData.gpHash = localGpHash;
      } else if (evalData.gpHash != localGpHash) {
        std::cerr << "[AUTOC_GP_HASH_MISMATCH] workerId=" << workerId
                  << " expected=0x" << std::hex << evalData.gpHash
                  << " got=0x" << localGpHash << std::dec
                  << " size=" << evalData.gp.size()
                  << std::endl;
      }

      ensureScenarioMetadata(evalData);

      // flip this back for return trip
      evalResults.gp = evalData.gp;
      evalResults.gpHash = localGpHash;
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

      // Always log GP program while debugging
      static const bool logGpStrings = true;
      if (logGpStrings) {
        if (isGPTreeData) {
          std::ostringstream oss;
          gp.printOn(oss);
          std::string treeStr = oss.str();
          const size_t maxLen = 800;
          if (treeStr.size() > maxLen) {
            treeStr.resize(maxLen);
            treeStr.append("...<truncated>");
          }
          std::cerr << "[AUTOC_GP_STRING] worker=" << workerId
                    << " eval=" << evalCounter
                    << " hash=0x" << std::hex << evalData.gpHash << std::dec
                    << " bytes=" << evalData.gp.size()
                    << " tree=\"" << treeStr << "\""
                    << std::endl;
        } else if (isBytecodeData) {
          std::cerr << "[AUTOC_GP_STRING] worker=" << workerId
                    << " eval=" << evalCounter
                    << " hash=0x" << std::hex << evalData.gpHash << std::dec
                    << " bytecode_bytes=" << evalData.gp.size()
                    << " (bytecode)"
                    << std::endl;
        }
      }
      
      // for each path, evaluate
      for (int i = 0; i < evalData.pathList.size(); i++) {
        std::vector<Path> path = evalData.pathList.at(i);

        // accumulate steps
        std::vector<AircraftState> aircraftStateSteps;

        // Fixed initial orientation and position for deterministic evaluation
        gp_quat aircraft_orientation;
        gp_vec3 initialPosition;
        {
          // Set orientation to be upright and level (180 deg yaw, 0 pitch, 0 roll)
          aircraft_orientation = gp_quat(Eigen::AngleAxis<gp_scalar>(static_cast<gp_scalar>(M_PI), gp_vec3::UnitZ())) *
            gp_quat(Eigen::AngleAxis<gp_scalar>(0, gp_vec3::UnitY())) *
            gp_quat(Eigen::AngleAxis<gp_scalar>(0, gp_vec3::UnitX()));

          // Set fixed starting position
          initialPosition = gp_vec3(static_cast<gp_scalar>(-2.19f), static_cast<gp_scalar>(5.49f), static_cast<gp_scalar>(-36.93f));
        }

        // compute initial velocity vector based on aircraft orientation (keeping minisim simple)
        gp_vec3 initial_velocity = aircraft_orientation * gp_vec3(SIM_INITIAL_VELOCITY, 0.0f, 0.0f);
        
        // reset sim state
        aircraftState = AircraftState{ 0, SIM_INITIAL_VELOCITY, initial_velocity, aircraft_orientation, initialPosition, 0.0f, 0.0f, SIM_INITIAL_THROTTLE, 0 };

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
          gp_scalar pre_roll = aircraftState.getRollCommand();
          gp_scalar pre_pitch = aircraftState.getPitchCommand();
          gp_scalar pre_throttle = aircraftState.getThrottleCommand();
          
          // run the controller (GP tree or bytecode interpreter) - BOTH NOW HAVE SAME BASELINE
          gp_scalar evaluation_result = 0.0f;
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
          gp_scalar distanceFromOrigin = std::sqrt(aircraftState.getPosition()[0] * aircraftState.getPosition()[0] +
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

      // prepare results for next cycle
      evalResults.pathList.clear();
      evalResults.aircraftStateList.clear();
      evalResults.crashReasonList.clear();
      evalResults.scenarioList.clear();
    }
  }

private:
  tcp::socket socket_;
  int workerPid = 0;
  int workerId = 0;
  int evalCounter = 0;
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

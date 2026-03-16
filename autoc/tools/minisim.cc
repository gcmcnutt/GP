/* minisim — lightweight simulator for NN controller evaluation */
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unistd.h>

#include "autoc/rpc/protocol.h"
#include "autoc/autoc.h"
#include "autoc/eval/sensor_math.h"
#include "autoc/nn/serialization.h"
#include "autoc/nn/evaluator.h"

using namespace std;

// Global aircraft state (declared extern in autoc.h)
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
}

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

class SimProcess {
public:
  SimProcess(int id, unsigned short port) {
    socket_.connect("localhost", port);
    workerPid = static_cast<int>(getpid());
    workerId = id;
  }

  void run() {
    while (true) {
      EvalResults evalResults;
      evalResults.workerPid = workerPid;
      evalResults.workerId = workerId;
      evalResults.workerEvalCounter = ++evalCounter;

      EvalData evalData = receiveRPC<EvalData>(socket_);

      const uint64_t localGpHash = hashByteVector(evalData.gp);
      if (evalData.gpHash == 0) {
        evalData.gpHash = localGpHash;
      } else if (evalData.gpHash != localGpHash) {
        std::cerr << "[HASH_MISMATCH] workerId=" << workerId
                  << " expected=0x" << std::hex << evalData.gpHash
                  << " got=0x" << localGpHash << std::dec
                  << " size=" << evalData.gp.size()
                  << std::endl;
      }

      ensureScenarioMetadata(evalData);

      evalResults.gp = evalData.gp;
      evalResults.gpHash = localGpHash;
      if (!evalData.scenarioList.empty()) {
        evalResults.scenario = evalData.scenarioList.front();
      } else {
        evalResults.scenario = evalData.scenario;
      }
      evalResults.scenarioList.clear();
      evalResults.scenarioList.reserve(evalData.scenarioList.size());

      // Deserialize NN genome
      NNGenome nnGenome;
      if (!nn_detect_format(reinterpret_cast<const uint8_t*>(evalData.gp.data()), evalData.gp.size())) {
        std::cerr << "[MINISIM] Payload missing NN01 magic bytes, size=" << evalData.gp.size() << std::endl;
        continue;
      }
      if (!nn_deserialize(reinterpret_cast<const uint8_t*>(evalData.gp.data()), evalData.gp.size(), nnGenome)) {
        std::cerr << "[MINISIM] Failed to deserialize NN genome" << std::endl;
        continue;
      }

      // Log topology once on first eval
      if (evalCounter == 1) {
        std::ostringstream topo;
        for (size_t i = 0; i < nnGenome.topology.size(); i++) {
          if (i > 0) topo << "x";
          topo << nnGenome.topology[i];
        }
        std::cerr << "[MINISIM] worker=" << workerId
                  << " topology=" << topo.str()
                  << " weights=" << nnGenome.weights.size()
                  << std::endl;
      }

      // Evaluate each path
      for (int i = 0; i < static_cast<int>(evalData.pathList.size()); i++) {
        std::vector<Path> path = evalData.pathList.at(i);
        std::vector<AircraftState> aircraftStateSteps;

        // Fixed initial orientation and position for deterministic evaluation
        gp_quat aircraft_orientation = gp_quat(Eigen::AngleAxis<gp_scalar>(static_cast<gp_scalar>(M_PI), gp_vec3::UnitZ())) *
          gp_quat(Eigen::AngleAxis<gp_scalar>(0, gp_vec3::UnitY())) *
          gp_quat(Eigen::AngleAxis<gp_scalar>(0, gp_vec3::UnitX()));
        gp_vec3 initialPosition(static_cast<gp_scalar>(-2.19f), static_cast<gp_scalar>(5.49f), static_cast<gp_scalar>(-36.93f));
        gp_vec3 initial_velocity = aircraft_orientation * gp_vec3(SIM_INITIAL_VELOCITY, 0.0f, 0.0f);

        aircraftState = AircraftState{ 0, SIM_INITIAL_VELOCITY, initial_velocity, aircraft_orientation, initialPosition, 0.0f, 0.0f, SIM_INITIAL_THROTTLE, 0 };
        aircraftState.clearHistory();

        aircraftStateSteps.push_back(aircraftState);

        unsigned long int duration_msec = 0;
        CrashReason crashReason = CrashReason::None;

        while (crashReason == CrashReason::None) {

          // Capture temporal history before NN evaluation
          {
            VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
            gp_scalar dPhi = executeGetDPhi(pathProvider, aircraftState, 0.0f);
            gp_scalar dTheta = executeGetDTheta(pathProvider, aircraftState, 0.0f);
            gp_vec3 targetPos = getInterpolatedTargetPosition(
                pathProvider, static_cast<int32_t>(aircraftState.getSimTimeMsec()), 0.0f);
            gp_scalar distance = (targetPos - aircraftState.getPosition()).norm();
            aircraftState.recordErrorHistory(dPhi, dTheta, distance, duration_msec);
          }

          // Run NN controller
          {
            VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
            NNControllerBackend nnBackend(nnGenome);
            nnBackend.evaluate(aircraftState, pathProvider);
          }

          // Advance aircraft state
          aircraftState.minisimAdvanceState(SIM_TIME_STEP_MSEC);
          duration_msec += SIM_TIME_STEP_MSEC;
          aircraftState.setSimTimeMsec(duration_msec);

          // Crash detection
          gp_scalar distanceFromOrigin = std::sqrt(aircraftState.getPosition()[0] * aircraftState.getPosition()[0] +
            aircraftState.getPosition()[1] * aircraftState.getPosition()[1]);
          if (aircraftState.getPosition()[2] < (SIM_MAX_ELEVATION) ||
            aircraftState.getPosition()[2] > (SIM_MIN_ELEVATION) ||
            distanceFromOrigin > SIM_PATH_RADIUS_LIMIT) {
            crashReason = CrashReason::Eval;
          }

          // Advance path index
          while (aircraftState.getThisPathIndex() < static_cast<int>(path.size()) - 2 &&
                 (path.at(aircraftState.getThisPathIndex()).simTimeMsec < static_cast<int32_t>(duration_msec))) {
            aircraftState.setThisPathIndex(aircraftState.getThisPathIndex() + 1);
          }

          aircraftStateSteps.push_back(aircraftState);

          if (duration_msec >= SIM_TOTAL_TIME_MSEC) {
            crashReason = CrashReason::Time;
          }
          if (aircraftState.getThisPathIndex() >= static_cast<int>(path.size()) - 2) {
            crashReason = CrashReason::Distance;
          }
        }

        {
          ScenarioMetadata meta = scenarioForPathIndex(evalData, static_cast<size_t>(i));
          evalResults.scenarioList.push_back(meta);
        }
        evalResults.pathList.push_back(path);
        evalResults.aircraftStateList.push_back(aircraftStateSteps);
        evalResults.crashReasonList.push_back(crashReason);
      }

      sendRPC(socket_, evalResults);

      evalResults.pathList.clear();
      evalResults.aircraftStateList.clear();
      evalResults.crashReasonList.clear();
      evalResults.scenarioList.clear();
    }
  }

private:
  TcpSocket socket_;
  int workerPid = 0;
  int workerId = 0;
  int evalCounter = 0;
};

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: minisim <base> <id> <port>" << std::endl;
    return 1;
  }

  int id = std::atoi(argv[2]);
  unsigned short port = std::atoi(argv[3]);
  SimProcess sim_process(id, port);
  sim_process.run();

  return 0;
}

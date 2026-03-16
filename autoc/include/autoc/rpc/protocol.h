/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <vector>
#include <iostream>
#include <cstdint>
#include <arpa/inet.h>

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include <cstdio>

// Cross-platform serialization notes:
// - cereal binary archive serializes primitives in native format
// - Works across x86/ARM if both are little-endian (typical for modern systems)
// - IEEE 754 float representation is standard
// - Custom Eigen serialization (gp_vec3, gp_quat) uses element-wise saves (portable)
// - Current implementation: deterministic within same architecture

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "autoc/util/socket_wrapper.h"
#include "autoc/eval/aircraft_state.h"

// FNV-1a hash for serialized GP programs/bytecode (matches dtest tracker)
inline uint64_t hashByteVector(const std::vector<char>& data) {
  const uint8_t* bytes = reinterpret_cast<const uint8_t*>(data.data());
  uint64_t hash = 0xcbf29ce484222325ULL;
  for (size_t i = 0; i < data.size(); ++i) {
    hash ^= bytes[i];
    hash *= 0x100000001b3ULL;
  }
  return hash;
}

// Hash a string using FNV-1a (deterministic, only depends on string content)
inline uint64_t hashString(const std::string& str) {
  const uint8_t* bytes = reinterpret_cast<const uint8_t*>(str.data());
  uint64_t hash = 0xcbf29ce484222325ULL;
  for (size_t i = 0; i < str.size(); ++i) {
    hash ^= bytes[i];
    hash *= 0x100000001b3ULL;
  }
  return hash;
}

namespace cereal {
    template<class Archive>
    void save(Archive& ar, const gp_vec3& v) {
        ar(v[0], v[1], v[2]);
    }
    template<class Archive>
    void load(Archive& ar, gp_vec3& v) {
        ar(v[0], v[1], v[2]);
    }
    template<class Archive>
    void save(Archive& ar, const gp_quat& q) {
        gp_scalar w = q.w(), x = q.x(), y = q.y(), z = q.z();
        ar(w, x, y, z);
    }
    template<class Archive>
    void load(Archive& ar, gp_quat& q) {
        gp_scalar w, x, y, z;
        ar(w, x, y, z);
        q = gp_quat(w, x, y, z);
    }
}


/*
 * here we send our requested paths to the sims
 */
struct ScenarioMetadata {
  int pathVariantIndex = -1;   // -1 = unset/aggregated
  int windVariantIndex = -1;   // -1 = unset/aggregated
  unsigned int windSeed = 0;
  uint64_t scenarioSequence = 0;
  uint64_t bakeoffSequence = 0;
  bool enableDeterministicLogging = false;

  // VARIATIONS1: Entry and wind direction offsets (computed by autoc, applied by crrcsim)
  // All angles in radians, speed as multiplier
  double entryHeadingOffset = 0.0;   // radians, offset from path tangent
  double entryRollOffset = 0.0;      // radians, initial roll attitude
  double entryPitchOffset = 0.0;     // radians, initial pitch attitude
  double entrySpeedFactor = 1.0;     // multiplier on reference speed
  double windDirectionOffset = 0.0;  // radians, offset from base wind direction

  // Entry position offsets (see specs/005-entry-fitness-ramp)
  double entryNorthOffset = 0.0;     // meters, NED North
  double entryEastOffset = 0.0;      // meters, NED East
  double entryAltOffset = 0.0;       // meters, NED Down (negative=up)

  template<class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {
    ar(pathVariantIndex, windVariantIndex, windSeed, scenarioSequence,
       bakeoffSequence, enableDeterministicLogging, entryHeadingOffset,
       entryRollOffset, entryPitchOffset, entrySpeedFactor,
       windDirectionOffset, entryNorthOffset, entryEastOffset, entryAltOffset);
  }
};
CEREAL_CLASS_VERSION(ScenarioMetadata, 1)

// Controller type tag for RPC wire format
enum class ControllerType : int {
  NEURAL_NET = 2   // NNGenome weight vector
};
// cereal doesn't serialize enum class directly - use int wrapper

struct EvalData {
  std::vector<char> gp;
  uint64_t gpHash = 0;  // FNV-1a hash of payload for verification
  bool isEliteReeval = false;
  ControllerType controllerType = ControllerType::NEURAL_NET;
  std::vector<std::vector<Path>> pathList;
  ScenarioMetadata scenario;
  std::vector<ScenarioMetadata> scenarioList;

  template<class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {
    ar(gp, gpHash, isEliteReeval);
    int ct = static_cast<int>(controllerType);
    ar(ct);
    controllerType = static_cast<ControllerType>(ct);
    ar(pathList, scenario, scenarioList);
  }

  void sanitizePaths() {
    for (auto& pathGroup : pathList) {
      for (auto& path : pathGroup) {
        path.sanitize();
      }
    }
  }
};
CEREAL_CLASS_VERSION(EvalData, 1)

enum class CrashReason {
  None,
  Boot,
  Sim,
  Eval,
  Time,
  Distance,
};

std::string crashReasonToString(CrashReason type);

struct DebugSample {
  int pathIndex = -1;
  int stepIndex = -1;
  gp_scalar simTimeMsec = static_cast<gp_scalar>(0.0f);
  gp_scalar dtSec = static_cast<gp_scalar>(0.0f);
  gp_scalar simSteps = static_cast<gp_scalar>(0.0f);
  gp_vec3 velRelGround = gp_vec3::Zero();   // FDM ground-relative velocity (m/s)
  gp_vec3 velRelAirmass = gp_vec3::Zero();  // FDM airmass-relative velocity (m/s)
  gp_vec3 position = gp_vec3::Zero();
  gp_vec3 velocity = gp_vec3::Zero();
  gp_vec3 acceleration = gp_vec3::Zero();
  gp_vec3 accelPast = gp_vec3::Zero();
  gp_vec3 angularRates = gp_vec3::Zero();  // body P/Q/R rad/s
  gp_vec3 angularAccelPast = gp_vec3::Zero();
  gp_quat quatDotPast = gp_quat::Identity(); // store e_dot_past as quat-like container
  gp_vec3 windBody = gp_vec3::Zero();      // body-frame wind (m/s)
  gp_quat orientation = gp_quat::Identity();
  gp_scalar pitchCommand = static_cast<gp_scalar>(0.0f);
  gp_scalar rollCommand = static_cast<gp_scalar>(0.0f);
  gp_scalar throttleCommand = static_cast<gp_scalar>(0.0f);
  gp_scalar elevatorSim = static_cast<gp_scalar>(0.0f);
  gp_scalar aileronSim = static_cast<gp_scalar>(0.0f);
  gp_scalar throttleSim = static_cast<gp_scalar>(0.0f);
  gp_scalar massKg = static_cast<gp_scalar>(0.0f);
  gp_scalar density = static_cast<gp_scalar>(0.0f);
  gp_scalar gravity = static_cast<gp_scalar>(0.0f);
  gp_scalar alpha = static_cast<gp_scalar>(0.0f);
  gp_scalar beta = static_cast<gp_scalar>(0.0f);
  gp_scalar vRelWind = static_cast<gp_scalar>(0.0f);
  gp_vec3 localAirmass = gp_vec3::Zero();   // local-frame air velocity
  gp_vec3 gustBody = gp_vec3::Zero();       // body-frame gust applied
  gp_vec3 forceBody = gp_vec3::Zero();      // last total force (body)
  gp_vec3 momentBody = gp_vec3::Zero();     // last moment at CG (body)
  gp_vec3 vLocal = gp_vec3::Zero();         // local-frame velocity (ft/s -> m/s)
  gp_vec3 vLocalDot = gp_vec3::Zero();      // local-frame acceleration (ft/s^2 -> m/s^2)
  gp_vec3 omegaBody = gp_vec3::Zero();      // body rates (rad/s)
  gp_vec3 omegaDotBody = gp_vec3::Zero();   // body angular accel (rad/s^2)
  gp_scalar latGeoc = static_cast<gp_scalar>(0.0f);
  gp_scalar lonGeoc = static_cast<gp_scalar>(0.0f);
  gp_scalar radiusToVehicle = static_cast<gp_scalar>(0.0f); // meters
  gp_scalar latDotPast = static_cast<gp_scalar>(0.0f);
  gp_scalar lonDotPast = static_cast<gp_scalar>(0.0f);
  gp_scalar radiusDotPast = static_cast<gp_scalar>(0.0f);
  uint32_t rngState16 = 0;
  uint32_t rngState32 = 0;

  template<class Archive>
  void serialize(Archive& ar, const std::uint32_t /*version*/) {
    ar(pathIndex, stepIndex, simTimeMsec, dtSec, simSteps,
       velRelGround, velRelAirmass, position, velocity, acceleration,
       accelPast, angularRates, angularAccelPast, quatDotPast, windBody,
       orientation, pitchCommand, rollCommand, throttleCommand,
       elevatorSim, aileronSim, throttleSim, massKg, density, gravity,
       alpha, beta, vRelWind, localAirmass, gustBody, forceBody,
       momentBody, vLocal, vLocalDot, omegaBody, omegaDotBody,
       latGeoc, lonGeoc, radiusToVehicle, latDotPast, lonDotPast,
       radiusDotPast, rngState16, rngState32);
  }
};

struct EvalResults {
  std::vector<char> gp;
  uint64_t gpHash = 0;  // FNV-1a hash of gp buffer for verification
  std::vector<CrashReason> crashReasonList;
  std::vector<std::vector<Path>> pathList;
  std::vector<std::vector<AircraftState>> aircraftStateList;
  ScenarioMetadata scenario;
  std::vector<ScenarioMetadata> scenarioList;
  std::vector<std::vector<DebugSample>> debugSamples;  // Debug snapshots per path (only populated for elite reeval)
  std::vector<std::vector<PhysicsTraceEntry>> physicsTrace;  // Full physics state trace per path (only populated for elite reeval)
  int workerId = -1;
  int workerPid = 0;
  int workerEvalCounter = 0;  // Incremented per evaluation on the worker

  template<class Archive>
  void serialize(Archive& ar, const std::uint32_t version) {
    ar(gp, gpHash, crashReasonList, pathList, aircraftStateList,
       scenario, scenarioList, debugSamples, physicsTrace,
       workerId, workerPid, workerEvalCounter);
  }

  void clear() {
    gp.clear();
    gpHash = 0;
    crashReasonList.clear();
    pathList.clear();
    aircraftStateList.clear();
    scenario = ScenarioMetadata();
    scenarioList.clear();
    debugSamples.clear();
    physicsTrace.clear();
    workerId = -1;
    workerPid = 0;
    workerEvalCounter = 0;
  }

  void dump(std::ostream& os) {
    char buf[512];
    snprintf(buf, sizeof(buf), "EvalResults: crash[%zu] paths[%zu] states[%zu]\n Paths:\n",
      crashReasonList.size(), pathList.size(), aircraftStateList.size());
    os << buf;

    snprintf(buf, sizeof(buf), "Scenario: pathVariant=%d windVariant=%d windSeed=%u seq=%llu bake=%llu\n",
      scenario.pathVariantIndex, scenario.windVariantIndex, scenario.windSeed,
      static_cast<unsigned long long>(scenario.scenarioSequence),
      static_cast<unsigned long long>(scenario.bakeoffSequence));
    os << buf;
    for (size_t i = 0; i < scenarioList.size(); ++i) {
      snprintf(buf, sizeof(buf), "  Scenario[%zu]: pathVariant=%d windVariant=%d windSeed=%u seq=%llu bake=%llu\n",
        i, scenarioList[i].pathVariantIndex, scenarioList[i].windVariantIndex, scenarioList[i].windSeed,
        static_cast<unsigned long long>(scenarioList[i].scenarioSequence),
        static_cast<unsigned long long>(scenarioList[i].bakeoffSequence));
      os << buf;
    }

    for (size_t i = 0; i < crashReasonList.size(); i++) {
      snprintf(buf, sizeof(buf), "  Crash %3zu: %s\n", i, crashReasonToString(crashReasonList.at(i)).c_str());
      os << buf;
    }

    for (size_t i = 0; i < pathList.size(); i++) {
      for (size_t j = 0; j < pathList.at(i).size(); j++) {
        Path path = pathList.at(i).at(j);
        snprintf(buf, sizeof(buf), "  Path %3zu:%3zu: start[%8.2f %8.2f %8.2f] orient[%8.2f %8.2f %8.2f] dist[%8.2f] rad[%8.2f]\n",
          i, j,
          path.start[0], path.start[1], path.start[2],
          path.orientation.x(), path.orientation.y(), path.orientation.z(),
          path.distanceFromStart, path.radiansFromStart);
        os << buf;
      }
    }
    os << " Aircraft States:\n";
    for (size_t i = 0; i < aircraftStateList.size(); i++) {
      for (size_t j = 0; j < aircraftStateList.at(i).size(); j++) {
        AircraftState aircraftState = aircraftStateList.at(i).at(j);
        gp_quat orientQuatF = aircraftState.getOrientation();
        if (!std::isnan(orientQuatF.norm()) && std::abs(orientQuatF.norm() - 1.0f) > 1e-6f) {
          orientQuatF.normalize();
        }

        Eigen::Matrix<gp_scalar, 3, 1> euler = orientQuatF.toRotationMatrix().eulerAngles(2, 1, 0);
        gp_vec3 eulerWrapped;
        for (int axis = 0; axis < 3; ++axis) {
          eulerWrapped[axis] = std::atan2(std::sin(euler[axis]), std::cos(euler[axis]));
        }
        snprintf(buf, sizeof(buf),
          "  Path %3zu:%3zu: Time %5ld Index %3d: pos[%8.2f %8.2f %8.2f] orientRPY[%8.2f %8.2f %8.2f] quat[%+7.4f %+7.4f %+7.4f %+7.4f] vel[%8.2f] pitch[%5.2f] roll[%5.2f] throttle[%5.2f]\n",
          i, j,
          aircraftState.getSimTimeMsec(),
          aircraftState.getThisPathIndex(),
          aircraftState.getPosition()[0], aircraftState.getPosition()[1], aircraftState.getPosition()[2],
          eulerWrapped.x(), eulerWrapped.y(), eulerWrapped.z(),
          orientQuatF.w(), orientQuatF.x(), orientQuatF.y(), orientQuatF.z(),
          aircraftState.getRelVel(),
          aircraftState.getPitchCommand(),
          aircraftState.getRollCommand(),
          aircraftState.getThrottleCommand());
        os << buf;
      }
    }
  }
};
CEREAL_CLASS_VERSION(EvalResults, 1)

struct WorkerContext {
  std::unique_ptr<TcpSocket> socket;
  pid_t childPid = 0;
  EvalResults evalResults;
  int workerId = -1;
};

/*
 * generic RPC wrappers
 */
template<typename T>
void sendRPC(TcpSocket& socket, const T& data) {
  std::ostringstream os(std::ios::binary);
  {
    cereal::BinaryOutputArchive archive(os);
    archive(data);
  }
  std::string outbound_data = os.str();
  uint32_t size = htonl(static_cast<uint32_t>(outbound_data.size()));
  socket.write(&size, sizeof(size));
  socket.write(outbound_data.data(), outbound_data.size());
}

template<typename T>
T receiveRPC(TcpSocket& socket) {
  uint32_t size;
  socket.read(&size, sizeof(size));
  size = ntohl(size);
  std::vector<char> buffer(size);
  socket.read(buffer.data(), buffer.size());
  std::istringstream is(std::string(buffer.begin(), buffer.end()), std::ios::binary);
  cereal::BinaryInputArchive archive(is);
  T data;
  archive(data);
  return data;
}


#endif

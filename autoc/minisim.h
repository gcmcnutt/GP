/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <vector>
#include <iostream>
#include <cstdint>
#include <arpa/inet.h>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/asio.hpp>
#include <boost/process.hpp>
#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

// Cross-platform serialization notes:
// - Boost binary_archive serializes primitives (int, float) in native format
// - Works across x86/ARM if both are little-endian (typical for modern systems)
// - IEEE 754 float representation is standard
// - Custom Eigen serialization (gp_vec3, gp_quat) uses element-wise saves (portable)
// - For true cross-endian support, consider eos::portable_oarchive
// - Current implementation: deterministic within same architecture

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using boost::asio::ip::tcp;
using boost::format;

#include "aircraft_state.h"

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

namespace boost {
  namespace serialization {

    // Serialization for float-based Eigen vector/quaternion used by GP eval
    // Serialize elements individually for cross-platform portability (avoiding
    // binary serialization of raw float arrays which is not portable)
    //
    // Note: We use split save/load to ensure proper portable serialization.
    // The save always writes in portable format, and load handles both old
    // (binary array) and new (element-wise) formats for backward compatibility.
    template<class Archive>
    void save(Archive& ar, const gp_vec3& v, const unsigned int version)
    {
      ar & v[0] & v[1] & v[2];
    }

    template<class Archive>
    void load(Archive& ar, gp_vec3& v, const unsigned int version)
    {
      ar & v[0] & v[1] & v[2];
    }

    // Serialization for Eigen::Quaternion<float>
    // Serialize coefficients individually for cross-platform portability
    template<class Archive>
    void save(Archive& ar, const gp_quat& q, const unsigned int version)
    {
      gp_scalar w = q.w();
      gp_scalar x = q.x();
      gp_scalar y = q.y();
      gp_scalar z = q.z();
      ar & w & x & y & z;
    }

    template<class Archive>
    void load(Archive& ar, gp_quat& q, const unsigned int version)
    {
      gp_scalar w, x, y, z;
      ar & w & x & y & z;
      q = gp_quat(w, x, y, z);
    }

  } // namespace serialization
} // namespace boost

// These macros tell boost to use the save/load functions we just defined
BOOST_SERIALIZATION_SPLIT_FREE(gp_vec3)
BOOST_SERIALIZATION_SPLIT_FREE(gp_quat)


/*
 * here we send our requested paths to the sims
 */
struct ScenarioMetadata {
  int pathVariantIndex = 0;
  int windVariantIndex = 0;
  unsigned int windSeed = 0;
  uint64_t scenarioSequence = 0;
  uint64_t bakeoffSequence = 0;
  bool enableDeterministicLogging = false;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& pathVariantIndex;
    ar& windVariantIndex;
    ar& windSeed;
    if (version > 2) {
      ar& scenarioSequence;
      ar& bakeoffSequence;
    }
    if (version > 3) {
      ar& enableDeterministicLogging;
    } else if (version > 1) {
      ar& scenarioSequence;
      if (Archive::is_loading::value) {
        bakeoffSequence = 0;
      }
    } else if (Archive::is_loading::value) {
      scenarioSequence = 0;
      bakeoffSequence = 0;
    }
  }
};
BOOST_CLASS_VERSION(ScenarioMetadata, 4)

struct EvalData {
  std::vector<char> gp;
  uint64_t gpHash = 0;  // FNV-1a hash of gp buffer for verification
  bool isEliteReeval = false;  // Enable detailed logging for elite re-evaluations
  std::vector<std::vector<Path>> pathList;
  ScenarioMetadata scenario;
  std::vector<ScenarioMetadata> scenarioList;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& gp;
    if (version > 3) {
      ar& gpHash;
      if (version > 4) {
        ar& isEliteReeval;
      } else if (Archive::is_loading::value) {
        isEliteReeval = false;
      }
    } else if (Archive::is_loading::value) {
      gpHash = 0;
      isEliteReeval = false;
    }
    ar& pathList;
    if (version > 1) {
      ar& scenario;
      if (version > 2) {
        ar& scenarioList;
      } else if (Archive::is_loading::value) {
        scenarioList.assign(pathList.size(), scenario);
        for (size_t idx = 0; idx < scenarioList.size(); ++idx) {
          if (scenarioList[idx].pathVariantIndex < 0) {
            scenarioList[idx].pathVariantIndex = static_cast<int>(idx);
          }
        }
      }
    } else if (Archive::is_loading::value) {
      scenario = ScenarioMetadata();
      scenarioList.assign(pathList.size(), scenario);
      for (size_t idx = 0; idx < scenarioList.size(); ++idx) {
        scenarioList[idx].pathVariantIndex = static_cast<int>(idx);
      }
    }
  }

  void sanitizePaths() {
    for (auto& pathGroup : pathList) {
      for (auto& path : pathGroup) {
        path.sanitize();
      }
    }
  }
};
BOOST_CLASS_VERSION(EvalData, 5)  // Bumped for isEliteReeval field

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
  void serialize(Archive& ar, const unsigned int /*version*/) {
    ar& pathIndex;
    ar& stepIndex;
    ar& simTimeMsec;
    ar& dtSec;
    ar& simSteps;
    ar& velRelGround;
    ar& velRelAirmass;
    ar& position;
    ar& velocity;
    ar& acceleration;
    ar& accelPast;
    ar& angularRates;
    ar& angularAccelPast;
    ar& quatDotPast;
    ar& windBody;
    ar& orientation;
    ar& pitchCommand;
    ar& rollCommand;
    ar& throttleCommand;
    ar& elevatorSim;
    ar& aileronSim;
    ar& throttleSim;
    ar& massKg;
    ar& density;
    ar& gravity;
    ar& alpha;
    ar& beta;
    ar& vRelWind;
    ar& localAirmass;
    ar& gustBody;
    ar& forceBody;
    ar& momentBody;
    ar& vLocal;
    ar& vLocalDot;
    ar& omegaBody;
    ar& omegaDotBody;
    ar& latGeoc;
    ar& lonGeoc;
    ar& radiusToVehicle;
    ar& latDotPast;
    ar& lonDotPast;
    ar& radiusDotPast;
    ar& rngState16;
    ar& rngState32;
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

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& gp;
    if (version > 7) {
      ar& gpHash;
    } else if (Archive::is_loading::value) {
      gpHash = 0;
    }
    ar& crashReasonList;
    ar& pathList;
    ar& aircraftStateList;
    if (version > 1) {
      ar& scenario;
      if (version > 2) {
        ar& scenarioList;
        if (version > 3) {
          ar& debugSamples;
          if (version > 8) {
            ar& physicsTrace;
          }
          ar& workerId;
          ar& workerPid;
          ar& workerEvalCounter;
        }
      } else if (Archive::is_loading::value) {
        scenarioList.assign(pathList.size(), scenario);
        for (size_t idx = 0; idx < scenarioList.size(); ++idx) {
          scenarioList[idx].pathVariantIndex = static_cast<int>(idx);
        }
      }
    } else if (Archive::is_loading::value) {
      scenario = ScenarioMetadata();
      scenarioList.assign(pathList.size(), scenario);
      for (size_t idx = 0; idx < scenarioList.size(); ++idx) {
        scenarioList[idx].pathVariantIndex = static_cast<int>(idx);
      }
      debugSamples.clear();
      physicsTrace.clear();
      workerId = -1;
      workerPid = 0;
      workerEvalCounter = 0;
    }
  }

  void dump(std::ostream& os) {
    os << format("EvalResults: crash[%d] paths[%d] states[%d]\n Paths:\n")
      % crashReasonList.size() % pathList.size() % aircraftStateList.size();

    os << format("GP: %s\n") % "TODO";
    os << format("Scenario: pathVariant=%d windVariant=%d windSeed=%u seq=%llu bake=%llu\n")
      % scenario.pathVariantIndex % scenario.windVariantIndex % scenario.windSeed
      % static_cast<unsigned long long>(scenario.scenarioSequence)
      % static_cast<unsigned long long>(scenario.bakeoffSequence);
    for (size_t i = 0; i < scenarioList.size(); ++i) {
      os << format("  Scenario[%zu]: pathVariant=%d windVariant=%d windSeed=%u seq=%llu bake=%llu\n")
        % i % scenarioList[i].pathVariantIndex % scenarioList[i].windVariantIndex % scenarioList[i].windSeed
        % static_cast<unsigned long long>(scenarioList[i].scenarioSequence)
        % static_cast<unsigned long long>(scenarioList[i].bakeoffSequence);
    }

    for (int i = 0; i < crashReasonList.size(); i++) {
      os << format("  Crash %3d: %s\n") % i % crashReasonToString(crashReasonList.at(i));
    }

    for (int i = 0; i < pathList.size(); i++) {
      for (int j = 0; j < pathList.at(i).size(); j++) {
        Path path = pathList.at(i).at(j);
        os << format("  Path %3d:%3d: start[%8.2f %8.2f %8.2f] orient[%8.2f %8.2f %8.2f] dist[%8.2f] rad[%8.2f]\n")
          % i
          % j
          % path.start[0] % path.start[1] % path.start[2]
          % path.orientation.x() % path.orientation.y() % path.orientation.z()
          % path.distanceFromStart % path.radiansFromStart;
      }
    }
    os << " Aircraft States:\n";
    for (int i = 0; i < aircraftStateList.size(); i++) {
      for (int j = 0; j < aircraftStateList.at(i).size(); j++) {
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
        // TODO distance from start
        os << format("  Path %3d:%3d: Time %5d Index %3d: pos[%8.2f %8.2f %8.2f] orientRPY[%8.2f %8.2f %8.2f] quat[%+7.4f %+7.4f %+7.4f %+7.4f] vel[%8.2f] pitch[%5.2f] roll[%5.2f] throttle[%5.2f]\n")
          % i
          % j
          % aircraftState.getSimTimeMsec()
          % aircraftState.getThisPathIndex()
          % aircraftState.getPosition()[0] % aircraftState.getPosition()[1] % aircraftState.getPosition()[2]
          % eulerWrapped.x() % eulerWrapped.y() % eulerWrapped.z()
          % orientQuatF.w() % orientQuatF.x() % orientQuatF.y() % orientQuatF.z()
          % aircraftState.getRelVel()
          % aircraftState.getPitchCommand()
          % aircraftState.getRollCommand()
          % aircraftState.getThrottleCommand();
      }
    }
  }
};
BOOST_CLASS_VERSION(EvalResults, 9)

struct WorkerContext {
  std::unique_ptr<boost::asio::ip::tcp::socket> socket;
  boost::process::child child_process;
  EvalResults evalResults;
  int workerId = -1;  // set by thread pool
};

/*
 * generic RPC wrappers
 */
template<typename T>
void sendRPC(tcp::socket& socket, const T& data) {
  std::ostringstream archive_stream(std::ios::binary);
  boost::archive::binary_oarchive archive(archive_stream);
  archive << data;
  std::string outbound_data = archive_stream.str();

  // Send size header (4 bytes, network byte order)
  uint32_t size = htonl(static_cast<uint32_t>(outbound_data.size()));
  boost::asio::write(socket, boost::asio::buffer(&size, sizeof(size)));

  // Send binary data
  boost::asio::write(socket, boost::asio::buffer(outbound_data));
};

template<typename T>
T receiveRPC(tcp::socket& socket) {
  // Read size header (4 bytes, network byte order)
  uint32_t size;
  boost::asio::read(socket, boost::asio::buffer(&size, sizeof(size)));
  size = ntohl(size);

  // Read binary data
  std::vector<char> buffer(size);
  boost::asio::read(socket, boost::asio::buffer(buffer));

  std::istringstream archive_stream(std::string(buffer.begin(), buffer.end()), std::ios::binary);
  boost::archive::binary_iarchive archive(archive_stream);
  T data;
  archive >> data;
  return data;
};


#endif

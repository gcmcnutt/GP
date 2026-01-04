/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <vector>
#include <iostream>
#include <cstdint>

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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>

using boost::asio::ip::tcp;
using boost::format;

#include "aircraft_state.h"

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
  std::vector<std::vector<Path>> pathList;
  ScenarioMetadata scenario;
  std::vector<ScenarioMetadata> scenarioList;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& gp;
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
BOOST_CLASS_VERSION(EvalData, 3)

enum class CrashReason {
  None,
  Boot,
  Sim,
  Eval,
  Time,
  Distance,
};

std::string crashReasonToString(CrashReason type);

struct EvalResults {
  std::vector<char> gp;
  std::vector<CrashReason> crashReasonList;
  std::vector<std::vector<Path>> pathList;
  std::vector<std::vector<AircraftState>> aircraftStateList;
  ScenarioMetadata scenario;
  std::vector<ScenarioMetadata> scenarioList;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& gp;
    ar& crashReasonList;
    ar& pathList;
    ar& aircraftStateList;
    if (version > 1) {
      ar& scenario;
      if (version > 2) {
        ar& scenarioList;
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
BOOST_CLASS_VERSION(EvalResults, 3)

struct WorkerContext {
  std::unique_ptr<boost::asio::ip::tcp::socket> socket;
  boost::process::child child_process;
  EvalResults evalResults;
};

/*
 * generic RPC wrappers
 */
template<typename T>
void sendRPC(tcp::socket& socket, const T& data) {
  std::ostringstream archive_stream;
  boost::archive::text_oarchive archive(archive_stream);
  archive << data;
  std::string outbound_data = archive_stream.str() + "\n";
  boost::asio::write(socket, boost::asio::buffer(outbound_data));
};

template<typename T>
T receiveRPC(tcp::socket& socket) {
  boost::asio::streambuf buf;
  boost::asio::read_until(socket, buf, '\n');
  std::string archive_data((std::istreambuf_iterator<char>(&buf)),
    std::istreambuf_iterator<char>());
  std::istringstream archive_stream(archive_data);
  boost::archive::text_iarchive archive(archive_stream);
  T data;
  archive >> data;
  return data;
};


#endif

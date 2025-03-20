/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <vector>
#include <iostream>

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

#include <Eigen/Dense>
#include <Eigen/Geometry>

using boost::asio::ip::tcp;
using boost::format;

#include "aircraft_state.h"

namespace boost {
  namespace serialization {

    // Serialization for Eigen::Vector3d
    template<class Archive>
    void serialize(Archive& ar, Eigen::Vector3d& v, const unsigned int version)
    {
      ar& boost::serialization::make_array(v.data(), 3);
    }

    // Serialization for Eigen::Quaterniond
    template<class Archive>
    void save(Archive& ar, const Eigen::Quaterniond& q, const unsigned int version)
    {
      ar& boost::serialization::make_array(q.coeffs().data(), 4);
    }

    template<class Archive>
    void load(Archive& ar, Eigen::Quaterniond& q, const unsigned int version)
    {
      Eigen::Vector4d coeffs;
      ar& boost::serialization::make_array(coeffs.data(), 4);
      q = Eigen::Quaterniond(coeffs[3], coeffs[0], coeffs[1], coeffs[2]);
    }

  } // namespace serialization
} // namespace boost

// This macro tells boost to use the save/load functions we just defined for Quaterniond
BOOST_SERIALIZATION_SPLIT_FREE(Eigen::Quaterniond)


/*
 * here we send our requested paths to the sims
 */
  struct EvalData {
  std::vector<char> gp;
  std::vector<std::vector<Path>> pathList;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& gp;
    ar& pathList;
  }
};
BOOST_CLASS_VERSION(EvalData, 1)

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

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& gp;
    ar& crashReasonList;
    ar& pathList;
    ar& aircraftStateList;
  }

  void dump(std::ostream& os) {
    os << format("EvalResults: crash[%d] paths[%d] states[%d]\n Paths:\n")
      % crashReasonList.size() % pathList.size() % aircraftStateList.size();

    os << format("GP: %s\n") % "TODO";

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
        Eigen::Vector3d euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
        // TODO distance from start
        os << format("  Path %3d:%3d: Time %5d Index %3d: pos[%8.2f %8.2f %8.2f] orient[%8.2f %8.2f %8.2f] vel[%8.2f] pitch[%5.2f] roll[%5.2f] throttle[%5.2f]\n")
          % i
          % j
          % aircraftState.getSimTimeMsec()
          % aircraftState.getThisPathIndex()
          % aircraftState.getPosition()[0] % aircraftState.getPosition()[1] % aircraftState.getPosition()[2]
          % euler.x() % euler.y() % euler.z()
          % aircraftState.getRelVel()
          % aircraftState.getPitchCommand()
          % aircraftState.getRollCommand()
          % aircraftState.getThrottleCommand();
      }
    }
  }
};
BOOST_CLASS_VERSION(EvalResults, 1)

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
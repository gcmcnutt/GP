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

#define SIM_MAX_ROLL_RATE_RADSEC (M_PI)
#define SIM_MAX_PITCH_RATE_RADSEC (M_PI)

#define SIM_INITIAL_VELOCITY 10.0
#define SIM_RABBIT_VELOCITY 12.0
#define SIM_THROTTLE_SCALE 5.0
#define SIM_CRASH_PENALTY (100.0 + 100.0 + 100.0)
#define SIM_INITIAL_ALTITUDE -10.0
#define SIM_INITIAL_THROTTLE 0.0
#define SIM_INITIAL_LOCATION_DITHER 30.0
#define SIM_PATH_BOUNDS 40.0
#define SIM_PATH_RADIUS_LIMIT 60.0
#define SIM_MIN_ELEVATION -7.0

#define SIM_TOTAL_TIME_MSEC (30 * 1000)
#define SIM_TIME_STEP_MSEC (200)

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
 * some generic path information about routes
 */
  class Path {
  public:
    Eigen::Vector3d start;
    Eigen::Vector3d orientation;
    double distanceFromStart;
    double radiansFromStart;

    void dump(std::ostream& os) {
      os << format("Path: (%f, %f, %f), Odometer: %f, Turnmeter: %f")
        % start[0] % start[1] % start[2]
        % distanceFromStart
        % radiansFromStart;
    }

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      ar& start;
      ar& orientation;
      ar& distanceFromStart;
      ar& radiansFromStart;
    }
};
BOOST_CLASS_VERSION(Path, 1)

/*
* portable aircraft state
 */
  struct AircraftState {
  public:

    AircraftState() {}
    AircraftState(int thisPathIndex, double relVel, Eigen::Quaterniond orientation,
      Eigen::Vector3d pos, double pc, double rc, double tc,
      unsigned long int timeMsec, bool crashed)
      : thisPathIndex(thisPathIndex), dRelVel(relVel), aircraft_orientation(orientation), position(pos), simTimeMsec(timeMsec), simCrashed(crashed),
      pitchCommand(pc), rollCommand(rc), throttleCommand(tc) {
    }

    // generate setters and getters
    int getThisPathIndex() const { return thisPathIndex; }
    void setThisPathIndex(int index) { thisPathIndex = index; }

    double getRelVel() const { return dRelVel; }
    void setRelVel(double relVel) { dRelVel = relVel; }

    Eigen::Quaterniond getOrientation() const { return aircraft_orientation; }
    void setOrientation(Eigen::Quaterniond orientation) { aircraft_orientation = orientation; }

    Eigen::Vector3d getPosition() const { return position; }
    void setPosition(Eigen::Vector3d pos) { position = pos; }

    unsigned long int getSimTimeMsec() const { return simTimeMsec; }
    void setSimTimeMsec(unsigned long int timeMsec) { simTimeMsec = timeMsec; }

    bool isSimCrashed() const { return simCrashed; }
    void setSimCrashed(bool crashed) { simCrashed = crashed; }

    double getPitchCommand() const { return pitchCommand; }
    double setPitchCommand(double pitch) { return (pitchCommand = std::clamp(pitch, -1.0, 1.0)); }
    double getRollCommand() const { return rollCommand; }
    double setRollCommand(double roll) { return (rollCommand = std::clamp(roll, -1.0, 1.0)); }
    double getThrottleCommand() const { return throttleCommand; }
    double setThrottleCommand(double throttle) { return (throttleCommand = std::clamp(throttle, -1.0, 1.0)); }

    void minisimAdvanceState(double dt) {
      double dtSec = dt / 1000.0;

      // get current roll state, compute left/right force (positive roll is right)
      double delta_roll = remainder(rollCommand * dtSec * SIM_MAX_ROLL_RATE_RADSEC, M_PI);

      // get current pitch state, compute up/down force (positive pitch is up)
      double delta_pitch = remainder(pitchCommand * dtSec * SIM_MAX_PITCH_RATE_RADSEC, M_PI);

      // adjust velocity as a function of throttle (-1:1)
      dRelVel = SIM_INITIAL_VELOCITY + (throttleCommand * SIM_THROTTLE_SCALE);

      // Convert pitch and roll updates to quaternions (in the body frame)
      Eigen::Quaterniond delta_roll_quat(Eigen::AngleAxisd(delta_roll, Eigen::Vector3d::UnitX()));
      Eigen::Quaterniond delta_pitch_quat(Eigen::AngleAxisd(delta_pitch, Eigen::Vector3d::UnitY()));

      // Apply the roll and pitch adjustments to the aircraft's orientation
      aircraft_orientation = aircraft_orientation * delta_roll_quat;
      aircraft_orientation = aircraft_orientation * delta_pitch_quat;

      // Normalize the resulting quaternion
      aircraft_orientation.normalize();

      // Define the initial velocity vector in the body frame
      Eigen::Vector3d velocity_body(dRelVel * dtSec, 0, 0);

      // Rotate the velocity vector using the updated quaternion
      Eigen::Vector3d velocity_world = aircraft_orientation * velocity_body;

      // adjust position
      position += velocity_world;
    }

  private:
    int thisPathIndex;
    double dRelVel;
    Eigen::Quaterniond aircraft_orientation;
    Eigen::Vector3d position;
    unsigned long int simTimeMsec;
    bool simCrashed;
    double pitchCommand;
    double rollCommand;
    double throttleCommand;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      ar& thisPathIndex;
      ar& dRelVel;
      ar& aircraft_orientation;
      ar& position;
      ar& pitchCommand;
      ar& rollCommand;
      ar& throttleCommand;
      ar& simTimeMsec;
      ar& simCrashed;
    }
};
BOOST_CLASS_VERSION(AircraftState, 1)


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
};


struct EvalResults {
  std::vector<CrashReason> crashReasonList;
  std::vector<std::vector<Path>> pathList;
  std::vector<std::vector<AircraftState>> aircraftStateList;

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& crashReasonList;
    ar& pathList;
    ar& aircraftStateList;
  }

  void dump(std::ostream& os) {
    os << format("EvalResults: crash[%d] paths[%d] states[%d]\n Paths:\n")
      % crashReasonList.size() % pathList.size() % aircraftStateList.size();

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
/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <boost/asio.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>

using boost::asio::ip::tcp;

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
#define SIM_MIN_ELEVATION -3.0

#define SIM_TOTAL_TIME_MSEC (50 * 1000)
#define SIM_TIME_STEP_MSEC (200)

namespace boost {
namespace serialization {

// Serialization for Eigen::Vector3d
template<class Archive>
void serialize(Archive & ar, Eigen::Vector3d & v, const unsigned int version)
{
    ar & boost::serialization::make_array(v.data(), 3);
}

// Serialization for Eigen::Quaterniond
template<class Archive>
void save(Archive & ar, const Eigen::Quaterniond & q, const unsigned int version)
{
    ar & boost::serialization::make_array(q.coeffs().data(), 4);
}

template<class Archive>
void load(Archive & ar, Eigen::Quaterniond & q, const unsigned int version)
{
    Eigen::Vector4d coeffs;
    ar & boost::serialization::make_array(coeffs.data(), 4);
    q = Eigen::Quaterniond(coeffs[3], coeffs[0], coeffs[1], coeffs[2]);
}

} // namespace serialization
} // namespace boost

// This macro tells boost to use the save/load functions we just defined for Quaterniond
BOOST_SERIALIZATION_SPLIT_FREE(Eigen::Quaterniond)

/*
 * this is always sent from sim to main
 * and occasionally sent as a rest from main to sim
 */
struct AircraftState {
  double dRelVel;
  Eigen::Quaterniond aircraft_orientation;
  Eigen::Vector3d position;
  double pitchCommand;
  double rollCommand;
  double throttleCommand;
  unsigned long int simTime;
  bool simCrashed;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & dRelVel;
    ar & aircraft_orientation;
    ar & position;
    ar & pitchCommand;
    ar & rollCommand;
    ar & throttleCommand;
    ar & simTime;
    ar & simCrashed;
  }
};

/*
 * here we send our control signals to the aircraft
 */
struct ControlSignal {
  double pitchCommand;
  double rollCommand;
  double throttleCommand;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & pitchCommand;
    ar & rollCommand;
    ar & throttleCommand;
  }
};

/*
 * generic command from sim to main
 */
enum class ControlType {
  CONTROL_SIGNAL,
  AIRCRAFT_STATE
};

struct MainToSim {
  ControlType controlType;
  union {
    ControlSignal controlSignal;
    AircraftState aircraftState;
  };

  // Default constructor
  MainToSim() : controlType(ControlType::CONTROL_SIGNAL), controlSignal() {}

  // Constructor for ControlSignal
  MainToSim(ControlType type, const ControlSignal& signal)
    : controlType(type), controlSignal(signal) {
    assert(type == ControlType::CONTROL_SIGNAL);
  }

  // Constructor for AircraftState
  MainToSim(ControlType type, const AircraftState& state)
    : controlType(type), aircraftState(state) {
    assert(type == ControlType::AIRCRAFT_STATE);
  }

  // Destructor to handle union cleanup if necessary
  ~MainToSim() {
    if (controlType == ControlType::CONTROL_SIGNAL) {
      controlSignal.~ControlSignal();
    } else if (controlType == ControlType::AIRCRAFT_STATE) {
      aircraftState.~AircraftState();
    }
  }

  // Copy constructor
  MainToSim(const MainToSim& other) : controlType(other.controlType) {
    if (controlType == ControlType::CONTROL_SIGNAL) {
      new (&controlSignal) ControlSignal(other.controlSignal);
    } else if (controlType == ControlType::AIRCRAFT_STATE) {
      new (&aircraftState) AircraftState(other.aircraftState);
    }
  }

  // Move constructor
  MainToSim(MainToSim&& other) noexcept : controlType(other.controlType) {
    if (controlType == ControlType::CONTROL_SIGNAL) {
      new (&controlSignal) ControlSignal(std::move(other.controlSignal));
    } else if (controlType == ControlType::AIRCRAFT_STATE) {
      new (&aircraftState) AircraftState(std::move(other.aircraftState));
    }
  }

  // Assignment operator
  MainToSim& operator=(const MainToSim& other) {
    if (this != &other) {
      this->~MainToSim();
      new (this) MainToSim(other);
    }
    return *this;
  }

  // Move assignment operator
  MainToSim& operator=(MainToSim&& other) noexcept {
    if (this != &other) {
      this->~MainToSim();
      new (this) MainToSim(std::move(other));
    }
    return *this;
  }

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & controlType;
    switch (controlType) {
      case ControlType::CONTROL_SIGNAL:
        ar & controlSignal;
        break;
      case ControlType::AIRCRAFT_STATE:
        ar & aircraftState;
        break;
    }
  }
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
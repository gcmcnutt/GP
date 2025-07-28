/* test sim for aircraft */
#ifndef AIRCRAFT_STATE_H
#define AIRCRAFT_STATE_H

#ifdef GP_BUILD
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <boost/format.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/version.hpp>
#define CLAMP_DEF(v, min, max) std::clamp(v, min, max) 
#define ATAN2_DEF(y, x) std::atan2(y, x)
#define ABS_DEF(v) std::abs(v)
#define SQRT_DEF(v) std::sqrt(v)
#define MIN_DEF(a, b) std::min(a, b)
#define MAX_DEF(a, b) std::max(a, b)
#else
#include <ArduinoEigenDense.h>
// #include <ArduinoEigen/ArduinoEigen/Eigen/Geometry>
#define CLAMP_DEF(v, min, max) ((v) < (min) ? (min) : ((v) > (max) ? (max) : (v)))
#define ATAN2_DEF(y, x) atan2(y, x)
#define ABS_DEF(v) ((v) < 0 ? -(v) : (v))
#define SQRT_DEF(v) sqrt(v)
#define MIN_DEF(a, b) ((a) < (b) ? (a) : (b))
#define MAX_DEF(a, b) ((a) > (b) ? (a) : (b))
#endif

#define SIM_MAX_ROLL_RATE_RADSEC (M_PI)
#define SIM_MAX_PITCH_RATE_RADSEC (M_PI)

#define SIM_INITIAL_VELOCITY 20.0
#define SIM_RABBIT_VELOCITY 22.0
#define SIM_THROTTLE_SCALE 10.0
#define SIM_CRASH_PENALTY (100.0 + 100.0 + 100.0)
#define SIM_INITIAL_ALTITUDE -10.0
#define SIM_INITIAL_THROTTLE 0.0
#define SIM_INITIAL_LOCATION_DITHER 30.0
#define SIM_PATH_BOUNDS 40.0
#define SIM_PATH_RADIUS_LIMIT 60.0
#define SIM_MIN_ELEVATION -7.0
#define SIM_MAX_ELEVATION -110.0

#define SIM_TOTAL_TIME_MSEC (30 * 1000)
#define SIM_TIME_STEP_MSEC (200)
#define SIM_MAX_INTERVAL_MSEC (SIM_TIME_STEP_MSEC * 5)

/*
 * some generic path information about routes
 */
class Path {
public:
  Eigen::Vector3d start;
  Eigen::Vector3d orientation;
  double distanceFromStart;
  double radiansFromStart;

#ifdef GP_BUILD
  void dump(std::ostream& os) {
    os << boost::format("Path: (%f, %f, %f), Odometer: %f, Turnmeter: %f")
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
#endif
};
#ifdef GP_BUILD
BOOST_CLASS_VERSION(Path, 1)
#endif

/*
* portable aircraft state
 */
struct AircraftState {
  public:

    AircraftState() {}
    AircraftState(int thisPathIndex, double relVel, Eigen::Vector3d vel, Eigen::Quaterniond orientation,
      Eigen::Vector3d pos, double pc, double rc, double tc,
      unsigned long int timeMsec)
      : thisPathIndex(thisPathIndex), dRelVel(relVel), velocity(vel), aircraft_orientation(orientation), position(pos), simTimeMsec(timeMsec),
      pitchCommand(pc), rollCommand(rc), throttleCommand(tc) {
    }

    // generate setters and getters
    int getThisPathIndex() const { return thisPathIndex; }
    void setThisPathIndex(int index) { thisPathIndex = index; }

    double getRelVel() const { return dRelVel; }
    void setRelVel(double relVel) { dRelVel = relVel; }
    
    Eigen::Vector3d getVelocity() const { return velocity; }
    void setVelocity(const Eigen::Vector3d& vel) { velocity = vel; }

    Eigen::Quaterniond getOrientation() const { return aircraft_orientation; }
    void setOrientation(Eigen::Quaterniond orientation) { aircraft_orientation = orientation; }

    Eigen::Vector3d getPosition() const { return position; }
    void setPosition(Eigen::Vector3d pos) { position = pos; }

    unsigned long int getSimTimeMsec() const { return simTimeMsec; }
    void setSimTimeMsec(unsigned long int timeMsec) { simTimeMsec = timeMsec; }

    double getPitchCommand() const { return pitchCommand; }
    double setPitchCommand(double pitch) { return (pitchCommand = CLAMP_DEF(pitch, -1.0, 1.0)); }
    double getRollCommand() const { return rollCommand; }
    double setRollCommand(double roll) { return (rollCommand = CLAMP_DEF(roll, -1.0, 1.0)); }
    double getThrottleCommand() const { return throttleCommand; }
    double setThrottleCommand(double throttle) { return (throttleCommand = CLAMP_DEF(throttle, -1.0, 1.0)); }

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

      // Store the actual velocity vector (convert from distance per timestep to velocity)
      velocity = velocity_world / dtSec;

      // adjust position
      position += velocity_world;
    }

  private:
    int thisPathIndex;
    double dRelVel;
    Eigen::Vector3d velocity;  // Actual velocity vector (north, east, down)
    Eigen::Quaterniond aircraft_orientation;
    Eigen::Vector3d position;
    unsigned long int simTimeMsec;
    double pitchCommand;
    double rollCommand;
    double throttleCommand;

#ifdef GP_BUILD
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      ar& thisPathIndex;
      ar& dRelVel;
      ar& velocity;
      ar& aircraft_orientation;
      ar& position;
      ar& pitchCommand;
      ar& rollCommand;
      ar& throttleCommand;
      ar& simTimeMsec;
    }
#endif
};
#ifdef GP_BUILD
BOOST_CLASS_VERSION(AircraftState, 1)
#endif

#endif
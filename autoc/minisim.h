/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <Eigen/Dense>
#include <Eigen/Geometry>

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

#define SIM_TOTAL_TIME 50.0
#define SIM_TIME_STEP 0.2

class Aircraft {
  public:
    // TODO should be private
    double dRelVel; // reltive forward velocity on +x airplane axis m/s

    // world frame for now
    Eigen::Quaterniond aircraft_orientation;

    // NED convention for location x+ north, y+ east, z+ down
    Eigen::Vector3d position;

    // not used yet
    double R_X;     // rotationX
    double R_Y;     // rotationY
    double R_Z;     // rotationZ

    Aircraft(double dRelVel, Eigen::Quaterniond aircraft_orientation, Eigen::Vector3d position, double R_X, double R_Y, double R_Z);
    
    double setPitchCommand(double pitchCommand);
    double getPitchCommand();
    double setRollCommand(double rollCommand);
    double getRollCommand();
    double setThrottleCommand(double throttleCommand);
    double getThrottleCommand();
    void advanceState(double dt);
    void toString(char * output);

  private:

    // aircraft command values
    double pitchCommand;  // -1:1
    double rollCommand;   // -1:1
    double throttleCommand; // -1:1
};

#endif
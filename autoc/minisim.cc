/* test sim for aircraft */
#include <iostream>
#include <cmath>
#include <array>

#include "minisim.h"

using namespace std;

Aircraft::Aircraft(double dRelVel, Eigen::Quaterniond aircraft_orientation, Eigen::Vector3d position, double R_X, double R_Y, double R_Z) {
      this->dRelVel = dRelVel;
      this->aircraft_orientation = aircraft_orientation;
      this->position = position;
      this->R_X = R_X;
      this->R_Y = R_Y;
      this->R_Z = R_Z;
      
      this->pitchCommand = 0;
      this->rollCommand = 0;
      this->throttleCommand = 0;
}

double Aircraft::setPitchCommand(double pitchCommand) {
  this->pitchCommand = pitchCommand;
  return pitchCommand;
}

double Aircraft::getPitchCommand() {
  return std::clamp(pitchCommand, -1.0, 1.0);
}

double Aircraft::setRollCommand(double rollCommand) {
  this->rollCommand = rollCommand;
  return rollCommand;
}

double Aircraft::getRollCommand() {
  return std::clamp(rollCommand, -1.0, 1.0);
}

double Aircraft::setThrottleCommand(double throttleCommand) {
  this->throttleCommand = throttleCommand;
  return throttleCommand;
}

double Aircraft::getThrottleCommand() {
  return std::clamp(throttleCommand, -1.0, +1.0);
}

void Aircraft::advanceState(double dt) {
  // get current roll state, compute left/right force (positive roll is right)
  double delta_roll = remainder(getRollCommand() * dt * SIM_MAX_ROLL_RATE_RADSEC, M_PI);

  // get current pitch state, compute up/down force (positive pitch is up)
  double delta_pitch = remainder(getPitchCommand() * dt * SIM_MAX_PITCH_RATE_RADSEC, M_PI);

  // adjust velocity as a function of throttle (-1:1)
  dRelVel = SIM_INITIAL_VELOCITY + (getThrottleCommand() * SIM_THROTTLE_SCALE);

  // Convert pitch and roll updates to quaternions (in the body frame)
  Eigen::Quaterniond delta_roll_quat(Eigen::AngleAxisd(delta_roll, Eigen::Vector3d::UnitX()));
  Eigen::Quaterniond delta_pitch_quat(Eigen::AngleAxisd(delta_pitch, Eigen::Vector3d::UnitY()));

  // Apply the roll and pitch adjustments to the aircraft's orientation
  aircraft_orientation = aircraft_orientation * delta_roll_quat;
  aircraft_orientation = aircraft_orientation * delta_pitch_quat;

  // Normalize the resulting quaternion
  aircraft_orientation.normalize();
 
  // Define the initial velocity vector in the body frame
  Eigen::Vector3d velocity_body(dRelVel * dt, 0, 0);

  // Rotate the velocity vector using the updated quaternion
  Eigen::Vector3d velocity_world = aircraft_orientation * velocity_body;

  // adjust position
  position += velocity_world;
}

void Aircraft::toString(char *output) {
  sprintf(output, "AircraftState: %f %f %f %f %f %f %f %f %f %f %f  Command: %f %f %f\n", dRelVel, 
    aircraft_orientation.w(), aircraft_orientation.x(), aircraft_orientation.y(), aircraft_orientation.z(),
    position[0], position[1], position[2], R_X, R_Y, R_Z,
    pitchCommand, rollCommand, throttleCommand);
}



/* test sim for aircraft */
#include <iostream>
#include <cmath>
#include <array>

#include "minisim.h"

using namespace std;
using boost::asio::ip::tcp;

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
  void toString(char* output);

private:

  // aircraft command values
  double pitchCommand;  // -1:1
  double rollCommand;   // -1:1
  double throttleCommand; // -1:1
};

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

void Aircraft::toString(char* output) {
  sprintf(output, "AircraftState: %f %f %f %f %f %f %f %f %f %f %f  Command: %f %f %f\n", dRelVel,
    aircraft_orientation.w(), aircraft_orientation.x(), aircraft_orientation.y(), aircraft_orientation.z(),
    position[0], position[1], position[2], R_X, R_Y, R_Z,
    pitchCommand, rollCommand, throttleCommand);
}


class SimProcess {
public:
  Aircraft aircraft = Aircraft(0.0, Eigen::Quaterniond(1, 0, 0, 0), Eigen::Vector3d(0, 0, 0), 0.0, 0.0, 0.0);
  unsigned long int simTime = 0;
  bool simCrashed = false;
  SimProcess(boost::asio::io_context& io_context, unsigned short port) : socket_(io_context) {
    tcp::resolver resolver(io_context);
    auto endpoints = resolver.resolve("localhost", std::to_string(port));
    boost::asio::connect(socket_, endpoints);
  }

  void run() {
    while (true) { // TODO have reply have a controlled loop exit
      AircraftState aircraftState{ aircraft.dRelVel, aircraft.aircraft_orientation, aircraft.position,
          aircraft.getPitchCommand(), aircraft.getRollCommand(), aircraft.getThrottleCommand(), simTime, simCrashed };
      // always send our state
      sendRPC(socket_, aircraftState);

      // ok what does main say to do
      MainToSim mainToSim = receiveRPC<MainToSim>(socket_);

      switch (mainToSim.controlType) {
        // here we get a reset from the main controller, just update local state (reset, etc)
      case ControlType::AIRCRAFT_STATE:
        aircraft.dRelVel = mainToSim.aircraftState.dRelVel;
        aircraft.aircraft_orientation = mainToSim.aircraftState.aircraft_orientation;
        aircraft.position = mainToSim.aircraftState.position;
        aircraft.setPitchCommand(mainToSim.aircraftState.pitchCommand);
        aircraft.setRollCommand(mainToSim.aircraftState.rollCommand);
        aircraft.setThrottleCommand(mainToSim.aircraftState.throttleCommand);
        simTime = mainToSim.aircraftState.simTime;
        simCrashed = false;
        break;

        // here we get some control signals, simulate
      case ControlType::CONTROL_SIGNAL:
        // update controls
        aircraft.setPitchCommand(mainToSim.controlSignal.pitchCommand);
        aircraft.setRollCommand(mainToSim.controlSignal.rollCommand);
        aircraft.setThrottleCommand(mainToSim.controlSignal.throttleCommand);

        // bump simulation a step
        aircraft.advanceState(SIM_TIME_STEP_MSEC / 1000.0);
        simTime += SIM_TIME_STEP_MSEC;

        // send our updated state to evaluator
        aircraftState.dRelVel = aircraft.dRelVel;
        aircraftState.aircraft_orientation = aircraft.aircraft_orientation;
        aircraftState.position = aircraft.position;
        aircraftState.pitchCommand = aircraft.getPitchCommand();
        aircraftState.rollCommand = aircraft.getRollCommand();
        aircraftState.throttleCommand = aircraft.getThrottleCommand();
        aircraftState.simTime = simTime;

        // for now simulate a simulator detected crash
        if (aircraft.position[2] > SIM_MIN_ELEVATION) {
          simCrashed = true;
        }
        else {
          simCrashed = false;
        }
        break;
      }
    }
  }

private:
  tcp::socket socket_;
};

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: subprocess <port>" << std::endl;
    return 1;
  }

  unsigned int instance = std::atoi(argv[1]);
  unsigned short port = std::atoi(argv[2]);
  boost::asio::io_context io_context;
  SimProcess sim_process(io_context, port);
  sim_process.run();

  return 0;
}

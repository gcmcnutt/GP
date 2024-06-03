/* test sim for aircraft */
#include <iostream>
#include "minisim.h"

AircraftState::AircraftState() {
      this->dRelVel = 0;
      this->dPhi = 0;
      this->dTheta = 0;
      this->dPsi = 0;
      this->X = 0;
      this->Y = 0;
      this->Z = 0;
      this->R_X = 0;
      this->R_Y = 0;
      this->R_Z = 0;
} 

void AircraftState::setState(double dRelVel, double dPhi, double dTheta, double dPsi, double X, double Y, double Z, double R_X, double R_Y, double R_Z) {
      this->dRelVel = dRelVel;
      this->dPhi = dPhi;
      this->dTheta = dTheta;
      this->dPsi = dPsi;
      this->X = X;
      this->Y = Y;
      this->Z = Z;
      this->R_X = R_X;
      this->R_Y = R_Y;
      this->R_Z = R_Z;
}


Aircraft::Aircraft(AircraftState *state) {
      this->state = state;
      this->pitchCommand = 0;
      this->rollCommand = 0;
      this->throttleCommand = 0;
}

void Aircraft::setState(AircraftState *state) {
      this->state = state;
}

AircraftState *Aircraft::getState() {
      return state;
}

double Aircraft::setPitchCommand(double pitchCommand) {
  this->pitchCommand = pitchCommand;
  return pitchCommand;
}

double Aircraft::getPitchCommand() {
  return pitchCommand;
}

double Aircraft::setRollCommand(double rollCommand) {
  this->rollCommand = rollCommand;
  return rollCommand;
}

double Aircraft::getRollCommand() {
  return rollCommand;
}

double Aircraft::setThrottleCommand(double throttleCommand) {
  this->throttleCommand = throttleCommand;
  return throttleCommand;
}

double Aircraft::getThrottleCommand() {
  return throttleCommand;
}

void Aircraft::advanceState(double dt) {
  // do some physics
}

void Aircraft::toString(char *output) {
  sprintf(output, "AircraftState: %f %f %f %f %f %f %f %f %f %f  Command: %f %f %f\n", state->dRelVel, state->dPhi, state->dTheta,
    state->dPsi, state->X, state->Y, state->Z, state->R_X, state->R_Y, state->R_Z,
    pitchCommand, rollCommand, throttleCommand);
}

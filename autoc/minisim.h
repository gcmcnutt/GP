/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

class AircraftState {
  public:
    void setState(double dRelVel, double dPhi, double dTheta, double dPsi, double X, double Y, double Z, double R_X, double R_Y, double R_Z);
    AircraftState(); 

    double dRelVel;
    double dPhi;
    double dTheta;
    double dPsi;
    double X;
    double Y;
    double Z;
    double R_X;
    double R_Y;
    double R_Z;
};

class Aircraft {
  public:
    Aircraft(AircraftState *state);
    
    void setState(AircraftState *state);
    AircraftState *getState();
    double setPitchCommand(double pitchCommand);
    double getPitchCommand();
    double setRollCommand(double rollCommand);
    double getRollCommand();
    double setThrottleCommand(double throttleCommand);
    double getThrottleCommand();
    void advanceState(double dt);
    void toString(char * output);

  private:
    AircraftState *state; 

    // aircraft command values
    double pitchCommand;
    double rollCommand;
    double throttleCommand;
};

#endif
/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include "pathgen.h"
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>

#define MAX_DELTA_ANGLE_RADSEC M_PI
#define MAX_SERVO_DEFLECTION 1000.0
#define MAX_THROTTLE_DEFLECTION 1000.0

#define SIM_INITIAL_VELOCITY 5.0
#define SIM_INITIAL_ALTITUDE 10.0
#define SIM_INITIAL_THROTTLE 0.5
#define SIM_PATH_BOUNDS 40.0
#define SIM_PATH_RADIUS_LIMIT 60.0

#define SIM_TOTAL_TIME 50.0
#define SIM_CRASH_FITNESS_PENALTY 10000.0

class AircraftState {
  public:
    AircraftState(double dRelVel, double dPhi, double dTheta, double dPsi, double X, double Y, double Z, double R_X, double R_Y, double R_Z);
    AircraftState(); 

    double dRelVel; // reltive velocity
    double dPhi;    // roll
    double dTheta;  // pitch
    double dPsi;    // yaw
    double X;       // positionX
    double Y;       // positionY
    double Z;       // positionZ
    double R_X;     // rotationX
    double R_Y;     // rotationY
    double R_Z;     // rotationZ
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

class Renderer {
  public:
    Renderer(std::vector<Point3D> path, std::vector<Point3D> actual);
    ~Renderer();
  private:
    vtkSmartPointer<vtkNamedColors> colors;
    vtkSmartPointer<vtkPolyDataMapper> mapper1;
    vtkSmartPointer<vtkPolyDataMapper> mapper2;
    vtkSmartPointer<vtkActor> actor1;
    vtkSmartPointer<vtkActor> actor2;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;
};

#endif
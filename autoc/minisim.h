/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <mutex>

#include "pathgen.h"

#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCommand.h>
#include <vtkCamera.h>
#include <vtkPlaneSource.h>
#include <vtkCellData.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLine.h>

#define MAX_ROLL_RATE_RADSEC M_PI
#define MAX_PITCH_RATE_RADSEC M_PI
#define MAX_YAW_RATE_RADSEC M_PI / 2.0

#define SIM_INITIAL_VELOCITY 4.0
#define SIM_INITIAL_ALTITUDE -10.0
#define SIM_INITIAL_HEADING 0.0
#define SIM_INITIAL_THROTTLE 0.5
#define SIM_PATH_BOUNDS 40.0
#define SIM_PATH_RADIUS_LIMIT 60.0
#define SIM_MIN_ELEVATION -3.0

#define SIM_TOTAL_TIME 75.0
#define SIM_TIME_STEP 0.1

#define SIM_FITNESS_EXPONENT 2.0
#define SIM_CRASH_FITNESS_PENALTY pow(1.0, SIM_FITNESS_EXPONENT)
#define SIM_HALF_DISTANCE SIM_PATH_BOUNDS

class AircraftState {
  public:
    AircraftState(double dRelVel, double dPhi, double dTheta, double dPsi, double X, double Y, double Z, double R_X, double R_Y, double R_Z);
    AircraftState(); 

    double dRelVel; // reltive forward velocity on +x airplane axis m/s

    // world frame for now
    double dPhi;    // roll+ right -pi:pi
    double dTheta;  // pitch+ up -pi:pi
    double dPsi;    // yaw+ clockwise -pi:pi -- ALWAYS 0

    // NED convention for location
    double X;       // positionX+ north
    double Y;       // positionY+ east
    double Z;       // positionZ+ down

    //not used
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
    double pitchCommand;  // -1:1
    double rollCommand;   // -1:1
    double throttleCommand; // 0:1
};

class Renderer : public vtkCommand {
  public:
    void update(std::vector<Point3D> path, std::vector<Point3D> actual);
    void start();
    virtual void Execute(vtkObject* caller, unsigned long eventId, void* vtkNotUsed(callData));

  private:
    
    // Shared resources
    std::mutex dataMutex;
    bool newDataAvailable = false;
    vtkSmartPointer<vtkPolyData> path;
    vtkSmartPointer<vtkPolyData> actual;
    vtkSmartPointer<vtkPolyData> segmentGap;
    vtkSmartPointer<vtkActor> actor1;
    vtkSmartPointer<vtkActor> actor2;
    vtkSmartPointer<vtkActor> actor3;
    vtkSmartPointer<vtkActor> planeActor;
    vtkSmartPointer<vtkPlaneSource> planeSource;

    int TimerCount;

    vtkSmartPointer<vtkPolyData> createPointSet(const std::vector<Point3D> points);
    vtkSmartPointer<vtkPolyData> createSegmentSet(const std::vector<Point3D> start, const std::vector<Point3D> end);
    void RenderInBackground(vtkSmartPointer<vtkRenderWindow> renderWindow);
};

#endif
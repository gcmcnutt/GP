/* test sim for aircraft */
#ifndef MINISIM_H
#define MINISIM_H

#include <mutex>
#include <Eigen/Dense>
#include <Eigen/Geometry>

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
#include <vtkAppendPolyData.h>

#define NUM_SEGMENTS_PER_PATH 18
#define FIELD_SIZE 100.0
#define FIELD_GAP 10.0

#define MAX_ROLL_RATE_RADSEC (M_PI / 2.0)
#define MAX_PITCH_RATE_RADSEC (M_PI / 3.0)

#define SIM_INITIAL_VELOCITY 5.0
#define SIM_THROTTLE_SCALE 3.0

#define SIM_INITIAL_ALTITUDE -10.0
#define SIM_INITIAL_HEADING 0.0
#define SIM_INITIAL_THROTTLE 0.5
#define SIM_PATH_BOUNDS 40.0
#define SIM_PATH_RADIUS_LIMIT 60.0
#define SIM_MIN_ELEVATION -3.0

#define SIM_TOTAL_TIME 75.0
#define SIM_TIME_STEP 0.2

class ExtraConfig {
  public:
    int simNumPathsPerGen = 1;
    
    // // Custom implementation of the << operator for the extraCfg type
    // std::ostream& operator << (std::ostream& os) {
    //   os << "simNumPathsPerGen: " + simNumPathsPerGen;
    //   return os;
    // }
};

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

class Renderer : public vtkCommand {
  public:
    Renderer(ExtraConfig &extraCfg) : extraCfg(extraCfg) {};

    void update();
    void start();
    bool isRunning();
    virtual void Execute(vtkObject* caller, unsigned long eventId, void* vtkNotUsed(callData));
    void addPathElementList(std::vector<Path> plan, std::vector<Path> actual);
    
    std::recursive_mutex dataMutex;
      
    // the path(s) a population will attempt
    std::vector<std::vector<Path>> generationPaths;

    // intermediate paths and results
    std::vector<std::vector<Path>> pathList;
    std::vector<std::vector<Path>> actualList;

    ExtraConfig &extraCfg;

  private:
    
    // Shared resources
    bool newDataAvailable = false;
    bool exitFlag = false;

    vtkSmartPointer<vtkAppendPolyData> paths;
    vtkSmartPointer<vtkAppendPolyData> actuals;
    vtkSmartPointer<vtkAppendPolyData> segmentGaps;
    vtkSmartPointer<vtkAppendPolyData> planeData;

    vtkSmartPointer<vtkActor> actor1;
    vtkSmartPointer<vtkActor> actor2;
    vtkSmartPointer<vtkActor> actor3;
    vtkSmartPointer<vtkActor> planeActor;

    Eigen::Vector3d renderingOffset(int i); // locate a coordinate offset for our rendering screen
    vtkSmartPointer<vtkPolyData> createPointSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points);
    vtkSmartPointer<vtkPolyData> createSegmentSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> start, const std::vector<Eigen::Vector3d> end);
    std::vector<Eigen::Vector3d> pathToVector(const std::vector<Path> path);
    void RenderInBackground(vtkSmartPointer<vtkRenderWindow> renderWindow);

};

#endif
#ifndef RENDERER_H
#define RENDERER_H

#include "minisim.h"

#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPolyLine.h>
#include <vtkObjectFactory.h>
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
#include <vtkInteractorStyleUser.h>
#include <vtkLine.h>
#include <vtkAppendPolyData.h>
#include <vtkTubeFilter.h>
#include <vtkRibbonFilter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#define FIELD_SIZE 100.0
#define FIELD_GAP 10.0

// Custom event IDs
enum {
  NextModelEvent = vtkCommand::UserEvent + 1,
  PreviousModelEvent
};

class CustomInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static CustomInteractorStyle* New();
  vtkTypeMacro(CustomInteractorStyle, vtkInteractorStyleTrackballCamera);

protected:
  void OnChar() override {
    vtkRenderWindowInteractor* rwi = this->Interactor;
    std::string key = rwi->GetKeySym();

    if (key == "n" || key == "N") {
      this->InvokeEvent(NextModelEvent, nullptr);
    }
    else if (key == "p" || key == "P") {
      this->InvokeEvent(PreviousModelEvent, nullptr);
    }
    else {
      vtkInteractorStyleTrackballCamera::OnChar();
    }
  }
};

vtkStandardNewMacro(CustomInteractorStyle);

class Renderer {
public:
  void initialize();
  bool isRunning();
  bool updateGenerationDisplay(int genNumber);

  int genNumber = 0;

  vtkSmartPointer<vtkRenderWindow> renderWindow;
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

private:
  vtkSmartPointer<vtkOrientationMarkerWidget> orientationWidget;
  vtkSmartPointer<vtkAppendPolyData> paths;
  vtkSmartPointer<vtkAppendPolyData> actuals;
  vtkSmartPointer<vtkAppendPolyData> segmentGaps;
  vtkSmartPointer<vtkAppendPolyData> planeData;
  vtkSmartPointer<vtkAppendPolyData> blackboxTapes;

  vtkSmartPointer<vtkActor> actor1;
  vtkSmartPointer<vtkActor> actor2;
  vtkSmartPointer<vtkActor> actor3;
  vtkSmartPointer<vtkActor> blackboxActor;

  Eigen::Vector3d renderingOffset(int i); // locate a coordinate offset for our rendering screen
  vtkSmartPointer<vtkPolyData> createPointSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points);
  vtkSmartPointer<vtkPolyData> createSegmentSet(Eigen::Vector3d offset, const std::vector<AircraftState> state, const std::vector<Eigen::Vector3d> end);
  vtkSmartPointer<vtkPolyData> createTapeSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points, const std::vector<Eigen::Vector3d> normals);
  std::vector<Eigen::Vector3d> pathToVector(const std::vector<Path> path);
  std::vector<Eigen::Vector3d> stateToVector(const std::vector<AircraftState> path);
  std::vector<Eigen::Vector3d> stateToOrientation(const std::vector<AircraftState> state);
};

#endif

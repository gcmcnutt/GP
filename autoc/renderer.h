#ifndef RENDERER_H
#define RENDERER_H

#include <mutex>
#include <thread>

#include "pathgen.h"
#include "autoc.h"

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
#include <vtkTubeFilter.h>
#include <vtkRibbonFilter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

class Renderer : public vtkCommand {
public:
  Renderer(std::string computedKeyName) : computedKeyName(computedKeyName) {};

  void update();
  void start();
  bool isRunning();
  virtual void Execute(vtkObject* caller, unsigned long eventId, void* vtkNotUsed(callData));
  void addPathElementList(std::vector<Path> plan, std::vector<Path> actual);

  std::recursive_mutex dataMutex;

  // intermediate paths and results
  std::vector<std::vector<Path>> pathList;
  std::vector<std::vector<Path>> actualList;

private:
  std::string computedKeyName;

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
  vtkSmartPointer<vtkPolyData> createTapeSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points,
    const std::vector<Eigen::Vector3d> normals);
  std::vector<Eigen::Vector3d> pathToVector(const std::vector<Path> path);
  std::vector<Eigen::Vector3d> pathToOrientation(const std::vector<Path> path);
  void RenderInBackground(vtkSmartPointer<vtkRenderWindow> renderWindow);
};

#endif

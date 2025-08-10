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
#include <vtkPolyDataMapper2D.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCommand.h>
#include <vtkCallbackCommand.h>
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
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkRendererCollection.h>
#include <chrono>

#define FIELD_SIZE 100.0
#define FIELD_GAP 10.0
#define MSPRCOVERRIDE_FLAG "MSPRCOVERRIDE"

// Forward declarations
class Renderer;

// Custom event IDs
enum {
  NextModelEvent = vtkCommand::UserEvent + 1,
  PreviousModelEvent,
  NewestModelEvent,
  OldestModelEvent,
  NextTestEvent,
  PreviousTestEvent,
  AllFlightEvent,
  PlaybackEvent
};

// CustomInteractorStyle forward declaration
class CustomInteractorStyle;

// Test span structure for MSPRCOVERRIDE segments
struct TestSpan {
  size_t startIndex;
  size_t endIndex;
  unsigned long startTime;
  unsigned long endTime;
};

class Renderer {
public:
  void initialize();
  bool isRunning();
  bool updateGenerationDisplay(int genNumber);
  void updateTextDisplay(int generation, double fitness);
  void jumpToNewestGeneration();
  void jumpToOldestGeneration();
  void nextTest();
  void previousTest();
  void showAllFlight();
  void extractTestSpans();
  void togglePlaybackAnimation();
  void updatePlaybackAnimation();
  void pausePlaybackAnimation();
  void resumePlaybackAnimation();
  void renderFullScene(); // Render complete scene without S3 fetch
  void hideStopwatch();

  int genNumber = 0;
  
  // Store current generation and fitness for resize updates
  int currentGeneration = 0;
  double currentFitness = 0.0;
  
  // Test span navigation state
  std::vector<TestSpan> testSpans;
  int currentTestIndex = 0;
  bool showingFullFlight = false;
  bool inDecodeMode = false;
  
  // Animation state
  bool isPlaybackActive = false;
  bool isPlaybackPaused = false;
  std::chrono::steady_clock::time_point animationStartTime;
  std::chrono::steady_clock::time_point pauseStartTime;
  std::chrono::duration<double> totalPausedTime = std::chrono::duration<double>::zero();
  double animationSpeed = 1.0; // seconds per animation second
  double totalAnimationDuration = 10.0; // total animation duration in seconds
  unsigned long animationTimerId = 0; // VTK timer ID for animation

  vtkSmartPointer<vtkRenderWindow> renderWindow;
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

private:
  vtkSmartPointer<vtkOrientationMarkerWidget> orientationWidget;
  vtkSmartPointer<vtkAppendPolyData> paths;
  vtkSmartPointer<vtkAppendPolyData> actuals;
  vtkSmartPointer<vtkAppendPolyData> segmentGaps;
  vtkSmartPointer<vtkAppendPolyData> planeData;
  vtkSmartPointer<vtkAppendPolyData> blackboxTapes;
  vtkSmartPointer<vtkAppendPolyData> blackboxHighlightTapes;  // For highlighted test spans

  vtkSmartPointer<vtkActor> actor1;
  vtkSmartPointer<vtkActor> actor2;
  vtkSmartPointer<vtkActor> actor3;
  vtkSmartPointer<vtkActor> blackboxActor;
  vtkSmartPointer<vtkActor> blackboxHighlightActor;  // For highlighted test spans
  
  vtkSmartPointer<vtkTextActor> generationTextActor;
  vtkSmartPointer<vtkTextActor> generationValueActor;
  vtkSmartPointer<vtkTextActor> fitnessTextActor;
  vtkSmartPointer<vtkTextActor> fitnessValueActor;
  vtkSmartPointer<vtkTextActor> testTextActor;
  vtkSmartPointer<vtkTextActor> testValueActor;
  
  // Stopwatch components
  vtkSmartPointer<vtkActor2D> stopwatchActor;
  vtkSmartPointer<vtkTextActor> stopwatchTimeActor;
  bool stopwatchVisible = false;
  double stopwatchTime = 0.0;

  Eigen::Vector3d renderingOffset(int i); // locate a coordinate offset for our rendering screen
  vtkSmartPointer<vtkPolyData> createPointSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points);
  vtkSmartPointer<vtkPolyData> createPointSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points, double timeProgress);
  vtkSmartPointer<vtkPolyData> createSegmentSet(Eigen::Vector3d offset, const std::vector<AircraftState> state, const std::vector<Eigen::Vector3d> end);
  vtkSmartPointer<vtkPolyData> createSegmentSet(Eigen::Vector3d offset, const std::vector<AircraftState> state, const std::vector<Eigen::Vector3d> end, double timeProgress);
  vtkSmartPointer<vtkPolyData> createTapeSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points, const std::vector<Eigen::Vector3d> normals);
  vtkSmartPointer<vtkPolyData> createTapeSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points, const std::vector<Eigen::Vector3d> normals, double timeProgress);
  std::vector<Eigen::Vector3d> pathToVector(const std::vector<Path> path);
  std::vector<Eigen::Vector3d> stateToVector(const std::vector<AircraftState> path);
  std::vector<Eigen::Vector3d> stateToOrientation(const std::vector<AircraftState> state);
  double extractFitnessFromGP(const std::vector<char>& gpData);
  void createHighlightedFlightTapes(Eigen::Vector3d offset);
  void createStopwatch();
  void updateStopwatch(double currentTime);
};

// CustomInteractorStyle implementation
class CustomInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
  static CustomInteractorStyle* New();
  vtkTypeMacro(CustomInteractorStyle, vtkInteractorStyleTrackballCamera);
  
  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

protected:
  void OnChar() override {
    vtkRenderWindowInteractor* rwi = this->Interactor;
    std::string key = rwi->GetKeySym();

    if (key == "n") {
      this->InvokeEvent(NextModelEvent, nullptr);
    }
    else if (key == "N") {
      this->InvokeEvent(NewestModelEvent, nullptr);
    }
    else if (key == "p") {
      this->InvokeEvent(PreviousModelEvent, nullptr);
    }
    else if (key == "P") {
      this->InvokeEvent(OldestModelEvent, nullptr);
    }
    else if (key == "t") {
      this->InvokeEvent(NextTestEvent, nullptr);
    }
    else if (key == "r") {
      this->InvokeEvent(PreviousTestEvent, nullptr);
    }
    else if (key == "a") {
      this->InvokeEvent(AllFlightEvent, nullptr);
    }
    else if (key == "space") {
      this->InvokeEvent(PlaybackEvent, nullptr);
    }
    else {
      vtkInteractorStyleTrackballCamera::OnChar();
    }
  }
  
  void OnLeftButtonDown() override {
    if (renderer_ && renderer_->isPlaybackActive) {
      renderer_->pausePlaybackAnimation();
    }
    vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
  }
  
  void OnLeftButtonUp() override {
    vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
    if (renderer_ && renderer_->isPlaybackActive) {
      renderer_->resumePlaybackAnimation();
    }
  }
  
  void OnMiddleButtonDown() override {
    if (renderer_ && renderer_->isPlaybackActive) {
      renderer_->pausePlaybackAnimation();
    }
    vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
  }
  
  void OnMiddleButtonUp() override {
    vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
    if (renderer_ && renderer_->isPlaybackActive) {
      renderer_->resumePlaybackAnimation();
    }
  }
  
  void OnRightButtonDown() override {
    if (renderer_ && renderer_->isPlaybackActive) {
      renderer_->pausePlaybackAnimation();
    }
    vtkInteractorStyleTrackballCamera::OnRightButtonDown();
  }
  
  void OnRightButtonUp() override {
    vtkInteractorStyleTrackballCamera::OnRightButtonUp();
    if (renderer_ && renderer_->isPlaybackActive) {
      renderer_->resumePlaybackAnimation();
    }
  }

private:
  Renderer* renderer_ = nullptr;
};

vtkStandardNewMacro(CustomInteractorStyle);

#endif

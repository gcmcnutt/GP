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
#include <vtkLineSource.h>
#include <vtkArrowSource.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkTubeFilter.h>
#include <vtkRibbonFilter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkRendererCollection.h>
#include <array>
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

// Forward declarations for timestamped structures (defined in renderer.cc)
struct TimestampedVec;

// Test span structure for MSPRCOVERRIDE segments
struct TestSpan {
  size_t startIndex;
  size_t endIndex;
  unsigned long startTime;
  unsigned long endTime;
  gp_vec3 origin;  // Test origin for xiao mode
  std::vector<TimestampedVec> vecPoints;  // Vec arrows for this span

  TestSpan() : startIndex(0), endIndex(0), startTime(0), endTime(0), origin(0.0f, 0.0f, 0.0f) {}
  TestSpan(size_t start, size_t end, unsigned long stime, unsigned long etime)
    : startIndex(start), endIndex(end), startTime(stime), endTime(etime), origin(0.0f, 0.0f, 0.0f) {}
  TestSpan(size_t start, size_t end, unsigned long stime, unsigned long etime, gp_vec3 orig)
    : startIndex(start), endIndex(end), startTime(stime), endTime(etime), origin(orig) {}
};

class Renderer {
public:
  void initialize();
  bool isRunning();
  bool updateGenerationDisplay(int genNumber);
  void updateTextDisplay(int generation, gp_scalar fitness);
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
  void updateStopwatchPosition();
  void toggleFocusMode();
  void adjustFocusArena(int delta);
  void focusMoveLeft();
  void focusMoveRight();
  void focusMoveUp();
  void focusMoveDown();

  int genNumber = 0;
  
  // Store current generation and fitness for resize updates
  int currentGeneration = 0;
  gp_scalar currentFitness = 0.0f;
  
  // Test span navigation state
  std::vector<TestSpan> testSpans;
  int currentTestIndex = 0;
  bool showingFullFlight = false;
  bool inDecodeMode = false;
  bool inXiaoMode = false;
  bool focusMode = false;
  int focusArenaIndex = 0;
  std::array<gp_scalar,3> focusCameraPosition{0.0f,0.0f,0.0f};
  std::array<gp_scalar,3> focusCameraFocalPoint{0.0f,0.0f,0.0f};
  std::array<gp_scalar,3> focusCameraViewUp{0.0f,0.0f,-1.0f};
  
  // Animation state
  bool isPlaybackActive = false;
  bool isPlaybackPaused = false;
  std::chrono::steady_clock::time_point animationStartTime;
  std::chrono::steady_clock::time_point pauseStartTime;
  std::chrono::duration<gp_scalar> totalPausedTime = std::chrono::duration<gp_scalar>::zero();
  gp_scalar animationSpeed = 1.0f; // seconds per animation second
  gp_scalar totalAnimationDuration = 10.0f; // total animation duration in seconds
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
  vtkSmartPointer<vtkAppendPolyData> xiaoVecArrows;  // For xiao vec vectors

  vtkSmartPointer<vtkActor> actor1;
  vtkSmartPointer<vtkActor> actor2;
  vtkSmartPointer<vtkActor> actor3;
  vtkSmartPointer<vtkActor> blackboxActor;
  vtkSmartPointer<vtkActor> blackboxHighlightActor;  // For highlighted test spans
  vtkSmartPointer<vtkActor> xiaoVecActor;  // For xiao vec arrows
  std::vector<vtkSmartPointer<vtkActor>> arenaLabelActors;
  
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
  gp_scalar stopwatchTime = 0.0f;

  // Control HUD components (stick + throttle)
  vtkSmartPointer<vtkActor2D> controlOutlineActor;
  vtkSmartPointer<vtkActor2D> controlStickActor;
  vtkSmartPointer<vtkActor2D> throttleFillActor;
  vtkSmartPointer<vtkTextActor> controlSourceActor;
  vtkSmartPointer<vtkActor2D> attitudeSkyActor;
  vtkSmartPointer<vtkActor2D> attitudeGroundActor;
  vtkSmartPointer<vtkActor2D> attitudeOutlineActor;
  vtkSmartPointer<vtkTextActor> velocityActor;
  bool controlsVisible = false;
  gp_scalar lastControlPitch = 0.0f;
  gp_scalar lastControlRoll = 0.0f;
  gp_scalar lastControlThrottle = 0.0f;
  gp_scalar lastControlsTime = 0.0f;

  gp_vec3 renderingOffset(int i); // locate a coordinate offset for our rendering screen
  vtkSmartPointer<vtkPolyData> createPointSet(gp_vec3 offset, const std::vector<gp_vec3> points);
  vtkSmartPointer<vtkPolyData> createPointSet(gp_vec3 offset, const std::vector<gp_vec3> points, gp_scalar timeProgress);
  vtkSmartPointer<vtkPolyData> createSegmentSet(gp_vec3 offset, const std::vector<AircraftState> state, const std::vector<gp_vec3> end);
  vtkSmartPointer<vtkPolyData> createSegmentSet(gp_vec3 offset, const std::vector<AircraftState> state, const std::vector<gp_vec3> end, gp_scalar timeProgress);
  vtkSmartPointer<vtkPolyData> createTapeSet(gp_vec3 offset, const std::vector<gp_vec3> points, const std::vector<gp_vec3> normals);
  vtkSmartPointer<vtkPolyData> createTapeSet(gp_vec3 offset, const std::vector<gp_vec3> points, const std::vector<gp_vec3> normals, gp_scalar timeProgress);
  std::vector<gp_vec3> pathToVector(const std::vector<Path> path);
  std::vector<gp_vec3> stateToVector(const std::vector<AircraftState> path);
  std::vector<gp_vec3> stateToOrientation(const std::vector<AircraftState> state);
  gp_scalar extractFitnessFromGP(const std::vector<char>& gpData);
  void createHighlightedFlightTapes(gp_vec3 offset);
  void createStopwatch();
  void updateStopwatch(gp_scalar currentTime);
  void createControlsOverlay();
  void updateControlsOverlay(gp_scalar currentTime);
  void updateControlsPosition();
  bool getControlStateAtTime(gp_scalar currentTime, gp_scalar& pitch, gp_scalar& roll, gp_scalar& throttle, int arenaIndex, bool& usedBlackbox, const AircraftState*& chosenState);
  void setFocusArena(int arenaIdx);
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
    else if (key == "f") {
      if (renderer_) {
        renderer_->toggleFocusMode();
      }
    }
    else if (key == "Left") {
      if (renderer_) {
        renderer_->focusMoveLeft();
      }
    }
    else if (key == "Right") {
      if (renderer_) {
        renderer_->focusMoveRight();
      }
    }
    else if (key == "Up") {
      if (renderer_) {
        renderer_->focusMoveUp();
      }
    }
    else if (key == "Down") {
      if (renderer_) {
        renderer_->focusMoveDown();
      }
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

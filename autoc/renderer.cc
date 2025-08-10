#include <regex>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iomanip>

#include "renderer.h"
#include "config_manager.h"
#include "autoc.h"

#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <gp.h>

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/ListObjectsV2Request.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/client/ClientConfiguration.h>

EvalResults evalResults;
std::string computedKeyName = "";
Renderer renderer;

// Blackbox-related global variables
std::string decoderCommand = "";
vtkSmartPointer<vtkAppendPolyData> blackboxTape;
std::vector<Eigen::Vector3d> blackboxPoints;
std::vector<Eigen::Vector3d> blackboxNormals;
std::vector<AircraftState> blackboxAircraftStates;
std::vector<AircraftState> fullBlackboxAircraftStates;  // Store full flight for span extraction
Eigen::Vector3d blackboxOrigin(0.0, 0.0, 0.0);
double blackboxTimeOffset = 0.0;
std::vector<std::string> csvLines;  // Store CSV lines for span analysis

// Forward declarations
bool parseBlackboxData(const std::string& csvData);
bool loadBlackboxData();
void updateBlackboxForCurrentTest();
void printUsage(const char* progName);

std::shared_ptr<Aws::S3::S3Client> getS3Client() {
  return ConfigManager::getS3Client();
}

void PrintPolyDataInfo(vtkPolyData* polyData)
{
  vtkPoints* points = polyData->GetPoints();
  vtkCellArray* cells = polyData->GetPolys();
  vtkIdList* idList = vtkIdList::New();

  if (points == NULL) {
    printf("No points in the polydata\n");
  }
  else {
    printf("Number of Points: %lld\n", static_cast<long long>(points->GetNumberOfPoints()));
    printf("Points:\n");
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i)
    {
      double p[3];
      points->GetPoint(i, p);
      printf("  Point %lld: (%f, %f, %f)\n", static_cast<long long>(i), p[0], p[1], p[2]);
    }
  }

  printf("Number of Cells: %lld\n", static_cast<long long>(cells->GetNumberOfCells()));
  printf("Cells:\n");
  cells->InitTraversal();
  while (cells->GetNextCell(idList))
  {
    printf("  Cell with %lld points: ", static_cast<long long>(idList->GetNumberOfIds()));
    for (vtkIdType j = 0; j < idList->GetNumberOfIds(); ++j)
    {
      printf("%lld ", static_cast<long long>(idList->GetId(j)));
    }
    printf("\n");
  }

  idList->Delete();
}

// Function to lay out the squares
// given i, and NUM_PATHS_PER_GEN, compute the offsets for this particular square
Eigen::Vector3d Renderer::renderingOffset(int i) {
  // Calculate the dimension of the larger square
  int sideLength = std::ceil(std::sqrt(evalResults.pathList.size()));

  int row = i / sideLength;
  int col = i % sideLength;

  // for now put them in a line
  double xOffset = col * (FIELD_SIZE + FIELD_GAP);
  double yOffset = row * (FIELD_SIZE + FIELD_GAP);

  return Eigen::Vector3d(xOffset, yOffset, 0.0);
}


vtkSmartPointer<vtkPolyData> Renderer::createPointSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points) {
  return createPointSet(offset, points, 1.0);
}

vtkSmartPointer<vtkPolyData> Renderer::createPointSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points, double timeProgress) {
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  
  // Calculate how many points to include based on time progress
  size_t numPointsToShow = static_cast<size_t>(points.size() * timeProgress);
  if (numPointsToShow == 0 && timeProgress > 0.0) numPointsToShow = 1; // Show at least one point if progress > 0
  if (numPointsToShow > points.size()) numPointsToShow = points.size();
  
  // Return empty polydata if no points to show
  if (numPointsToShow == 0 || points.empty()) {
    vtkSmartPointer<vtkPoints> emptyPoints = vtkSmartPointer<vtkPoints>::New();
    polyData->SetPoints(emptyPoints);
    return polyData;
  }
  
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  for (size_t i = 0; i < numPointsToShow; ++i) {
    Eigen::Vector3d rPoint = points[i] + offset;
    vtp->InsertNextPoint(rPoint[0], rPoint[1], rPoint[2]);
  }

  vtkSmartPointer<vtkPolyLine> lines = vtkSmartPointer<vtkPolyLine>::New();
  lines->GetPointIds()->SetNumberOfIds(numPointsToShow);
  for (size_t i = 0; i < numPointsToShow; ++i) {
    lines->GetPointIds()->SetId(i, i);
  }

  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(lines);

  polyData->SetPoints(vtp);
  polyData->SetLines(cells);

  return polyData;
}

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(Eigen::Vector3d offset, const std::vector<AircraftState> state, const std::vector<Eigen::Vector3d> end) {
  return createSegmentSet(offset, state, end, 1.0);
}

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(Eigen::Vector3d offset, const std::vector<AircraftState> state, const std::vector<Eigen::Vector3d> end, double timeProgress) {
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  
  // Calculate how many segments to show based on time progress
  size_t numStatesToShow = static_cast<size_t>(state.size() * timeProgress);
  if (numStatesToShow == 0 && timeProgress > 0.0) numStatesToShow = 1;
  if (numStatesToShow > state.size()) numStatesToShow = state.size();
  
  // Return empty polydata if no states to show
  if (numStatesToShow == 0 || state.empty()) {
    vtkSmartPointer<vtkPoints> emptyPoints = vtkSmartPointer<vtkPoints>::New();
    polyData->SetPoints(emptyPoints);
    return polyData;
  }
  
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (size_t i = 0; i < numStatesToShow; i++) {
    auto& s = state.at(i);
    Eigen::Vector3d rStart = Eigen::Vector3d{ s.getPosition()[0], s.getPosition()[1], s.getPosition()[2] } + offset;
    Eigen::Vector3d rEnd = end[s.getThisPathIndex()] + offset;
    points->InsertNextPoint(rStart[0], rStart[1], rStart[2]);
    points->InsertNextPoint(rEnd[0], rEnd[1], rEnd[2]);
  }

  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  for (int i = 0; i < points->GetNumberOfPoints(); i += 2) {
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, i);
    line->GetPointIds()->SetId(1, i + 1);
    lines->InsertNextCell(line);
  }

  polyData->SetPoints(points);
  polyData->SetLines(lines);

  return polyData;
}

/*
 ** actual data is rendered as a tape with a top and bottom
 */
vtkSmartPointer<vtkPolyData> Renderer::createTapeSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points,
  const std::vector<Eigen::Vector3d> normals) {
  return createTapeSet(offset, points, normals, 1.0);
}

vtkSmartPointer<vtkPolyData> Renderer::createTapeSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points,
  const std::vector<Eigen::Vector3d> normals, double timeProgress) {
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkDoubleArray> orientations = vtkSmartPointer<vtkDoubleArray>::New();
  orientations->SetNumberOfComponents(3);
  orientations->SetName("Orientations");

  // Calculate how many points to show based on time progress
  size_t numPointsToShow = static_cast<size_t>(points.size() * timeProgress);
  if (numPointsToShow == 0 && timeProgress > 0.0) numPointsToShow = 1;
  if (numPointsToShow > points.size()) numPointsToShow = points.size();
  
  // Need at least 2 points for ribbon filter to work properly
  if (numPointsToShow < 2 && points.size() >= 2) {
    if (timeProgress > 0.5) numPointsToShow = 2;
    else if (numPointsToShow < 1) numPointsToShow = 0;
  }

  if (numPointsToShow == 0) {
    // Return empty polydata
    vtkSmartPointer<vtkPolyData> emptyPolyData = vtkSmartPointer<vtkPolyData>::New();
    return emptyPolyData;
  }

  vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
  polyLine->GetPointIds()->SetNumberOfIds(numPointsToShow);

  for (size_t i = 0; i < numPointsToShow; ++i) {
    Eigen::Vector3d point = points[i] + offset;
    vtp->InsertNextPoint(point[0], point[1], point[2]);
    polyLine->GetPointIds()->SetId(i, i);

    // Use the path's orientation as the normal for the ribbon
    orientations->InsertNextTuple(normals[i].data());
  }

  lines->InsertNextCell(polyLine);

  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(vtp);
  polyData->SetLines(lines);
  polyData->GetPointData()->SetNormals(orientations);

  // Create a ribbon filter only if we have enough points
  if (numPointsToShow >= 2) {
    vtkSmartPointer<vtkRibbonFilter> ribbonFilter = vtkSmartPointer<vtkRibbonFilter>::New();
    ribbonFilter->SetInputData(polyData);
    ribbonFilter->SetWidth(0.5);
    ribbonFilter->Update();
    return ribbonFilter->GetOutput();
  } else {
    return polyData;
  }
}

/*
 * re-render the generation number asked for
 */
bool Renderer::updateGenerationDisplay(int newGen) {
  // now do initial fetch
  Aws::S3::Model::GetObjectRequest request;
  request.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
  std::string keyName = computedKeyName + "gen" + std::to_string(newGen) + ".dmp";
  request.SetKey(keyName);
  auto outcome = getS3Client()->GetObject(request);
  std::cerr << "Fetched " << keyName << " result " << outcome.IsSuccess() << std::endl;
  if (outcome.IsSuccess()) {
    std::ostringstream oss;
    oss << outcome.GetResult().GetBody().rdbuf();
    std::string retrievedData = oss.str();

    // Deserialize the data
    try {
      std::istringstream iss(retrievedData);
      boost::archive::text_iarchive ia(iss);
      ia >> evalResults;
    }
    catch (const std::exception& e) {
      std::cerr << "Error during deserialization: " << e.what() << std::endl;
      return false;
    }
  }
  else {
    std::cerr << "Error retrieving object " << keyName << " from S3: " << outcome.GetError().GetMessage() << std::endl;
    return false;
  }

  // Clear the existing data
  this->paths->RemoveAllInputs();
  this->actuals->RemoveAllInputs();
  this->segmentGaps->RemoveAllInputs();
  this->planeData->RemoveAllInputs();
  this->blackboxTapes->RemoveAllInputs();
  this->blackboxHighlightTapes->RemoveAllInputs();
  
  // Always ensure highlight tapes has at least empty data to prevent VTK pipeline errors
  vtkNew<vtkPolyData> emptyHighlightData;
  vtkNew<vtkPoints> emptyHighlightPoints;
  emptyHighlightData->SetPoints(emptyHighlightPoints);
  this->blackboxHighlightTapes->AddInputData(emptyHighlightData);

  int maxArenas = (!blackboxAircraftStates.empty()) ? 1 : evalResults.pathList.size();
  for (int i = 0; i < maxArenas; i++) {
    Eigen::Vector3d offset = renderingOffset(i);

    std::vector<Eigen::Vector3d> p = pathToVector(evalResults.pathList[i]);
    std::vector<Eigen::Vector3d> a = stateToVector(evalResults.aircraftStateList[i]);

    if (!p.empty()) {
      this->paths->AddInputData(createPointSet(offset, p));
    }
    if (!a.empty()) {
      this->actuals->AddInputData(createTapeSet(offset, a, stateToOrientation(evalResults.aircraftStateList[i])));
    }
    if (!a.empty() && !p.empty()) {
      this->segmentGaps->AddInputData(createSegmentSet(offset, evalResults.aircraftStateList[i], p));
    }

    // Create a plane source at z = 0
    vtkNew<vtkPlaneSource> planeSource;

    double width = FIELD_SIZE;
    double height = FIELD_SIZE;
    int resolution = FIELD_SIZE / 10.0;
    planeSource->SetOrigin(offset[0] - width / 2.0, offset[1] - height / 2.0, 0.0);
    planeSource->SetPoint1(offset[0] + width / 2.0, offset[1] - height / 2.0, 0.0);
    planeSource->SetPoint2(offset[0] - width / 2.0, offset[1] + height / 2.0, 0.0);
    planeSource->SetXResolution(resolution);
    planeSource->SetYResolution(resolution);
    planeSource->Update();

    // Create cell data.
    vtkNew<vtkUnsignedCharArray> cellData;
    cellData->SetNumberOfComponents(4);
    cellData->SetNumberOfTuples(planeSource->GetOutput()->GetNumberOfCells());

    // checkerboard
    for (int i = 0; i < planeSource->GetOutput()->GetNumberOfCells(); i++) {
      if (i % 2 ^ (i / 10) % 2) {
        double rgb[4] = { 255.0, 255.0, 255.0, 100.0 };
        cellData->InsertTuple(i, rgb);
      }
      else {
        double rgb[4] = { 0.0, 0.0, 0.0, 100.0 };
        cellData->InsertTuple(i, rgb);
      }
    }
    planeSource->GetOutput()->GetCellData()->SetScalars(cellData);
    planeSource->Update();
    planeData->AddInputConnection(planeSource->GetOutputPort());
    
    // Add blackbox data to first arena only
    if (i == 0 && !blackboxAircraftStates.empty()) {
      // Center blackbox data in the arena by adding the arena offset
      // The blackbox data is already origin-centered, so we add the arena offset to position it properly
      Eigen::Vector3d blackboxOffset = offset;
      
      // Use the same rendering pipeline as the blue/yellow tape
      std::vector<Eigen::Vector3d> a = stateToVector(blackboxAircraftStates);
      
      // Only create tape if we have enough points
      if (a.size() >= 2) {
        if (inDecodeMode && !testSpans.empty() && !showingFullFlight) {
          // Create regular tape for current test span only
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a, stateToOrientation(blackboxAircraftStates)));
          // Reset opacity to full brightness for single span display
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        } else if (inDecodeMode && showingFullFlight && !testSpans.empty()) {
          // Show full flight with highlighted test spans
          createHighlightedFlightTapes(blackboxOffset);
        } else {
          // Regular blackbox rendering (no test spans or not in decode mode)
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a, stateToOrientation(blackboxAircraftStates)));
          // Reset opacity to full brightness for regular display
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        }
      }
    }
  }

  this->planeData->Update();
  this->paths->Update();
  this->actuals->Update();
  this->segmentGaps->Update();
  
  // Only update blackbox tapes if there's blackbox data
  if (!blackboxAircraftStates.empty()) {
    this->blackboxTapes->Update();
    this->blackboxHighlightTapes->Update();
  }

  // Update the window title
  std::string title = keyName + " - " + std::to_string(10000 - newGen);
  renderWindow->SetWindowName(title.c_str());

  // Extract fitness from GP data and update text displays
  double fitness = extractFitnessFromGP(evalResults.gp);
  updateTextDisplay(newGen, fitness);

  // Render the updated scene
  renderWindow->Render();
  return true;
}

// Custom command classes
class NextModelCommand : public vtkCommand {
public:
  static NextModelCommand* New() {
    return new NextModelCommand();
  }
  vtkTypeMacro(NextModelCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Load next results..." << std::endl;
      if (renderer_->updateGenerationDisplay(renderer_->genNumber - 1)) {
        renderer_->genNumber--;
      }
    }
  }

protected:
  NextModelCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class PreviousModelCommand : public vtkCommand {
public:
  static PreviousModelCommand* New() {
    return new PreviousModelCommand();
  }
  vtkTypeMacro(PreviousModelCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Load previous results..." << std::endl;
      if (renderer_->updateGenerationDisplay(renderer_->genNumber + 1)) {
        renderer_->genNumber++;
      }
    }
  }

protected:
  PreviousModelCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class ExitCommand : public vtkCommand
{
public:
  static ExitCommand* New() { return new ExitCommand; }
  virtual void Execute(vtkObject*, unsigned long, void*)
  {
    std::cout << "Exit requested. Stopping interactor." << std::endl;
    if (this->Interactor)
      this->Interactor->TerminateApp();
  }
  vtkRenderWindowInteractor* Interactor;
};

class WindowResizeCommand : public vtkCommand {
public:
  static WindowResizeCommand* New() {
    return new WindowResizeCommand();
  }
  vtkTypeMacro(WindowResizeCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      // Update text display with current values when window is resized
      renderer_->updateTextDisplay(renderer_->currentGeneration, renderer_->currentFitness);
    }
  }

protected:
  WindowResizeCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class NewestModelCommand : public vtkCommand {
public:
  static NewestModelCommand* New() {
    return new NewestModelCommand();
  }
  vtkTypeMacro(NewestModelCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Jump to newest generation..." << std::endl;
      // Find the newest generation by getting the lowest generation number (highest actual generation)
      renderer_->jumpToNewestGeneration();
    }
  }

protected:
  NewestModelCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class OldestModelCommand : public vtkCommand {
public:
  static OldestModelCommand* New() {
    return new OldestModelCommand();
  }
  vtkTypeMacro(OldestModelCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Jump to oldest generation..." << std::endl;
      // Jump to generation 1 (stored as gen9999.dmp)
      renderer_->jumpToOldestGeneration();
    }
  }

protected:
  OldestModelCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class NextTestCommand : public vtkCommand {
public:
  static NextTestCommand* New() {
    return new NextTestCommand();
  }
  vtkTypeMacro(NextTestCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Next test..." << std::endl;
      renderer_->nextTest();
    }
  }

protected:
  NextTestCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class PreviousTestCommand : public vtkCommand {
public:
  static PreviousTestCommand* New() {
    return new PreviousTestCommand();
  }
  vtkTypeMacro(PreviousTestCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Previous test..." << std::endl;
      renderer_->previousTest();
    }
  }

protected:
  PreviousTestCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class AllFlightCommand : public vtkCommand {
public:
  static AllFlightCommand* New() {
    return new AllFlightCommand();
  }
  vtkTypeMacro(AllFlightCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Show all flight..." << std::endl;
      renderer_->showAllFlight();
    }
  }

protected:
  AllFlightCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};

class PlaybackCommand : public vtkCommand {
public:
  static PlaybackCommand* New() {
    return new PlaybackCommand();
  }
  vtkTypeMacro(PlaybackCommand, vtkCommand);

  void SetRenderer(Renderer* renderer) { renderer_ = renderer; }

  void Execute(vtkObject* caller, unsigned long eventId, void* callData) override {
    if (renderer_) {
      std::cout << "Toggle playback animation..." << std::endl;
      renderer_->togglePlaybackAnimation();
    }
  }

protected:
  PlaybackCommand() : renderer_(nullptr) {}

private:
  Renderer* renderer_;
};


void Renderer::initialize() {
  // Create a renderer and render window interactor
  vtkNew<vtkRenderer> renderer;
  renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(1080, 900);

  renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

  // Set the interactor style to extended trackball camera
  vtkNew<CustomInteractorStyle> interactorStyle;
  interactorStyle->SetRenderer(this);
  renderWindowInteractor->SetInteractorStyle(interactorStyle);
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Create command objects for the callbacks
  vtkNew<NextModelCommand> nextModelCommand;
  nextModelCommand->SetRenderer(this);
  vtkNew<PreviousModelCommand> previousModelCommand;
  previousModelCommand->SetRenderer(this);
  vtkNew<NewestModelCommand> newestModelCommand;
  newestModelCommand->SetRenderer(this);
  vtkNew<OldestModelCommand> oldestModelCommand;
  oldestModelCommand->SetRenderer(this);
  vtkNew<NextTestCommand> nextTestCommand;
  nextTestCommand->SetRenderer(this);
  vtkNew<PreviousTestCommand> previousTestCommand;
  previousTestCommand->SetRenderer(this);
  vtkNew<AllFlightCommand> allFlightCommand;
  allFlightCommand->SetRenderer(this);
  vtkNew<PlaybackCommand> playbackCommand;
  playbackCommand->SetRenderer(this);

  // Add observers for the custom events
  interactorStyle->AddObserver(NextModelEvent, nextModelCommand);
  interactorStyle->AddObserver(PreviousModelEvent, previousModelCommand);
  interactorStyle->AddObserver(NewestModelEvent, newestModelCommand);
  interactorStyle->AddObserver(OldestModelEvent, oldestModelCommand);
  interactorStyle->AddObserver(NextTestEvent, nextTestCommand);
  interactorStyle->AddObserver(PreviousTestEvent, previousTestCommand);
  interactorStyle->AddObserver(AllFlightEvent, allFlightCommand);
  interactorStyle->AddObserver(PlaybackEvent, playbackCommand);

  // Add window resize observer
  vtkNew<WindowResizeCommand> resizeCommand;
  resizeCommand->SetRenderer(this);
  renderWindow->AddObserver(vtkCommand::WindowResizeEvent, resizeCommand);

  // Configure the camera
  vtkNew<vtkCamera> camera;
  camera->SetPosition(-50, 0, -2);        // behind the action
  camera->SetFocalPoint(0, 0, -10);       // TODO center this on arenas
  camera->SetViewUp(0, 0, -1);           // Set the view up vector
  camera->SetViewAngle(60);             // Set the field of view (FOV) in degrees

  // Apply the camera settings to the renderer
  renderer->SetActiveCamera(camera);
  renderer->ResetCameraClippingRange(); // Adjust clipping range based on the scene

  // Create the orientation marker widget
  vtkNew<vtkAxesActor> axes;
  orientationWidget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  orientationWidget->SetOutlineColor(0.9300, 0.5700, 0.1300);
  orientationWidget->SetOrientationMarker(axes);
  orientationWidget->SetInteractor(renderWindowInteractor);
  orientationWidget->SetViewport(0.0, 0.0, 0.2, 0.2);
  orientationWidget->SetEnabled(1);
  orientationWidget->InteractiveOn();

  vtkNew<vtkNamedColors> colors;

  actor1 = vtkSmartPointer<vtkActor>::New();
  actor1->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
  actor1->GetProperty()->SetPointSize(4);

  actor2 = vtkSmartPointer<vtkActor>::New();
  actor2->GetProperty()->SetColor(colors->GetColor3d("Yellow").GetData());
  actor2->GetProperty()->SetPointSize(4);

  actor3 = vtkSmartPointer<vtkActor>::New();
  actor3->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
  actor3->GetProperty()->SetPointSize(2);

  renderer->SetBackground(colors->GetColor3d("Black").GetData());

  // Create empty append poly data filters
  paths = vtkSmartPointer<vtkAppendPolyData>::New();
  actuals = vtkSmartPointer<vtkAppendPolyData>::New();
  segmentGaps = vtkSmartPointer<vtkAppendPolyData>::New();
  planeData = vtkSmartPointer<vtkAppendPolyData>::New();
  blackboxTapes = vtkSmartPointer<vtkAppendPolyData>::New();
  blackboxHighlightTapes = vtkSmartPointer<vtkAppendPolyData>::New();

  // Temporary until first update
  vtkNew<vtkPolyData> emptyPolyData;
  vtkNew<vtkPoints> emptyPoints;
  emptyPolyData->SetPoints(emptyPoints);

  // Add the empty polydata to elements
  paths->AddInputData(emptyPolyData);
  actuals->AddInputData(emptyPolyData);
  segmentGaps->AddInputData(emptyPolyData);
  planeData->AddInputData(emptyPolyData);
  blackboxTapes->AddInputData(emptyPolyData);
  blackboxHighlightTapes->AddInputData(emptyPolyData);

  // Update planeData to ensure it has an input port
  planeData->Update();

  // Create an actor for the plane
  vtkNew<vtkPolyDataMapper> planeMapper;
  planeMapper->SetInputConnection(planeData->GetOutputPort());
  vtkNew<vtkActor> planeActor;
  planeActor->SetMapper(planeMapper);

  // Update mappers
  vtkNew<vtkPolyDataMapper> mapper1;
  vtkNew<vtkTubeFilter> tubeFilter1;
  tubeFilter1->SetInputConnection(paths->GetOutputPort());
  tubeFilter1->SetRadius(0.2);
  tubeFilter1->SetNumberOfSides(10);
  mapper1->SetInputConnection(tubeFilter1->GetOutputPort());

  vtkNew<vtkPolyDataMapper> mapper2;
  mapper2->SetInputConnection(actuals->GetOutputPort());

  vtkNew<vtkPolyDataMapper> mapper3;
  vtkNew<vtkTubeFilter> tubeFilter3;
  tubeFilter3->SetInputConnection(segmentGaps->GetOutputPort());
  tubeFilter3->SetRadius(0.1);
  tubeFilter3->SetNumberOfSides(10);
  mapper3->SetInputConnection(tubeFilter3->GetOutputPort());

  // Create blackbox actor (dimmed for full flight display)
  blackboxActor = vtkSmartPointer<vtkActor>::New();
  vtkNew<vtkPolyDataMapper> blackboxMapper;
  blackboxMapper->SetInputConnection(blackboxTapes->GetOutputPort());
  blackboxActor->SetMapper(blackboxMapper);
  
  // Set properties for the blackbox tape (top/bottom coloring)
  vtkProperty* blackboxProperty = blackboxActor->GetProperty();
  blackboxProperty->SetColor(0.0, 1.0, 0.0);  // Green for top face
  blackboxProperty->SetLighting(true);
  blackboxProperty->SetInterpolation(VTK_FLAT);
  blackboxProperty->SetBackfaceCulling(false);
  blackboxProperty->SetFrontfaceCulling(false);
  blackboxProperty->SetAmbient(0.1);
  blackboxProperty->SetDiffuse(0.8);
  blackboxProperty->SetSpecular(0.1);
  blackboxProperty->SetSpecularPower(10);
  blackboxProperty->SetOpacity(1.0);
  
  // Create a new property for the back face (bottom)
  vtkNew<vtkProperty> blackboxBackProperty;
  blackboxBackProperty->SetColor(1.0, 0.5, 0.0);  // Orange for bottom face
  blackboxBackProperty->SetAmbient(0.1);
  blackboxBackProperty->SetDiffuse(0.8);
  blackboxBackProperty->SetSpecular(0.1);
  blackboxBackProperty->SetSpecularPower(10);
  blackboxBackProperty->SetOpacity(1.0);
  
  // Set the back face property
  blackboxActor->SetBackfaceProperty(blackboxBackProperty);
  
  // Create blackbox highlight actor (bright for test spans)
  blackboxHighlightActor = vtkSmartPointer<vtkActor>::New();
  vtkNew<vtkPolyDataMapper> blackboxHighlightMapper;
  blackboxHighlightMapper->SetInputConnection(blackboxHighlightTapes->GetOutputPort());
  blackboxHighlightActor->SetMapper(blackboxHighlightMapper);
  
  // Set properties for highlighted test spans (brighter)
  vtkProperty* highlightProperty = blackboxHighlightActor->GetProperty();
  highlightProperty->SetColor(0.0, 1.0, 0.0);  // Bright green for top face
  highlightProperty->SetLighting(true);
  highlightProperty->SetInterpolation(VTK_FLAT);
  highlightProperty->SetBackfaceCulling(false);
  highlightProperty->SetFrontfaceCulling(false);
  highlightProperty->SetAmbient(0.2);
  highlightProperty->SetDiffuse(1.0);
  highlightProperty->SetSpecular(0.2);
  highlightProperty->SetSpecularPower(20);
  highlightProperty->SetOpacity(1.0);
  
  // Create a new property for the back face (bottom) of highlights
  vtkNew<vtkProperty> highlightBackProperty;
  highlightBackProperty->SetColor(1.0, 0.6, 0.0);  // Bright orange for bottom face
  highlightBackProperty->SetAmbient(0.2);
  highlightBackProperty->SetDiffuse(1.0);
  highlightBackProperty->SetSpecular(0.2);
  highlightBackProperty->SetSpecularPower(20);
  highlightBackProperty->SetOpacity(1.0);
  
  // Set the back face property for highlights
  blackboxHighlightActor->SetBackfaceProperty(highlightBackProperty);

  // Update actors
  actor1->SetMapper(mapper1);
  actor2->SetMapper(mapper2);
  actor3->SetMapper(mapper3);

  // Set properties for the tape (actor2)
  vtkProperty* property = actor2->GetProperty();

  // Enable double-sided surface rendering
  property->SetLighting(true);
  property->SetInterpolation(VTK_FLAT);
  property->SetBackfaceCulling(false);
  property->SetFrontfaceCulling(false);

  // Set different colors for front and back faces
  property->SetColor(1.0, 0.7, 0.0);  // Orange for front face (top)
  property->SetAmbient(0.1);
  property->SetDiffuse(0.8);
  property->SetSpecular(0.1);
  property->SetSpecularPower(10);

  // Set opacity to 1 (fully opaque)
  property->SetOpacity(1.0);

  // Create a new property for the back face
  vtkNew<vtkProperty> backProperty;
  backProperty->SetColor(0.0, 1.0, 1.0);  // Blue for back face (bottom)
  backProperty->SetAmbient(0.1);
  backProperty->SetDiffuse(0.8);
  backProperty->SetSpecular(0.1);
  backProperty->SetSpecularPower(10);
  backProperty->SetOpacity(1.0);

  // Set the back face property
  actor2->SetBackfaceProperty(backProperty);

  renderer->AddActor(planeActor);
  renderer->AddActor(actor1);
  renderer->AddActor(actor2);
  renderer->AddActor(actor3);
  
  // Only add blackbox actors if there's blackbox data
  if (!blackboxAircraftStates.empty()) {
    renderer->AddActor(blackboxActor);
    renderer->AddActor(blackboxHighlightActor);
  }

  // Create text actors for generation and fitness display
  generationTextActor = vtkSmartPointer<vtkTextActor>::New();
  generationTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // White text
  generationTextActor->GetTextProperty()->SetVerticalJustificationToBottom();
  generationTextActor->GetTextProperty()->SetJustificationToRight();
  generationTextActor->SetInput("Generation:");
  renderer->AddActor2D(generationTextActor);

  generationValueActor = vtkSmartPointer<vtkTextActor>::New();
  generationValueActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // White text
  generationValueActor->GetTextProperty()->SetVerticalJustificationToBottom();
  generationValueActor->GetTextProperty()->SetJustificationToLeft();
  generationValueActor->SetInput("0");
  renderer->AddActor2D(generationValueActor);

  fitnessTextActor = vtkSmartPointer<vtkTextActor>::New();
  fitnessTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // White text
  fitnessTextActor->GetTextProperty()->SetVerticalJustificationToBottom();
  fitnessTextActor->GetTextProperty()->SetJustificationToRight();
  fitnessTextActor->SetInput("Fitness:");
  renderer->AddActor2D(fitnessTextActor);

  fitnessValueActor = vtkSmartPointer<vtkTextActor>::New();
  fitnessValueActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // White text
  fitnessValueActor->GetTextProperty()->SetVerticalJustificationToBottom();
  fitnessValueActor->GetTextProperty()->SetJustificationToLeft();
  fitnessValueActor->SetInput("0.000");
  renderer->AddActor2D(fitnessValueActor);

  testTextActor = vtkSmartPointer<vtkTextActor>::New();
  testTextActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // White text
  testTextActor->GetTextProperty()->SetVerticalJustificationToBottom();
  testTextActor->GetTextProperty()->SetJustificationToRight();
  testTextActor->SetInput("Test:");
  renderer->AddActor2D(testTextActor);

  testValueActor = vtkSmartPointer<vtkTextActor>::New();
  testValueActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0); // White text
  testValueActor->GetTextProperty()->SetVerticalJustificationToBottom();
  testValueActor->GetTextProperty()->SetJustificationToLeft();
  testValueActor->SetInput("-");
  renderer->AddActor2D(testValueActor);

  // Set initial text display
  updateTextDisplay(0, 0.0);

  // Enable anti-aliasing (multi-sampling)
  renderWindow->SetMultiSamples(4); // Use 4x MSAA

  // Enable depth peeling for proper transparency rendering
  renderWindow->SetAlphaBitPlanes(1);
  renderWindow->SetMultiSamples(0);
  renderer->SetUseDepthPeeling(1);
  renderer->SetMaximumNumberOfPeels(100);  // Maximum number of depth peels
  renderer->SetOcclusionRatio(0.1);        // Occlusion ratio

  // render
  renderWindow->Render();

  // exit command
  vtkNew<ExitCommand> exitCommand;
  exitCommand->Interactor = renderWindowInteractor;
  renderWindowInteractor->AddObserver(vtkCommand::ExitEvent, exitCommand);
  
  // Timer callback for animation
  vtkNew<vtkCallbackCommand> timerCallback;
  timerCallback->SetCallback([](vtkObject* caller, unsigned long eid, void* clientdata, void* calldata) {
    Renderer* renderer = static_cast<Renderer*>(clientdata);
    // Double-check that animation is still active and timer ID is valid
    if (renderer && renderer->isPlaybackActive && renderer->animationTimerId != 0) {
      // Additional safety check - only proceed if we have a valid renderWindowInteractor
      if (renderer->renderWindowInteractor) {
        renderer->updatePlaybackAnimation();
      }
    }
  });
  timerCallback->SetClientData(this);
  renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, timerCallback);
}

std::vector<Eigen::Vector3d> Renderer::pathToVector(std::vector<Path> path) {
  std::vector<Eigen::Vector3d> points;
  for (const auto& p : path) {
    points.push_back(p.start);
  }
  return points;
}

std::vector<Eigen::Vector3d> Renderer::stateToVector(std::vector<AircraftState> state) {
  std::vector<Eigen::Vector3d> points;
  for (const auto& s : state) {
    points.push_back({ s.getPosition()[0], s.getPosition()[1], s.getPosition()[2] });
  }
  return points;
}

std::vector<Eigen::Vector3d> Renderer::stateToOrientation(std::vector<AircraftState> state) {
  std::vector<Eigen::Vector3d> points;
  for (const auto& s : state) {
    // Use the aircraft body Z-axis (up direction) for ribbon orientation
    points.push_back(s.getOrientation() * -Eigen::Vector3d::UnitZ());
  }
  return points;
}

double Renderer::extractFitnessFromGP(const std::vector<char>& gpData) {
  if (gpData.empty()) {
    return 0.0;
  }
  
  try {
    // Create stream from the char vector
    boost::iostreams::stream<boost::iostreams::array_source> inStream(gpData.data(), gpData.size());
    
    // Create and load a base GP object
    GP gp;
    gp.load(inStream);
    
    return gp.getFitness();
  }
  catch (const std::exception& e) {
    std::cerr << "Error extracting fitness from GP: " << e.what() << std::endl;
    return 0.0;
  }
}

void Renderer::updateTextDisplay(int generation, double fitness) {
  // Store current values for resize updates
  currentGeneration = generation;
  currentFitness = fitness;
  
  // Get window size for responsive positioning
  int* size = renderWindow->GetSize();
  int windowWidth = size[0];
  int windowHeight = size[1];
  
  // Calculate font size based on window size (min 8, max 20)
  int fontSize = std::max(8, std::min(20, static_cast<int>(windowHeight * 0.02)));
  
  // Position calculations: bottom-right of text at 95% width and 95% height
  int anchorX = static_cast<int>(windowWidth * 0.95);
  // In VTK, Y=0 is at bottom, so 95% height means 5% from bottom
  int anchorY = static_cast<int>(windowHeight * 0.05); // 5% from bottom = lower right
  
  // Calculate line spacing
  int lineSpacing = fontSize + 6;
  
  // Position text so bottom-right corner is at anchor point
  // In VTK coordinates: test is lowest, fitness is middle, generation is highest
  int testY = anchorY;
  int fitnessY = anchorY + lineSpacing;
  int generationY = anchorY + (2 * lineSpacing);
  
  // Labels and values positioned relative to anchor
  int labelX = anchorX - 80; // Labels positioned left of anchor for right-justification
  int valueX = anchorX - 75; // Values start slightly right of labels
  
  // Update font sizes for all text actors
  generationTextActor->GetTextProperty()->SetFontSize(fontSize);
  generationValueActor->GetTextProperty()->SetFontSize(fontSize);
  fitnessTextActor->GetTextProperty()->SetFontSize(fontSize);
  fitnessValueActor->GetTextProperty()->SetFontSize(fontSize);
  testTextActor->GetTextProperty()->SetFontSize(fontSize);
  testValueActor->GetTextProperty()->SetFontSize(fontSize);
  
  // Update positions - labels right-justified, values left-justified
  generationTextActor->SetPosition(labelX, generationY);
  generationValueActor->SetPosition(valueX, generationY);
  fitnessTextActor->SetPosition(labelX, fitnessY);
  fitnessValueActor->SetPosition(valueX, fitnessY);
  testTextActor->SetPosition(labelX, testY);
  testValueActor->SetPosition(valueX, testY);
  
  // Update text content - use same formula as title bar (10000 - generation)
  int realGeneration = 10000 - generation;
  generationValueActor->SetInput(std::to_string(realGeneration).c_str());
  
  std::ostringstream fitnessStream;
  fitnessStream << std::fixed << std::setprecision(3) << fitness;
  fitnessValueActor->SetInput(fitnessStream.str().c_str());
  
  // Update test display based on decode mode and current state
  if (inDecodeMode && !testSpans.empty()) {
    // Make test actors visible
    testTextActor->SetVisibility(1);
    testValueActor->SetVisibility(1);
    
    if (showingFullFlight) {
      testValueActor->SetInput("All");
    } else if (currentTestIndex >= 0 && currentTestIndex < testSpans.size()) {
      testValueActor->SetInput(std::to_string(currentTestIndex + 1).c_str());
    } else {
      testValueActor->SetInput("-");
    }
  } else {
    // Hide test actors when not in decode mode or no test spans
    testTextActor->SetVisibility(0);
    testValueActor->SetVisibility(0);
  }
}

// Extract the generation number from the key
int extractGenNumber(const std::string& input) {
  std::regex pattern("autoc-.*\\/gen(\\d+)\\.dmp");
  std::smatch matches;

  if (std::regex_search(input, matches, pattern) && matches.size() > 1) {
    return std::stoi(matches[1].str());
  }

  return -1;
}

// main method
int main(int argc, char** argv) {
  // Parse command line arguments
  static struct option long_options[] = {
    {"keyname", required_argument, 0, 'k'},
    {"decoder", required_argument, 0, 'd'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  int option_index = 0;
  int c;
  
  while ((c = getopt_long(argc, argv, "k:d:h", long_options, &option_index)) != -1) {
    switch (c) {
      case 'k':
        computedKeyName = optarg;
        break;
      case 'd':
        decoderCommand = optarg;
        break;
      case 'h':
        printUsage(argv[0]);
        return 0;
      case '?':
        printUsage(argv[0]);
        return 1;
      default:
        break;
    }
  }
  
  // Handle positional arguments for backward compatibility
  if (computedKeyName.empty() && optind < argc) {
    computedKeyName = argv[optind];
  }
  
  std::string keyName = "";
  
  // Initialize configuration
  ConfigManager::initialize("autoc.ini");

  // AWS setup
  Aws::SDKOptions options;
  Aws::InitAPI(options);

  std::shared_ptr<Aws::S3::S3Client> s3_client = getS3Client();
  
  // Load blackbox data if decoder command is specified
  if (!decoderCommand.empty()) {
    if (!loadBlackboxData()) {
      std::cerr << "Failed to load blackbox data" << std::endl;
      return 1;
    }
    if (!blackboxAircraftStates.empty()) {
      std::cout << "Blackbox data loaded: " << blackboxAircraftStates.size() << " states" << std::endl;
      renderer.inDecodeMode = true;
      renderer.extractTestSpans();
      
      if (!renderer.testSpans.empty()) {
        std::cout << "Found " << renderer.testSpans.size() << " test spans" << std::endl;
        renderer.currentTestIndex = 0;
        renderer.showingFullFlight = false;
        // Filter blackbox data to show first test span immediately
        updateBlackboxForCurrentTest();
        std::cout << "Showing test 1: " << blackboxAircraftStates.size() << " states" << std::endl;
      } else {
        std::cout << "No MSPRCOVERRIDE test spans found, showing full flight" << std::endl;
        renderer.showingFullFlight = true;
      }
    }
  }

  // should we look up the latest run?
  if (computedKeyName.empty()) {
    Aws::S3::Model::ListObjectsV2Request listFolders;
    listFolders.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
    listFolders.SetPrefix("autoc-");
    listFolders.SetDelimiter("/");

    bool isTruncated = false;
    do {
      auto outcome = s3_client->ListObjectsV2(listFolders);
      if (outcome.IsSuccess()) {
        const auto& result = outcome.GetResult();

        // Process common prefixes (these are our 'folders')
        for (const auto& commonPrefix : result.GetCommonPrefixes()) {
          computedKeyName = commonPrefix.GetPrefix();
          break;
        }

        // Check if there are more results to retrieve
        isTruncated = result.GetIsTruncated();
        if (isTruncated) {
          listFolders.SetContinuationToken(result.GetNextContinuationToken());
        }
      }
      else {
        std::cerr << "Error listing objects: " << outcome.GetError().GetMessage() << std::endl;
        break;
      }
    } while (isTruncated);
  }

  // ok for this run, look for the last generation
  Aws::S3::Model::ListObjectsV2Request listItem;
  listItem.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
  listItem.SetPrefix(computedKeyName + "gen");
  bool isTruncated = false;
  do {
    auto outcome = s3_client->ListObjectsV2(listItem);
    if (outcome.IsSuccess()) {
      // S3 returns objects in lexicographical order (gen9900, gen9901, gen9902...)
      // Due to reverse numbering (10000-gen), the FIRST object is the LATEST generation
      for (const auto& object : outcome.GetResult().GetContents()) {
        keyName = object.GetKey();
        break; // Take the first (latest) generation
      }

      // Check if the response is truncated
      isTruncated = outcome.GetResult().GetIsTruncated();
      if (isTruncated) {
        // Set the continuation token for the next request
        listItem.SetContinuationToken(outcome.GetResult().GetNextContinuationToken());
      }
    }
    else {
      std::cerr << "Error listing objects: " << outcome.GetError().GetMessage() << std::endl;
      break;
    }
  } while (isTruncated);

  // Create a renderer
  renderer.initialize();

  // display initial data
  renderer.genNumber = extractGenNumber(keyName);
  renderer.updateGenerationDisplay(renderer.genNumber);
  
  // Print keyboard controls
  std::cout << "\nKeyboard Controls:" << std::endl;
  std::cout << "  n - Next generation" << std::endl;
  std::cout << "  p - Previous generation" << std::endl;
  std::cout << "  N - Jump to newest generation" << std::endl;
  std::cout << "  P - Jump to oldest generation (generation 1)" << std::endl;
  std::cout << "  SPACE - Toggle playback animation" << std::endl;
  if (!decoderCommand.empty() && !renderer.testSpans.empty()) {
    std::cout << "  t - Next test segment" << std::endl;
    std::cout << "  r - Previous test segment" << std::endl;
    std::cout << "  a - Show all flight (with test highlights)" << std::endl;
  }
  std::cout << "  q - Quit" << std::endl;

  // kick it off
  renderer.renderWindowInteractor->Start();

  return 0;
}

// Parse blackbox decode CSV output
bool parseBlackboxData(const std::string& csvData) {
  std::istringstream stream(csvData);
  std::string line;
  std::vector<std::string> headers;
  bool headerParsed = false;
  
  blackboxPoints.clear();
  blackboxNormals.clear();
  blackboxAircraftStates.clear();
  fullBlackboxAircraftStates.clear();
  csvLines.clear();
  
  // Find column indices
  int latIndex = -1, lonIndex = -1, altIndex = -1, timeIndex = -1;
  int quatWIndex = -1, quatXIndex = -1, quatYIndex = -1, quatZIndex = -1;
  
  while (std::getline(stream, line)) {
    csvLines.push_back(line);  // Store all lines for span analysis
    
    if (line.empty() || line.find("End of log") != std::string::npos) {
      continue;
    }
    
    // Parse header
    if (!headerParsed) {
      std::istringstream headerStream(line);
      std::string header;
      int index = 0;
      
      while (std::getline(headerStream, header, ',')) {
        // Trim whitespace from header
        header.erase(0, header.find_first_not_of(" \t"));
        header.erase(header.find_last_not_of(" \t") + 1);
        
        headers.push_back(header);
        if (header == "navPos[0]") latIndex = index;
        else if (header == "navPos[1]") lonIndex = index;
        else if (header == "navPos[2]") altIndex = index;
        else if (header == "time (us)") timeIndex = index;
        else if (header == "quaternion[0]") quatWIndex = index;
        else if (header == "quaternion[1]") quatXIndex = index;
        else if (header == "quaternion[2]") quatYIndex = index;
        else if (header == "quaternion[3]") quatZIndex = index;
        index++;
      }
      headerParsed = true;
      continue;
    }
    
    // Parse data rows
    std::istringstream dataStream(line);
    std::string cell;
    std::vector<std::string> row;
    
    while (std::getline(dataStream, cell, ',')) {
      // Trim whitespace from cell
      cell.erase(0, cell.find_first_not_of(" \t"));
      cell.erase(cell.find_last_not_of(" \t") + 1);
      row.push_back(cell);
    }
    
    if (row.size() < headers.size()) continue;
    
    // Sample at 50ms intervals (50,000 microseconds)
    static unsigned long int lastSampleTime = 0;
    const unsigned long int SAMPLE_INTERVAL_US = 50000; // 50ms in microseconds
    
    if (timeIndex >= 0) {
      unsigned long int currentTime = std::stoul(row[timeIndex]);
      if (lastSampleTime > 0 && (currentTime - lastSampleTime) < SAMPLE_INTERVAL_US) {
        continue; // Skip this record, not enough time has passed
      }
      lastSampleTime = currentTime;
    }
    
    // Extract navPos coordinates and convert to meters
    if (latIndex >= 0 && lonIndex >= 0 && altIndex >= 0) {
      double navX = std::stod(row[latIndex]);  // navPos[0] in cm
      double navY = std::stod(row[lonIndex]);  // navPos[1] in cm
      double navZ = std::stod(row[altIndex]);  // navPos[2] in cm
      
      // Convert from centimeters to meters and center on origin
      static bool firstPoint = true;
      static double originX, originY, originZ;
      
      if (firstPoint) {
        originX = navX;
        originY = navY;
        originZ = navZ;
        firstPoint = false;
      }
      
      // Convert to meters and translate to origin
      // NED coordinates: North=X, East=Y, Down=Z
      // For VTK visualization, flip Z to make "up" positive (NED Down -> VTK Up)
      double x = (navX - originX) / 100.0;   // cm to m, North
      double y = (navY - originY) / 100.0;   // cm to m, East  
      double z = -(navZ - originZ) / 100.0;  // cm to m, Down->Up (flip sign)
      
      Eigen::Vector3d newPoint(x, y, z);
      
      // Filter out coincident points to avoid VTK ribbon filter issues
      if (blackboxPoints.empty() || (newPoint - blackboxPoints.back()).norm() > 0.01) {
        
        blackboxPoints.push_back(newPoint);
        
        // Create AircraftState object
        if (quatWIndex >= 0 && quatXIndex >= 0 && quatYIndex >= 0 && quatZIndex >= 0) {
          // Parse normalized quaternion values (stored as integers * 10000)
          double qw = std::stod(row[quatWIndex]) / 10000.0;
          double qx = std::stod(row[quatXIndex]) / 10000.0;
          double qy = std::stod(row[quatYIndex]) / 10000.0;
          double qz = std::stod(row[quatZIndex]) / 10000.0;
          
          // Create quaternion directly from normalized values
          Eigen::Quaterniond q(qw, qx, qy, qz);
          
          // Ensure quaternion is normalized (should already be, but safety check)
          q.normalize();
          
          // Get time if available, otherwise use index
          unsigned long int timeMs = (timeIndex >= 0) ? std::stoul(row[timeIndex]) : blackboxAircraftStates.size() * 100;
          
          // Create simple velocity vector for blackbox data (assuming forward flight)
          Eigen::Vector3d velocity_vector = q * Eigen::Vector3d(20.0, 0, 0);
          
          // Create AircraftState with blackbox data
          AircraftState state(blackboxAircraftStates.size(), 20.0, velocity_vector, q, newPoint, 0.0, 0.0, 0.0, timeMs);
          blackboxAircraftStates.push_back(state);
          fullBlackboxAircraftStates.push_back(state);  // Store full flight data
          
        } else {
          // Error: quaternion information is required for proper tape orientation
          std::cerr << "Error: Blackbox data missing quaternion information (quaternion[0], quaternion[1], quaternion[2], quaternion[3])" << std::endl;
          std::cerr << "Cannot render tape with proper orientation without aircraft quaternion data" << std::endl;
          return false;
        }
      }
    }
  }
  
  return !blackboxAircraftStates.empty();
}

// Method to update blackbox rendering based on current test selection
void updateBlackboxForCurrentTest() {
  extern Renderer renderer;
  
  if (!renderer.inDecodeMode || fullBlackboxAircraftStates.empty()) {
    return;
  }
  
  // Always update the blackbox data, even if renderer isn't fully initialized yet
  // The safety check for renderWindow is only needed for text display updates
  
  if (renderer.showingFullFlight || renderer.testSpans.empty()) {
    // Show full flight
    blackboxAircraftStates = fullBlackboxAircraftStates;
  } else if (renderer.currentTestIndex >= 0 && renderer.currentTestIndex < renderer.testSpans.size()) {
    // Show only the selected test span
    const TestSpan& span = renderer.testSpans[renderer.currentTestIndex];
    blackboxAircraftStates.clear();
    
    size_t startIdx = std::min(span.startIndex, fullBlackboxAircraftStates.size());
    size_t endIdx = std::min(span.endIndex + 1, fullBlackboxAircraftStates.size());
    
    // Ensure we have valid indices
    if (startIdx < endIdx && startIdx < fullBlackboxAircraftStates.size()) {
      for (size_t i = startIdx; i < endIdx; i++) {
        blackboxAircraftStates.push_back(fullBlackboxAircraftStates[i]);
      }
      
      std::cout << "Showing test " << (renderer.currentTestIndex + 1) << ": " 
                << blackboxAircraftStates.size() << " states (indices " 
                << startIdx << "-" << (endIdx-1) << ")" << std::endl;
    } else {
      std::cerr << "Invalid test span indices: " << startIdx << "-" << endIdx << std::endl;
      // Fallback to full flight
      blackboxAircraftStates = fullBlackboxAircraftStates;
      renderer.showingFullFlight = true;
    }
  }
  
  // Update text display (only if renderer is initialized)
  if (renderer.renderWindow) {
    renderer.updateTextDisplay(renderer.currentGeneration, renderer.currentFitness);
  }
}


// Load CSV data from decoder command or stdin
bool loadBlackboxData() {
  std::string csvData;
  
  if (decoderCommand == "-") {
    // Read from stdin
    std::ostringstream result;
    std::string line;
    while (std::getline(std::cin, line)) {
      result << line << "\n";
    }
    csvData = result.str();
  } else {
    // Execute the decoder command
    FILE* pipe = popen(decoderCommand.c_str(), "r");
    if (!pipe) {
      std::cerr << "Failed to execute decoder command: " << decoderCommand << std::endl;
      return false;
    }
    
    std::ostringstream result;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
      result << buffer;
    }
    
    int exitCode = pclose(pipe);
    if (exitCode != 0) {
      std::cerr << "Decoder command failed with exit code: " << exitCode << std::endl;
      return false;
    }
    
    csvData = result.str();
  }
  
  return parseBlackboxData(csvData);
}

void Renderer::jumpToNewestGeneration() {
  // The newest generation has the lowest generation number (because of 10000-gen numbering)
  // We need to find the lowest numbered generation file in S3
  
  extern std::string computedKeyName;
  auto s3_client = getS3Client();
  
  Aws::S3::Model::ListObjectsV2Request listItem;
  listItem.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
  listItem.SetPrefix(computedKeyName + "gen");
  
  std::string newestKeyName = "";
  int lowestGenNumber = 10000; // Start with highest possible
  
  bool isTruncated = false;
  do {
    auto outcome = s3_client->ListObjectsV2(listItem);
    if (outcome.IsSuccess()) {
      for (const auto& object : outcome.GetResult().GetContents()) {
        int genNum = extractGenNumber(object.GetKey());
        if (genNum != -1 && genNum < lowestGenNumber) {
          lowestGenNumber = genNum;
          newestKeyName = object.GetKey();
        }
      }
      
      isTruncated = outcome.GetResult().GetIsTruncated();
      if (isTruncated) {
        listItem.SetContinuationToken(outcome.GetResult().GetNextContinuationToken());
      }
    }
    else {
      std::cerr << "Error listing objects: " << outcome.GetError().GetMessage() << std::endl;
      return;
    }
  } while (isTruncated);
  
  if (!newestKeyName.empty()) {
    genNumber = lowestGenNumber;
    updateGenerationDisplay(genNumber);
  }
}

void Renderer::jumpToOldestGeneration() {
  // Generation 1 is stored as gen9999.dmp
  genNumber = 9999;
  updateGenerationDisplay(genNumber);
}

// Print usage information
void Renderer::extractTestSpans() {
  testSpans.clear();
  
  if (csvLines.empty() || fullBlackboxAircraftStates.empty()) {
    std::cerr << "No CSV data or aircraft states available for test span extraction" << std::endl;
    return;
  }
  
  // Find the flightModeFlags column index
  int flightModeFlagsIndex = -1;
  int timeIndex = -1;
  
  // Parse header to find column indices
  if (!csvLines.empty()) {
    std::istringstream headerStream(csvLines[0]);
    std::string header;
    int index = 0;
    
    while (std::getline(headerStream, header, ',')) {
      // Trim whitespace from header
      header.erase(0, header.find_first_not_of(" \t"));
      header.erase(header.find_last_not_of(" \t") + 1);
      
      if (header == "flightModeFlags (flags)") {
        flightModeFlagsIndex = index;
      } else if (header == "time (us)" || header == "time (us) ") {
        timeIndex = index;
      }
      index++;
    }
  }
  
  if (flightModeFlagsIndex == -1) {
    std::cerr << "flightModeFlags column not found in CSV" << std::endl;
    return;
  }
  
  // Scan through data rows to find MSPRCOVERRIDE spans
  // We need to correlate CSV rows with aircraft state indices based on sampling
  bool inSpan = false;
  TestSpan currentSpan;
  size_t aircraftStateIndex = 0;
  
  for (size_t lineIndex = 1; lineIndex < csvLines.size(); lineIndex++) {
    const std::string& line = csvLines[lineIndex];
    if (line.empty() || line.find("End of log") != std::string::npos) {
      continue;
    }
    
    // Parse data row
    std::istringstream dataStream(line);
    std::string cell;
    std::vector<std::string> row;
    
    while (std::getline(dataStream, cell, ',')) {
      // Trim whitespace from cell
      cell.erase(0, cell.find_first_not_of(" \t"));
      cell.erase(cell.find_last_not_of(" \t") + 1);
      row.push_back(cell);
    }
    
    if (row.size() <= flightModeFlagsIndex) continue;
    
    // Check if this row would have been sampled (same logic as parseBlackboxData)
    static unsigned long int lastSampleTime = 0;
    const unsigned long int SAMPLE_INTERVAL_US = 50000; // 50ms in microseconds
    
    bool wouldBeSampled = false;
    if (timeIndex >= 0 && timeIndex < row.size()) {
      unsigned long int currentTime = std::stoul(row[timeIndex]);
      if (lastSampleTime == 0 || (currentTime - lastSampleTime) >= SAMPLE_INTERVAL_US) {
        wouldBeSampled = true;
        lastSampleTime = currentTime;
      }
    } else {
      wouldBeSampled = true; // If no time column, assume all rows are sampled
    }
    
    if (!wouldBeSampled) continue;
    
    // This row corresponds to an aircraft state
    if (aircraftStateIndex >= fullBlackboxAircraftStates.size()) break;
    
    bool hasMSPRCOVERRIDE = row[flightModeFlagsIndex].find(MSPRCOVERRIDE_FLAG) != std::string::npos;
    unsigned long currentTime = (timeIndex >= 0 && timeIndex < row.size()) ? std::stoul(row[timeIndex]) : 0;
    
    if (hasMSPRCOVERRIDE && !inSpan) {
      // Start of new span
      currentSpan.startIndex = aircraftStateIndex;
      currentSpan.startTime = currentTime;
      inSpan = true;
    } else if (!hasMSPRCOVERRIDE && inSpan) {
      // End of current span
      currentSpan.endIndex = aircraftStateIndex - 1;
      currentSpan.endTime = currentTime;
      testSpans.push_back(currentSpan);
      inSpan = false;
    }
    
    aircraftStateIndex++;
  }
  
  // Handle case where span continues to end of file
  if (inSpan) {
    currentSpan.endIndex = aircraftStateIndex - 1;
    currentSpan.endTime = 0; // Use 0 if no time available
    testSpans.push_back(currentSpan);
  }
  
  // Validate all spans have valid indices
  std::vector<TestSpan> validSpans;
  for (const TestSpan& span : testSpans) {
    if (span.startIndex < fullBlackboxAircraftStates.size() && 
        span.endIndex < fullBlackboxAircraftStates.size() && 
        span.startIndex <= span.endIndex) {
      validSpans.push_back(span);
    } else {
      std::cerr << "Invalid span removed: indices " << span.startIndex << "-" << span.endIndex 
                << " (max: " << fullBlackboxAircraftStates.size() - 1 << ")" << std::endl;
    }
  }
  testSpans = validSpans;
  
  std::cout << "Extracted " << testSpans.size() << " valid test spans from " << aircraftStateIndex << " aircraft states" << std::endl;
  for (size_t i = 0; i < testSpans.size(); i++) {
    std::cout << "  Test " << (i+1) << ": indices " << testSpans[i].startIndex << "-" << testSpans[i].endIndex 
              << " (" << (testSpans[i].endIndex - testSpans[i].startIndex + 1) << " states)" << std::endl;
  }
}

void Renderer::nextTest() {
  if (testSpans.empty()) {
    std::cout << "No test spans available" << std::endl;
    return;
  }
  
  if (showingFullFlight) {
    // Switch from full flight to first test
    currentTestIndex = 0;
    showingFullFlight = false;
  } else {
    // Move to next test
    currentTestIndex = (currentTestIndex + 1) % testSpans.size();
  }
  
  // Update the blackbox rendering to show only the selected test span
  updateBlackboxForCurrentTest();
  
  // Re-render using in-memory data (no S3 fetch)
  renderFullScene();
}

void Renderer::previousTest() {
  if (testSpans.empty()) {
    std::cout << "No test spans available" << std::endl;
    return;
  }
  
  if (showingFullFlight) {
    // Switch from full flight to last test
    currentTestIndex = testSpans.size() - 1;
    showingFullFlight = false;
  } else {
    // Move to previous test
    currentTestIndex = (currentTestIndex - 1 + testSpans.size()) % testSpans.size();
  }
  
  // Update the blackbox rendering to show only the selected test span
  updateBlackboxForCurrentTest();
  
  // Re-render using in-memory data (no S3 fetch)
  renderFullScene();
}

void Renderer::showAllFlight() {
  showingFullFlight = true;
  
  // Restore full blackbox data
  blackboxAircraftStates = fullBlackboxAircraftStates;
  
  // Update the blackbox rendering
  updateBlackboxForCurrentTest();
  
  // Re-render using in-memory data (no S3 fetch)
  renderFullScene();
}

void Renderer::createHighlightedFlightTapes(Eigen::Vector3d offset) {
  if (fullBlackboxAircraftStates.empty() || testSpans.empty()) {
    return;
  }
  
  // Remove the empty data we added earlier
  this->blackboxHighlightTapes->RemoveAllInputs();
  
  // Create dimmed version of full flight (25% brightness)
  std::vector<Eigen::Vector3d> fullPoints = stateToVector(fullBlackboxAircraftStates);
  if (fullPoints.size() >= 2) {
    this->blackboxTapes->AddInputData(createTapeSet(offset, fullPoints, stateToOrientation(fullBlackboxAircraftStates)));
    // Adjust opacity to 25% for dimmed effect
    blackboxActor->GetProperty()->SetOpacity(0.25);
    blackboxActor->GetBackfaceProperty()->SetOpacity(0.25);
  }
  
  // Create bright highlighted segments for each test span
  for (const TestSpan& span : testSpans) {
    size_t startIdx = std::min(span.startIndex, fullBlackboxAircraftStates.size());
    size_t endIdx = std::min(span.endIndex + 1, fullBlackboxAircraftStates.size());
    
    if (startIdx < endIdx && startIdx < fullBlackboxAircraftStates.size()) {
      std::vector<AircraftState> spanStates;
      for (size_t i = startIdx; i < endIdx; i++) {
        spanStates.push_back(fullBlackboxAircraftStates[i]);
      }
      
      if (spanStates.size() >= 2) {
        std::vector<Eigen::Vector3d> spanPoints = stateToVector(spanStates);
        this->blackboxHighlightTapes->AddInputData(createTapeSet(offset, spanPoints, stateToOrientation(spanStates)));
      }
    }
  }
  
  // Ensure highlight actor is at full brightness
  blackboxHighlightActor->GetProperty()->SetOpacity(1.0);
  blackboxHighlightActor->GetBackfaceProperty()->SetOpacity(1.0);
}

void Renderer::togglePlaybackAnimation() {
  if (isPlaybackActive) {
    // Stop animation
    isPlaybackActive = false;
    isPlaybackPaused = false;
    totalPausedTime = std::chrono::duration<double>::zero();
    if (animationTimerId != 0) {
      renderWindowInteractor->DestroyTimer(animationTimerId);
      animationTimerId = 0;
    }
    std::cout << "Playback animation stopped" << std::endl;
    // Re-render full scene using existing data (no S3 fetch)
    renderFullScene();
  } else {
    // Start animation
    isPlaybackActive = true;
    isPlaybackPaused = false;
    totalPausedTime = std::chrono::duration<double>::zero();
    animationStartTime = std::chrono::steady_clock::now();
    std::cout << "Real-time playback animation started" << std::endl;
    // Start with first animation frame
    updatePlaybackAnimation();
  }
}

void Renderer::pausePlaybackAnimation() {
  if (isPlaybackActive && !isPlaybackPaused) {
    isPlaybackPaused = true;
    pauseStartTime = std::chrono::steady_clock::now();
    std::cout << "Animation paused for camera interaction" << std::endl;
  }
}

void Renderer::resumePlaybackAnimation() {
  if (isPlaybackActive && isPlaybackPaused) {
    auto pauseDuration = std::chrono::steady_clock::now() - pauseStartTime;
    totalPausedTime += pauseDuration;
    isPlaybackPaused = false;
    std::cout << "Animation resumed after camera interaction" << std::endl;
    // Restart the animation timer
    updatePlaybackAnimation();
  }
}

void Renderer::updatePlaybackAnimation() {
  if (!isPlaybackActive) return;
  
  // Skip updates while paused
  if (isPlaybackPaused) return;
  
  auto currentTime = std::chrono::steady_clock::now();
  auto totalElapsed = currentTime - animationStartTime;
  auto effectiveElapsed = totalElapsed - totalPausedTime;
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(effectiveElapsed).count() / 1000.0;
  
  // Path always drives the timing - find longest simulation duration across all paths
  double primaryDuration = 0.0;
  int maxArenas = evalResults.pathList.size();
  for (int i = 0; i < maxArenas; i++) {
    if (!evalResults.aircraftStateList[i].empty()) {
      double pathMaxTime = evalResults.aircraftStateList[i].back().getSimTimeMsec() / 1000.0;
      primaryDuration = std::max(primaryDuration, pathMaxTime);
    }
  }
  
  // Use real-time playback - elapsed seconds = simulation seconds
  double currentSimTime = elapsed;
  
  // Check if animation is complete
  if (primaryDuration > 0.0 && currentSimTime >= primaryDuration) {
    currentSimTime = primaryDuration;
    isPlaybackActive = false;
    animationTimerId = 0;
    std::cout << "Playback animation completed (duration: " << primaryDuration << "s)" << std::endl;
  }
  
  // Clear existing data
  this->paths->RemoveAllInputs();
  this->actuals->RemoveAllInputs();
  this->segmentGaps->RemoveAllInputs();
  this->blackboxTapes->RemoveAllInputs();
  this->blackboxHighlightTapes->RemoveAllInputs();
  
  // Always ensure highlight tapes has at least empty data
  vtkNew<vtkPolyData> emptyHighlightData;
  vtkNew<vtkPoints> emptyHighlightPoints;
  emptyHighlightData->SetPoints(emptyHighlightPoints);
  this->blackboxHighlightTapes->AddInputData(emptyHighlightData);
  
  // Render with synchronized time progress
  int renderArenas = (!blackboxAircraftStates.empty()) ? 1 : evalResults.pathList.size();
  for (int i = 0; i < renderArenas; i++) {
    Eigen::Vector3d offset = renderingOffset(i);
    
    std::vector<Eigen::Vector3d> p = pathToVector(evalResults.pathList[i]);
    std::vector<Eigen::Vector3d> a = stateToVector(evalResults.aircraftStateList[i]);
    
    // Calculate progress for this path based on real simulation time
    double simProgress = 1.0; // Default to show all
    if (!evalResults.aircraftStateList[i].empty()) {
      double simDuration = evalResults.aircraftStateList[i].back().getSimTimeMsec() / 1000.0;
      if (simDuration > 0.0) {
        simProgress = std::min(1.0, currentSimTime / simDuration);
      }
    }
    
    if (!p.empty()) {
      this->paths->AddInputData(createPointSet(offset, p, simProgress));
    }
    if (!a.empty()) {
      this->actuals->AddInputData(createTapeSet(offset, a, stateToOrientation(evalResults.aircraftStateList[i]), simProgress));
    }
    if (!a.empty() && !p.empty()) {
      this->segmentGaps->AddInputData(createSegmentSet(offset, evalResults.aircraftStateList[i], p, simProgress));
    }
    
    // Add blackbox data to first arena only (with animation)
    if (i == 0 && !blackboxAircraftStates.empty()) {
      Eigen::Vector3d blackboxOffset = offset;
      std::vector<Eigen::Vector3d> a_bb = stateToVector(blackboxAircraftStates);
      
      if (a_bb.size() >= 2) {
        // Calculate blackbox progress to sync with path timing
        double blackboxProgress = 1.0; // Default to show all
        if (!blackboxAircraftStates.empty()) {
          // The blackbox test data starts at some offset but should sync with path time 0
          // Simply use the blackbox duration and sync with currentSimTime
          double blackboxStartTime = blackboxAircraftStates.front().getSimTimeMsec();
          double blackboxEndTime = blackboxAircraftStates.back().getSimTimeMsec();
          double blackboxDuration = (blackboxEndTime - blackboxStartTime) / 1000000.0; // Convert microseconds to seconds
          
          if (blackboxDuration > 0.0) {
            // Blackbox animates synchronized with path: both start at currentSimTime=0
            blackboxProgress = std::min(1.0, currentSimTime / blackboxDuration);
          }
          
        }
        
        if (inDecodeMode && !testSpans.empty() && !showingFullFlight) {
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates), blackboxProgress));
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        } else if (inDecodeMode && showingFullFlight && !testSpans.empty()) {
          // For full flight mode, still animate the dimmed background tape
          std::vector<Eigen::Vector3d> fullPoints = stateToVector(fullBlackboxAircraftStates);
          if (fullPoints.size() >= 2) {
            double fullFlightProgress = 1.0;
            if (!fullBlackboxAircraftStates.empty()) {
              double fullStartTime = fullBlackboxAircraftStates.front().getSimTimeMsec() / 1000.0;
              double fullEndTime = fullBlackboxAircraftStates.back().getSimTimeMsec() / 1000.0;
              double fullDuration = fullEndTime - fullStartTime;
              
              if (fullDuration > 0.0) {
                // Simple progress for full flight too
                fullFlightProgress = std::min(1.0, currentSimTime / fullDuration);
              }
            }
            
            this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, fullPoints, stateToOrientation(fullBlackboxAircraftStates), fullFlightProgress));
            blackboxActor->GetProperty()->SetOpacity(0.25);
            blackboxActor->GetBackfaceProperty()->SetOpacity(0.25);
            
            // Animate highlighted test spans
            for (const TestSpan& span : testSpans) {
              size_t startIdx = std::min(span.startIndex, fullBlackboxAircraftStates.size());
              size_t endIdx = std::min(span.endIndex + 1, fullBlackboxAircraftStates.size());
              
              if (startIdx < endIdx && startIdx < fullBlackboxAircraftStates.size()) {
                std::vector<AircraftState> spanStates;
                for (size_t j = startIdx; j < endIdx; j++) {
                  spanStates.push_back(fullBlackboxAircraftStates[j]);
                }
                
                if (spanStates.size() >= 2) {
                  std::vector<Eigen::Vector3d> spanPoints = stateToVector(spanStates);
                  this->blackboxHighlightTapes->AddInputData(createTapeSet(blackboxOffset, spanPoints, stateToOrientation(spanStates), fullFlightProgress));
                }
              }
            }
          }
        } else {
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates), blackboxProgress));
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        }
      }
    }
  }
  
  // Update all pipelines
  this->paths->Update();
  this->actuals->Update();
  this->segmentGaps->Update();
  if (!blackboxAircraftStates.empty()) {
    this->blackboxTapes->Update();
    this->blackboxHighlightTapes->Update();
  }
  
  // Render the scene
  renderWindow->Render();
  
  // Continue animation if still active (use one-shot timers like original)
  if (isPlaybackActive) {
    animationTimerId = renderWindowInteractor->CreateOneShotTimer(33); // 33ms = ~30fps
  }
}

void Renderer::renderFullScene() {
  // Clear existing data
  this->paths->RemoveAllInputs();
  this->actuals->RemoveAllInputs();
  this->segmentGaps->RemoveAllInputs();
  this->blackboxTapes->RemoveAllInputs();
  this->blackboxHighlightTapes->RemoveAllInputs();
  
  // Always ensure highlight tapes has at least empty data
  vtkNew<vtkPolyData> emptyHighlightData;
  vtkNew<vtkPoints> emptyHighlightPoints;
  emptyHighlightData->SetPoints(emptyHighlightPoints);
  this->blackboxHighlightTapes->AddInputData(emptyHighlightData);
  
  // Render with full progress (1.0) using existing data
  int maxArenas = (!blackboxAircraftStates.empty()) ? 1 : evalResults.pathList.size();
  for (int i = 0; i < maxArenas; i++) {
    Eigen::Vector3d offset = renderingOffset(i);
    
    std::vector<Eigen::Vector3d> p = pathToVector(evalResults.pathList[i]);
    std::vector<Eigen::Vector3d> a = stateToVector(evalResults.aircraftStateList[i]);
    
    if (!p.empty()) {
      this->paths->AddInputData(createPointSet(offset, p)); // Full progress (no timeProgress param)
    }
    if (!a.empty()) {
      this->actuals->AddInputData(createTapeSet(offset, a, stateToOrientation(evalResults.aircraftStateList[i]))); // Full progress
    }
    if (!a.empty() && !p.empty()) {
      this->segmentGaps->AddInputData(createSegmentSet(offset, evalResults.aircraftStateList[i], p)); // Full progress
    }
    
    // Add blackbox data to first arena only (full progress)
    if (i == 0 && !blackboxAircraftStates.empty()) {
      Eigen::Vector3d blackboxOffset = offset;
      std::vector<Eigen::Vector3d> a_bb = stateToVector(blackboxAircraftStates);
      
      if (a_bb.size() >= 2) {
        if (inDecodeMode && !testSpans.empty() && !showingFullFlight) {
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates)));
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        } else if (inDecodeMode && showingFullFlight && !testSpans.empty()) {
          createHighlightedFlightTapes(blackboxOffset);
        } else {
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates)));
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        }
      }
    }
  }
  
  // Update all pipelines
  this->paths->Update();
  this->actuals->Update();
  this->segmentGaps->Update();
  if (!blackboxAircraftStates.empty()) {
    this->blackboxTapes->Update();
    this->blackboxHighlightTapes->Update();
  }
  
  // Render the scene
  renderWindow->Render();
}

void printUsage(const char* progName) {
  std::cout << "Usage: " << progName << " [OPTIONS]\n";
  std::cout << "Options:\n";
  std::cout << "  -k, --keyname KEYNAME    Specify GP log key name\n";
  std::cout << "  -d, --decoder COMMAND    Specify shell command to generate CSV data\n";
  std::cout << "                           Use '-' to read CSV data from stdin\n";
  std::cout << "                           If not specified, no blackbox data is loaded\n";
  std::cout << "  -h, --help               Show this help message\n";
  std::cout << "\n";
  std::cout << "Examples:\n";
  std::cout << "  " << progName << "                                    # Render all arenas without blackbox data\n";
  std::cout << "  " << progName << " -d 'blackbox_decode file.bbl'     # Pipe through decoder command\n";
  std::cout << "  " << progName << " -d -                              # Read CSV from stdin\n";
  std::cout << "  cat data.csv | " << progName << " -d -               # Read CSV from stdin (alternative)\n";
}

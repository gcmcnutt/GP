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
#include <algorithm>

#include "renderer.h"
#include "config_manager.h"
#include "autoc.h"
#include "aircraft_state.h"

#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkCamera.h>
#include <vtkVectorText.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolygon.h>

#include <gp.h>

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/ListObjectsV2Request.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/client/ClientConfiguration.h>

using scalar = gp_scalar;
using vec3 = gp_vec3;
using quat = gp_quat;
using mat3 = Eigen::Matrix<gp_scalar, 3, 3>;

EvalResults evalResults;
std::string computedKeyName = "";
Renderer renderer;

// Blackbox-related global variables
std::string decoderCommand = "";
vtkSmartPointer<vtkAppendPolyData> blackboxTape;
std::vector<vec3> blackboxPoints;
std::vector<vec3> blackboxNormals;
std::vector<AircraftState> blackboxAircraftStates;
std::vector<AircraftState> fullBlackboxAircraftStates;  // Store full flight for span extraction
vec3 blackboxOrigin(0.0f, 0.0f, 0.0f);
scalar blackboxTimeOffset = 0.0f;
std::vector<std::string> csvLines;  // Store CSV lines for span analysis

// Xiao-related global variables
std::string xiaoLogFile = "";

// Timestamped vec structure for craft-to-target vectors
struct TimestampedVec {
  vec3 vector;
  unsigned long timestampUs;  // Timestamp in microseconds

  TimestampedVec() : vector(0.0f, 0.0f, 0.0f), timestampUs(0) {}
  TimestampedVec(vec3 vec, unsigned long ts) : vector(vec), timestampUs(ts) {}
};

// SpanData structure for organizing xiao data by test span
struct SpanData {
  vec3 origin;
  std::vector<TimestampedVec> vecs;
  size_t startStateIdx;
  size_t endStateIdx;
};

std::vector<TimestampedVec> craftToTargetVectors;  // Vec field with timestamps (all spans combined)
std::vector<vec3> xiaoVirtualPositions;  // Virtual (pos) coordinates for individual test spans
std::vector<std::string> xiaoLines;  // Store xiao log lines for span analysis
std::vector<SpanData> xiaoSpanData;  // Rabbit and vec data organized by span

// Forward declarations
bool parseBlackboxData(const std::string& csvData);
bool loadBlackboxData();
bool parseXiaoData(const std::string& xiaoLogPath);
bool loadXiaoData();
void extractXiaoTestSpans();
void updateBlackboxForCurrentTest();
void printUsage(const char* progName);

std::shared_ptr<Aws::S3::S3Client> getS3Client() {
  return ConfigManager::getS3Client();
}

// Function to lay out the squares
// given i, and NUM_PATHS_PER_GEN, compute the offsets for this particular square
vec3 Renderer::renderingOffset(int i) {
  // Calculate the dimension of the larger square
  int sideLength = static_cast<int>(std::ceil(std::sqrt(static_cast<scalar>(evalResults.pathList.size()))));

  int row = i / sideLength;
  int col = i % sideLength;

  // for now put them in a line
  scalar xOffset = static_cast<scalar>(col) * static_cast<scalar>(FIELD_SIZE + FIELD_GAP);
  scalar yOffset = static_cast<scalar>(row) * static_cast<scalar>(FIELD_SIZE + FIELD_GAP);

  return vec3(xOffset, yOffset, 0.0f);
}


vtkSmartPointer<vtkPolyData> Renderer::createPointSet(vec3 offset, const std::vector<vec3> points) {
  return createPointSet(offset, points, static_cast<scalar>(1.0f));
}

vtkSmartPointer<vtkPolyData> Renderer::createPointSet(vec3 offset, const std::vector<vec3> points, scalar timeProgress) {
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  
  // Calculate how many points to include based on time progress
  size_t numPointsToShow = static_cast<size_t>(points.size() * timeProgress);
  if (numPointsToShow == 0 && timeProgress > static_cast<scalar>(0.0f)) numPointsToShow = 1; // Show at least one point if progress > 0
  if (numPointsToShow > points.size()) numPointsToShow = points.size();
  
  // Return empty polydata if no points to show
  if (numPointsToShow == 0 || points.empty()) {
    vtkSmartPointer<vtkPoints> emptyPoints = vtkSmartPointer<vtkPoints>::New();
    polyData->SetPoints(emptyPoints);
    return polyData;
  }
  
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  for (size_t i = 0; i < numPointsToShow; ++i) {
    vec3 rPoint = points[i] + offset;
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

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(vec3 offset, const std::vector<AircraftState> state, const std::vector<vec3> end) {
  return createSegmentSet(offset, state, end, static_cast<scalar>(1.0f));
}

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(vec3 offset, const std::vector<AircraftState> state, const std::vector<vec3> end, scalar timeProgress) {
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  
  // Calculate how many segments to show based on time progress
  size_t numStatesToShow = static_cast<size_t>(state.size() * timeProgress);
  if (numStatesToShow == 0 && timeProgress > static_cast<scalar>(0.0f)) numStatesToShow = 1;
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
    vec3 rStart = vec3{ s.getPosition()[0], s.getPosition()[1], s.getPosition()[2] } + offset;
    vec3 rEnd = end[s.getThisPathIndex()] + offset;
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
vtkSmartPointer<vtkPolyData> Renderer::createTapeSet(vec3 offset, const std::vector<vec3> points,
  const std::vector<vec3> normals) {
  return createTapeSet(offset, points, normals, static_cast<scalar>(1.0f));
}

vtkSmartPointer<vtkPolyData> Renderer::createTapeSet(vec3 offset, const std::vector<vec3> points,
  const std::vector<vec3> normals, scalar timeProgress) {
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkFloatArray> orientations = vtkSmartPointer<vtkFloatArray>::New();
  orientations->SetNumberOfComponents(3);
  orientations->SetName("Orientations");

  // Calculate how many points to show based on time progress
  size_t numPointsToShow = static_cast<size_t>(points.size() * timeProgress);
  if (numPointsToShow == 0 && timeProgress > static_cast<scalar>(0.0f)) numPointsToShow = 1;
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
    vec3 point = points[i] + offset;
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
  this->xiaoVecArrows->RemoveAllInputs();

  // Always ensure these have at least empty data to prevent VTK pipeline errors
  vtkNew<vtkPolyData> emptyHighlightData;
  vtkNew<vtkPoints> emptyHighlightPoints;
  emptyHighlightData->SetPoints(emptyHighlightPoints);
  this->blackboxHighlightTapes->AddInputData(emptyHighlightData);

  vtkNew<vtkPolyData> emptyVecData;
  vtkNew<vtkPoints> emptyVecPoints;
  emptyVecData->SetPoints(emptyVecPoints);
  this->xiaoVecArrows->AddInputData(emptyVecData);


  vtkRenderer* activeRenderer = renderWindow->GetRenderers()->GetFirstRenderer();
  for (auto& label : arenaLabelActors) {
    if (activeRenderer) {
      activeRenderer->RemoveActor(label);
    }
  }
  arenaLabelActors.clear();

  std::cout << "Arena summary (" << evalResults.pathList.size() << " entries)" << std::endl;
  bool hasPerPathMetadata = evalResults.scenarioList.size() == evalResults.pathList.size();
  for (int i = 0; i < evalResults.pathList.size(); ++i) {
    ScenarioMetadata meta = evalResults.scenario;
    if (hasPerPathMetadata) {
      meta = evalResults.scenarioList[i];
    } else {
      meta.pathVariantIndex = i;
      if (meta.windSeed == 0 && meta.windVariantIndex == 0 && evalResults.scenario.windSeed != 0) {
        meta.windSeed = evalResults.scenario.windSeed;
        meta.windVariantIndex = evalResults.scenario.windVariantIndex;
      }
    }

    vec3 finalPos = vec3::Zero();
    if (i < evalResults.aircraftStateList.size() && !evalResults.aircraftStateList[i].empty()) {
      finalPos = evalResults.aircraftStateList[i].back().getPosition();
    } else if (i < evalResults.pathList.size() && !evalResults.pathList[i].empty()) {
      finalPos = evalResults.pathList[i].back().start;
    }

    std::cout << "  Arena " << i
              << " pathIdx=" << meta.pathVariantIndex
              << " windIdx=" << meta.windVariantIndex
              << " windSeed=" << meta.windSeed
              << " finalXYZ=" << finalPos.transpose() << std::endl;
  }

  int maxArenas = evalResults.pathList.size();
  for (int i = 0; i < maxArenas; i++) {
    vec3 offset = renderingOffset(i);

    std::vector<vec3> p = pathToVector(evalResults.pathList[i]);
    std::vector<vec3> a = stateToVector(evalResults.aircraftStateList[i]);

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

    scalar width = static_cast<scalar>(FIELD_SIZE);
    scalar height = static_cast<scalar>(FIELD_SIZE);
    int resolution = static_cast<int>(FIELD_SIZE / 10.0f);
    planeSource->SetOrigin(offset[0] - width / 2.0f, offset[1] - height / 2.0f, 0.0);
    planeSource->SetPoint1(offset[0] + width / 2.0f, offset[1] - height / 2.0f, 0.0);
    planeSource->SetPoint2(offset[0] - width / 2.0f, offset[1] + height / 2.0f, 0.0);
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
        unsigned char rgb[4] = { 255, 255, 255, 100 };
        cellData->InsertTypedTuple(i, rgb);
      }
      else {
        unsigned char rgb[4] = { 0, 0, 0, 100 };
        cellData->InsertTypedTuple(i, rgb);
      }
    }
    planeSource->GetOutput()->GetCellData()->SetScalars(cellData);
    planeSource->Update();
    planeData->AddInputConnection(planeSource->GetOutputPort());

    ScenarioMetadata meta = evalResults.scenario;
    if (evalResults.scenarioList.size() == evalResults.pathList.size()) {
      meta = evalResults.scenarioList[i];
    } else {
      meta.pathVariantIndex = i;
    }

    std::ostringstream labelStream;
    labelStream << "Arena " << i;
    if (meta.pathVariantIndex >= 0) {
      labelStream << "  Path " << meta.pathVariantIndex;
    }
    labelStream << "  Wind " << meta.windVariantIndex;
    labelStream << "\nSeed " << meta.windSeed;

    scalar labelScale = static_cast<scalar>(6.0f * 0.5f);
    scalar textOffsetX = static_cast<scalar>(45.0f);
    scalar textOffsetY = static_cast<scalar>(-45.0f);
    scalar anchorX = offset[0] + textOffsetX;
    scalar anchorY = offset[1] + textOffsetY;
    scalar anchorZ = offset[2] + static_cast<scalar>(0.1f);

    vtkNew<vtkVectorText> textSource;
    textSource->SetText(labelStream.str().c_str());
    textSource->Update();
    const auto* bounds = textSource->GetOutput()->GetBounds();

    vtkNew<vtkTransform> textTransform;
    textTransform->PostMultiply();
    textTransform->Translate(-bounds[0], -bounds[3], -bounds[4]);
    textTransform->Scale(labelScale, -labelScale, labelScale);
    textTransform->RotateZ(90.0);
    textTransform->Translate(anchorX, anchorY, anchorZ);

    vtkNew<vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetTransform(textTransform);
    transformFilter->SetInputConnection(textSource->GetOutputPort());

    vtkNew<vtkPolyDataMapper> textMapper;
    textMapper->SetInputConnection(transformFilter->GetOutputPort());

    vtkSmartPointer<vtkActor> labelActor = vtkSmartPointer<vtkActor>::New();
    labelActor->SetMapper(textMapper);
    labelActor->GetProperty()->SetColor(0.9, 0.9, 0.9);
    labelActor->GetProperty()->SetOpacity(0.5);
    labelActor->GetProperty()->SetLighting(false);

    if (activeRenderer) {
      activeRenderer->AddActor(labelActor);
    }
    arenaLabelActors.push_back(labelActor);

    // Add blackbox data to first arena only
    if (i == 0 && !blackboxAircraftStates.empty()) {
      // Center blackbox data in the arena by adding the arena offset
      // The blackbox data is already origin-centered, so we add the arena offset to position it properly
      vec3 blackboxOffset = offset;
      
      // Use the same rendering pipeline as the blue/yellow tape
      std::vector<vec3> a = stateToVector(blackboxAircraftStates);
      
      // Only create tape if we have enough points
      if (a.size() >= 2) {
        if ((inDecodeMode || inXiaoMode) && !testSpans.empty() && !showingFullFlight) {
          // Create regular tape for current test span only
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a, stateToOrientation(blackboxAircraftStates)));
          // Reset opacity to full brightness for single span display
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        } else if ((inDecodeMode || inXiaoMode) && showingFullFlight && !testSpans.empty()) {
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

      // Render vec arrows (craft-to-target vectors) for xiao mode
      if (inXiaoMode && !craftToTargetVectors.empty() && !blackboxAircraftStates.empty()) {
        size_t arrowCount = 0;
        size_t skippedCount = 0;

        for (const auto& timestampedVec : craftToTargetVectors) {
          vec3 vecDir = timestampedVec.vector;
          if (vecDir.norm() > 0.01f) {
            // Find closest aircraft state by timestamp to get position
            scalar vecTime = static_cast<scalar>(timestampedVec.timestampUs);
            scalar closestTimeDiff = std::numeric_limits<scalar>::max();
            vec3 craftPos = blackboxAircraftStates[0].getPosition();

            for (const auto& state : blackboxAircraftStates) {
              scalar stateTime = static_cast<scalar>(state.getSimTimeMsec());
              scalar timeDiff = std::abs(stateTime - vecTime);
              if (timeDiff < closestTimeDiff) {
                closestTimeDiff = timeDiff;
                craftPos = state.getPosition();
              }
            }

            craftPos = craftPos + blackboxOffset;
            scalar vecLength = vecDir.norm();

            vtkNew<vtkArrowSource> arrowSource;
            arrowSource->SetTipResolution(8);
            arrowSource->SetShaftResolution(8);
            arrowSource->SetShaftRadius(0.1);    // Match blue line thickness (0.1 radius)
            arrowSource->SetTipRadius(0.2);      // Proportional tip thickness
            arrowSource->SetTipLength(0.12);     // Stubbier tip
            arrowSource->Update();

            vtkNew<vtkTransform> transform;
            transform->Translate(craftPos[0], craftPos[1], craftPos[2]);

            vec3 vecNorm = vecDir / vecLength;
            vec3 defaultDir(1.0, 0.0, 0.0);
            vec3 rotAxis = defaultDir.cross(vecNorm);
            scalar rotAxisLen = rotAxis.norm();
            if (rotAxisLen > 0.001f) {
              rotAxis = rotAxis / rotAxisLen;
              scalar dotProd = std::max(-1.0f, std::min(1.0f, defaultDir.dot(vecNorm)));
              scalar angle = std::acos(dotProd) * 180.0f / M_PI;
              transform->RotateWXYZ(angle, rotAxis[0], rotAxis[1], rotAxis[2]);
            }
            // Scale only in X direction (length) to keep shaft/tip thickness absolute
            transform->Scale(vecLength, 1.0, 1.0);

            vtkNew<vtkTransformPolyDataFilter> transformFilter;
            transformFilter->SetInputConnection(arrowSource->GetOutputPort());
            transformFilter->SetTransform(transform);
            transformFilter->Update();

            this->xiaoVecArrows->AddInputData(transformFilter->GetOutput());
            arrowCount++;
          } else {
            skippedCount++;
          }
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

  // Update xiao-specific actors
  if (inXiaoMode) {
    if (this->xiaoVecArrows->GetNumberOfInputConnections(0) > 0) {
      this->xiaoVecArrows->Update();
    }
  }

  // Update the window title
  std::string title = keyName + " - " + std::to_string(10000 - newGen);
  renderWindow->SetWindowName(title.c_str());

  // Extract fitness from GP data and update text displays
  gp_scalar fitness = extractFitnessFromGP(evalResults.gp);
  updateTextDisplay(newGen, fitness);

  // Render the updated scene
  renderWindow->Render();
  focusMode = false;
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
      renderer_->hideStopwatch();
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
      renderer_->hideStopwatch();
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
    // Stop any active playback animation before exiting to prevent crashes
    if (this->RendererPtr && this->RendererPtr->isPlaybackActive) {
      this->RendererPtr->togglePlaybackAnimation(); // This will stop the animation
    }
    if (this->Interactor)
      this->Interactor->TerminateApp();
  }
  vtkRenderWindowInteractor* Interactor;
  ::Renderer* RendererPtr;
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
      // Update stopwatch position when window is resized
      renderer_->updateStopwatchPosition();
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
  xiaoVecArrows = vtkSmartPointer<vtkAppendPolyData>::New();

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
  xiaoVecArrows->AddInputData(emptyPolyData);

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
  blackboxProperty->SetColor(1.0, 0.5, 0.0);  // Orange for top face (canopy)
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
  blackboxBackProperty->SetColor(0.0, 1.0, 0.0);  // Green for bottom face (belly)
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
  highlightProperty->SetColor(1.0, 0.6, 0.0);  // Bright orange for top face (canopy)
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
  highlightBackProperty->SetColor(0.0, 1.0, 0.0);  // Bright green for bottom face (belly)
  highlightBackProperty->SetAmbient(0.2);
  highlightBackProperty->SetDiffuse(1.0);
  highlightBackProperty->SetSpecular(0.2);
  highlightBackProperty->SetSpecularPower(20);
  highlightBackProperty->SetOpacity(1.0);
  
  // Set the back face property for highlights
  blackboxHighlightActor->SetBackfaceProperty(highlightBackProperty);

  // Create xiao vec arrow actor for visualizing craft-to-target vectors
  xiaoVecActor = vtkSmartPointer<vtkActor>::New();
  vtkNew<vtkPolyDataMapper> xiaoVecMapper;
  xiaoVecMapper->SetInputConnection(xiaoVecArrows->GetOutputPort());
  xiaoVecActor->SetMapper(xiaoVecMapper);
  xiaoVecActor->GetProperty()->SetColor(0.5, 0.0, 0.8);  // Purple arrows
  xiaoVecActor->GetProperty()->SetLineWidth(1.0);

  // Update actors
  actor1->SetMapper(mapper1);
  actor2->SetMapper(mapper2);
  actor3->SetMapper(mapper3);

  // Set properties for the tape (actor2)
  vtkProperty* property = actor2->GetProperty();

  // Enable two-sided surface rendering
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
  renderer->AddActor(actor1);  // Red path line
  renderer->AddActor(actor2);  // Yellow flight tape
  renderer->AddActor(actor3);  // Blue delta lines
  
  // Only add blackbox actors if there's blackbox data
  if (!blackboxAircraftStates.empty()) {
    renderer->AddActor(blackboxActor);
    renderer->AddActor(blackboxHighlightActor);
    // Add xiao-specific actors if in xiao mode
    if (inXiaoMode) {
      renderer->AddActor(xiaoVecActor);
    }
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

  // Create stopwatch (initially hidden)
  createStopwatch();
  // Create control HUD (initially hidden)
  createControlsOverlay();

  // render
  renderWindow->Render();

  // exit command
  vtkNew<ExitCommand> exitCommand;
  exitCommand->Interactor = renderWindowInteractor;
  exitCommand->RendererPtr = this;
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

std::vector<vec3> Renderer::pathToVector(std::vector<Path> path) {
  std::vector<vec3> points;
  for (const auto& p : path) {
    points.push_back(p.start);
  }
  return points;
}

std::vector<vec3> Renderer::stateToVector(std::vector<AircraftState> state) {
  std::vector<vec3> points;
  for (const auto& s : state) {
    points.push_back(s.getPosition());
  }
  return points;
}

std::vector<vec3> Renderer::stateToOrientation(std::vector<AircraftState> state) {
  std::vector<vec3> points;
  for (const auto& s : state) {
    // Use the aircraft body Z-axis (up direction) for ribbon orientation
    points.push_back(s.getOrientation() * -vec3::UnitZ());
  }
  return points;
}

gp_scalar Renderer::extractFitnessFromGP(const std::vector<char>& gpData) {
  if (gpData.empty()) {
    return 0.0f;
  }
  
  try {
    // Create stream from the char vector
    boost::iostreams::stream<boost::iostreams::array_source> inStream(gpData.data(), gpData.size());
    
    // Create and load a base GP object
    GP gp;
    gp.load(inStream);
    
    return static_cast<gp_scalar>(gp.getFitness());
  }
  catch (const std::exception& e) {
    std::cerr << "Error extracting fitness from GP: " << e.what() << std::endl;
    return 0.0f;
  }
}

void Renderer::updateTextDisplay(int generation, gp_scalar fitness) {
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
  
  // Update test display based on decode mode or xiao mode and current state
  if ((inDecodeMode || inXiaoMode) && !testSpans.empty()) {
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
    // Hide test actors when not in decode/xiao mode or no test spans
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
    {"xiaofile", required_argument, 0, 'x'},
    {"config", required_argument, 0, 'i'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  std::string configFile = "autoc.ini";
  int option_index = 0;
  int c;

  while ((c = getopt_long(argc, argv, "k:d:x:i:h", long_options, &option_index)) != -1) {
    switch (c) {
      case 'k':
        computedKeyName = optarg;
        break;
      case 'd':
        decoderCommand = optarg;
        break;
      case 'x':
        xiaoLogFile = optarg;
        break;
      case 'i':
        configFile = optarg;
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

  // Check mutual exclusivity
  if (!decoderCommand.empty() && !xiaoLogFile.empty()) {
    std::cerr << "Error: -d and -x options are mutually exclusive" << std::endl;
    printUsage(argv[0]);
    return 1;
  }

  // Handle positional arguments for backward compatibility
  if (computedKeyName.empty() && optind < argc) {
    computedKeyName = argv[optind];
  }
  
  std::string keyName = "";
  
  // Initialize configuration
  ConfigManager::initialize(configFile);

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

  // Load xiao log data if specified
  if (!xiaoLogFile.empty()) {
    if (!loadXiaoData()) {
      std::cerr << "Failed to load xiao log data" << std::endl;
      return 1;
    }
    if (!blackboxAircraftStates.empty()) {
      std::cout << "Xiao log data loaded: " << blackboxAircraftStates.size() << " states" << std::endl;
      renderer.inXiaoMode = true;
      extractXiaoTestSpans();

      if (!renderer.testSpans.empty()) {
        std::cout << "Found " << renderer.testSpans.size() << " test spans" << std::endl;
        renderer.currentTestIndex = 0;
        renderer.showingFullFlight = false;
        // Filter xiao data to show first test span immediately
        updateBlackboxForCurrentTest();
        std::cout << "Showing test 1: " << blackboxAircraftStates.size() << " states" << std::endl;
      } else {
        std::cout << "No autoc=Y test spans found, showing full flight" << std::endl;
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
  std::cout << "  f - Focus camera on current arena" << std::endl;
  std::cout << "  Arrow keys - Move focus between arenas" << std::endl;
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
  int rcRollIndex = -1, rcPitchIndex = -1, rcThrottleIndex = -1;
  
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
        else if (header == "rcCommand[0]") rcRollIndex = index;
        else if (header == "rcCommand[1]") rcPitchIndex = index;
        else if (header == "rcCommand[3]") rcThrottleIndex = index;
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
      scalar navX = std::stof(row[latIndex]);  // navPos[0] in cm
      scalar navY = std::stof(row[lonIndex]);  // navPos[1] in cm
      scalar navZ = std::stof(row[altIndex]);  // navPos[2] in cm
      
      // Convert from centimeters to meters and center on origin
      static bool firstPoint = true;
      static scalar originX, originY, originZ;
      
      if (firstPoint) {
        originX = navX;
        originY = navY;
        originZ = navZ;
        firstPoint = false;
      }
      
      // Convert to meters and translate to origin
      // NED coordinates: North=X, East=Y, Down=Z
      // For VTK visualization, flip Z to make "up" positive (NED Down -> VTK Up)
      scalar x = (navX - originX) / static_cast<scalar>(100.0f);   // cm to m, North
      scalar y = (navY - originY) / static_cast<scalar>(100.0f);   // cm to m, East  
      scalar z = -(navZ - originZ) / static_cast<scalar>(100.0f);  // cm to m, Down->Up (flip sign)
      
      vec3 newPoint(x, y, z);
      
      // Filter out coincident points to avoid VTK ribbon filter issues
      if (blackboxPoints.empty() || (newPoint - blackboxPoints.back()).norm() > static_cast<scalar>(0.01f)) {
        
        blackboxPoints.push_back(newPoint);
        
        // Create AircraftState object
        if (quatWIndex >= 0 && quatXIndex >= 0 && quatYIndex >= 0 && quatZIndex >= 0) {
          // Parse normalized quaternion values (stored as integers * 10000)
          scalar qw = std::stof(row[quatWIndex]) / static_cast<scalar>(10000.0f);
          scalar qx = std::stof(row[quatXIndex]) / static_cast<scalar>(10000.0f);
          scalar qy = std::stof(row[quatYIndex]) / static_cast<scalar>(10000.0f);
          scalar qz = std::stof(row[quatZIndex]) / static_cast<scalar>(10000.0f);
          
          // INAV blackbox logs raw body->earth quaternion in NED frame
          // (same as what MSP sends before the conjugate fix in msplink.cpp)
          quat inavQuat(qw, qx, qy, qz);
          inavQuat.normalize();

          // Convert body->earth to earth->body to match renderer/GP contract
          // Apply conjugate: (w, x, y, z) -> (w, -x, -y, -z)
          quat earthToBody(inavQuat.w(), -inavQuat.x(), -inavQuat.y(), -inavQuat.z());
          earthToBody.normalize();
          
          // Get time if available, otherwise use index
          unsigned long int timeMs = (timeIndex >= 0) ? std::stoul(row[timeIndex]) : blackboxAircraftStates.size() * 100;

          // Create simple velocity vector for blackbox data (assuming forward flight)
          // earth->body quaternion rotates body vectors to earth frame
          vec3 velocity_vector = earthToBody * vec3(static_cast<scalar>(20.0f), 0, 0);

          gp_scalar pitchCmd = 0.0f;
          gp_scalar rollCmd = 0.0f;
          gp_scalar throttleCmd = 0.0f;
          if (rcRollIndex >= 0 && rcPitchIndex >= 0 && rcThrottleIndex >= 0) {
            rollCmd = CLAMP_DEF(static_cast<gp_scalar>(std::stof(row[rcRollIndex]) / static_cast<scalar>(500.0f)), -1.0f, 1.0f);
            pitchCmd = CLAMP_DEF(static_cast<gp_scalar>(std::stof(row[rcPitchIndex]) / static_cast<scalar>(500.0f)), -1.0f, 1.0f);
            // rcCommand[3] is centered around 1500us; scale to -1..1
            throttleCmd = CLAMP_DEF(static_cast<gp_scalar>((std::stof(row[rcThrottleIndex]) - static_cast<scalar>(1500.0f)) / static_cast<scalar>(500.0f)), -1.0f, 1.0f);
          }

          // Create AircraftState with earth->body quaternion (matches CRRCSim/GP contract)
          AircraftState state(static_cast<int>(blackboxAircraftStates.size()), static_cast<scalar>(20.0f), velocity_vector, earthToBody, newPoint, pitchCmd, rollCmd, throttleCmd, timeMs);
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

  if ((!renderer.inDecodeMode && !renderer.inXiaoMode) || fullBlackboxAircraftStates.empty()) {
    return;
  }
  
  // Always update the blackbox data, even if renderer isn't fully initialized yet
  // The safety check for renderWindow is only needed for text display updates
  
  if (renderer.showingFullFlight || renderer.testSpans.empty()) {
    // Show full flight - use all data
    blackboxAircraftStates = fullBlackboxAircraftStates;
    // Keep all vec data
    // (craftToTargetVectors already has all data)
  } else if (renderer.currentTestIndex >= 0 && renderer.currentTestIndex < renderer.testSpans.size()) {
    // Show only the selected test span
    const TestSpan& span = renderer.testSpans[renderer.currentTestIndex];

    // Update vec data for current span
    craftToTargetVectors = span.vecPoints;
    blackboxAircraftStates.clear();

    size_t startIdx = std::min(span.startIndex, fullBlackboxAircraftStates.size());
    size_t endIdx = std::min(span.endIndex + 1, fullBlackboxAircraftStates.size());

    // Ensure we have valid indices
    if (startIdx < endIdx && startIdx < fullBlackboxAircraftStates.size()) {
      // For xiao mode, use virtual positions (pos field) directly without additional offsetting
      // For decode mode, calculate offset to align test start with path origin
      bool useVirtualPositions = renderer.inXiaoMode && !xiaoVirtualPositions.empty() &&
                                  xiaoVirtualPositions.size() == fullBlackboxAircraftStates.size();

      vec3 positionOffset(0.0f, 0.0f, 0.0f);
      if (!useVirtualPositions) {
        // Decode mode: calculate offset to align test start with path origin
        vec3 testStartPosition = fullBlackboxAircraftStates[startIdx].getPosition();
        vec3 pathOrigin(0.0f, 0.0f, 0.0f); // Path origin is at (0,0,0) for horizontal position

        positionOffset[0] = pathOrigin[0] - testStartPosition[0]; // North offset
        positionOffset[1] = pathOrigin[1] - testStartPosition[1]; // East offset
        positionOffset[2] = SIM_INITIAL_ALTITUDE - testStartPosition[2]; // Offset Z to SIM_INITIAL_ALTITUDE
      }

      for (size_t i = startIdx; i < endIdx; i++) {
        AircraftState offsetState = fullBlackboxAircraftStates[i];

        if (useVirtualPositions) {
          // Xiao mode: use virtual position directly (already offset in xiao log)
          offsetState.setPosition(xiaoVirtualPositions[i]);
        } else {
          // Decode mode: apply calculated offset
          vec3 originalPosition = offsetState.getPosition();
          vec3 offsetPosition = originalPosition + positionOffset;
          offsetState.setPosition(offsetPosition);
        }

        // Keep original attitude (yaw, pitch, roll) - no attitude offset needed
        blackboxAircraftStates.push_back(offsetState);
      }
      
      std::cout << "Showing test " << (renderer.currentTestIndex + 1) << ": " 
                << blackboxAircraftStates.size() << " states (indices " 
                << startIdx << "-" << (endIdx-1) << ") with position offset ("
                << positionOffset[0] << ", " << positionOffset[1] << ", " << positionOffset[2] << ")" << std::endl;
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

// Parse xiao log file
bool parseXiaoData(const std::string& xiaoLogPath) {
  std::ifstream file(xiaoLogPath);
  if (!file.is_open()) {
    std::cerr << "Failed to open xiao log file: " << xiaoLogPath << std::endl;
    return false;
  }

  blackboxPoints.clear();
  blackboxNormals.clear();
  blackboxAircraftStates.clear();
  fullBlackboxAircraftStates.clear();
  xiaoLines.clear();
  craftToTargetVectors.clear();
  xiaoVirtualPositions.clear();

  // Regex patterns for xiao log parsing
  // Format: #<seqnum> <xiao_ms> <inav_ms> <level> GP ...
  std::regex timestampRe(R"(^#\d+\s+(\d+)\s+(\d+)\s+\w)");  // Capture xiao_ms (2nd) and inav_ms (3rd)
  std::regex controlEnableRe(R"(GP Control: Switch enabled - test origin NED=\[([-0-9\.]+),\s*([-0-9\.]+),\s*([-0-9\.]+)\])");
  std::regex controlDisableRe(R"(GP Control: Switch disabled)");
  std::regex inputRe(R"(GP Input:.*idx=(\d+).*rabbit=\[([-0-9\.]+),([-0-9\.]+),([-0-9\.]+)\].*vec=\[([-0-9\.]+),([-0-9\.]+),([-0-9\.]+)\].*relvel=([-0-9\.]+))");
  std::regex stateRe(R"(GP State:.*pos_raw=\[([-0-9\.]+),([-0-9\.]+),([-0-9\.]+)\].*pos=\[([-0-9\.]+),([-0-9\.]+),([-0-9\.]+)\].*vel=\[([-0-9\.]+),([-0-9\.]+),([-0-9\.]+)\].*quat=\[([-0-9\.]+),([-0-9\.]+),([-0-9\.]+),([-0-9\.]+)\])");
  std::regex autocFlagRe(R"(autoc=(Y|N))");
  std::regex outputRe(R"(GP Output: rc=\[(\d+),(\d+),(\d+)\])");

  std::string line;
  bool inSpan = false;
  unsigned long spanStartTimeMs = 0;
  vec3 currentOrigin(0.0f, 0.0f, 0.0f);
  std::vector<TimestampedVec> currentVecPoints;
  std::smatch matches;

  // Clear and prepare to collect span data
  xiaoSpanData.clear();
  SpanData currentSpanData;

  // Map to store vec vectors per state index for later rendering
  std::map<size_t, vec3> stateVecMap;

  size_t currentStateIdx = 0;
  std::map<unsigned long, vec3> timestampVecMap;  // Map timestamp to vec for later association
  std::map<unsigned long, scalar> timestampRelVelMap;  // Map timestamp to relvel

  // Store RC commands by timestamp
  struct RCCommand {
    scalar roll;    // microseconds 1000-2000
    scalar pitch;   // microseconds 1000-2000
    scalar throttle; // microseconds 1000-2000
  };
  std::map<unsigned long, RCCommand> timestampRCMap;

  // FIRST PASS: Collect all RC commands
  while (std::getline(file, line)) {
    xiaoLines.push_back(line);

    // Extract INAV timestamp
    unsigned long inavMs = 0;
    std::smatch tsMatch;
    if (std::regex_search(line, tsMatch, timestampRe)) {
      inavMs = std::stoul(tsMatch[2].str());
    }

    // Parse GP Output (RC commands)
    if (std::regex_search(line, matches, outputRe)) {
      scalar rc_roll = std::stof(matches[1].str());
      scalar rc_pitch = std::stof(matches[2].str());
      scalar rc_throttle = std::stof(matches[3].str());

      RCCommand rc;
      rc.roll = rc_roll;
      rc.pitch = rc_pitch;
      rc.throttle = rc_throttle;
      timestampRCMap[inavMs] = rc;
    }
  }

  // SECOND PASS: Process all state data with RC commands available
  for (const std::string& line : xiaoLines) {

    // Extract timestamps from line
    unsigned long xiaoMs = 0;
    unsigned long inavMs = 0;
    std::smatch tsMatch;
    if (std::regex_search(line, tsMatch, timestampRe)) {
      xiaoMs = std::stoul(tsMatch[1].str());
      inavMs = std::stoul(tsMatch[2].str());
    }

    // Check for control enable (span start)
    if (std::regex_search(line, matches, controlEnableRe)) {
      // Save previous span before starting new one
      if (inSpan && !currentSpanData.vecs.empty()) {
        currentSpanData.endStateIdx = currentStateIdx > 0 ? currentStateIdx - 1 : 0;
        xiaoSpanData.push_back(currentSpanData);
      }

      scalar origin_n = std::stof(matches[1].str());
      scalar origin_e = std::stof(matches[2].str());
      scalar origin_d = std::stof(matches[3].str());
      currentOrigin = vec3(origin_n, origin_e, origin_d);
      inSpan = true;
      spanStartTimeMs = xiaoMs;  // Record span start time (still use xiao for relative time)
      // Start new span tracking
      currentSpanData.origin = currentOrigin;
      currentSpanData.vecs.clear();
      currentSpanData.startStateIdx = currentStateIdx;
      continue;
    }

    // Check for control disable (span end)
    if (std::regex_search(line, matches, controlDisableRe)) {
      inSpan = false;
      // Save span data
      currentSpanData.endStateIdx = currentStateIdx - 1;
      if (!currentSpanData.vecs.empty()) {
        xiaoSpanData.push_back(currentSpanData);
      }
      continue;
    }

    if (!inSpan) continue;

    // Parse GP Input (vec and relvel)
    if (std::regex_search(line, matches, inputRe)) {
      int idx = std::stoi(matches[1].str());
      scalar vec_n = std::stof(matches[5].str());
      scalar vec_e = std::stof(matches[6].str());
      scalar vec_d = std::stof(matches[7].str());
      scalar relVel = std::stof(matches[8].str());

      vec3 vecDir(vec_n, vec_e, vec_d);

      // Calculate relative timestamp in microseconds (same as aircraft states)
      unsigned long relativeTimeUs = (xiaoMs - spanStartTimeMs) * 1000;

      TimestampedVec timestampedVec(vecDir, relativeTimeUs);
      currentVecPoints.push_back(timestampedVec);
      currentSpanData.vecs.push_back(timestampedVec);  // Also add to span-specific data

      // Store vec and relvel with inavMs timestamp for matching with GP State
      timestampVecMap[inavMs] = vecDir;
      timestampRelVelMap[inavMs] = relVel;
    }

    // Parse GP State (position and attitude)
    if (std::regex_search(line, matches, stateRe)) {
      // Check autoc flag
      std::smatch autocMatch;
      std::string autocFlag = "N";
      if (std::regex_search(line, autocMatch, autocFlagRe)) {
        autocFlag = autocMatch[1].str();
      }

      // Extract positions
      scalar pos_raw_n = std::stof(matches[1].str());
      scalar pos_raw_e = std::stof(matches[2].str());
      scalar pos_raw_d = std::stof(matches[3].str());
      scalar pos_n = std::stof(matches[4].str());
      scalar pos_e = std::stof(matches[5].str());
      scalar pos_d = std::stof(matches[6].str());

      // Extract velocities
      scalar vel_n = std::stof(matches[7].str());
      scalar vel_e = std::stof(matches[8].str());
      scalar vel_d = std::stof(matches[9].str());

      // Extract quaternion (earth->body in xiao logs)
      scalar qw = std::stof(matches[10].str());
      scalar qx = std::stof(matches[11].str());
      scalar qy = std::stof(matches[12].str());
      scalar qz = std::stof(matches[13].str());

      quat earthToBody(qw, qx, qy, qz);
      earthToBody.normalize();

      // Use pos_raw (raw NED coordinates) for 'a' all flight display
      // Store pos (virtual coordinates) separately for individual test span display
      // In xiao logs, D is already in meters with negative = up, same as VTK convention
      vec3 position(pos_raw_n, pos_raw_e, pos_raw_d);  // Raw position (for 'a' all mode)
      vec3 virtualPosition(pos_n, pos_e, pos_d);        // Virtual position (for individual span mode)

      // Filter out coincident points
      if (blackboxPoints.empty() || (position - blackboxPoints.back()).norm() > static_cast<scalar>(0.01f)) {
        blackboxPoints.push_back(position);
        xiaoVirtualPositions.push_back(virtualPosition);

        // Calculate velocity vector
        vec3 velocity_vector(vel_n, vel_e, vel_d);

        // Use relVel from GP Input if available, otherwise compute from velocity vector
        scalar speed = velocity_vector.norm();
        if (timestampRelVelMap.find(inavMs) != timestampRelVelMap.end()) {
          speed = timestampRelVelMap[inavMs];
        }

        // Calculate relative time in microseconds (to match blackbox convention)
        unsigned long relativeTimeUs = (xiaoMs - spanStartTimeMs) * 1000;  // Convert ms to us

        // Get RC commands from INAV timestamp if available
        scalar pitchCmd = 0.0f;
        scalar rollCmd = 0.0f;
        scalar throttleCmd = 0.0f;

        if (timestampRCMap.find(inavMs) != timestampRCMap.end()) {
          const RCCommand& rc = timestampRCMap[inavMs];
          // Convert from microseconds (1000-2000) to normalized (-1 to +1)
          // Center is 1500us, range is ±500us
          rollCmd = (rc.roll - 1500.0f) / 500.0f;
          pitchCmd = (rc.pitch - 1500.0f) / 500.0f;
          throttleCmd = (rc.throttle - 1500.0f) / 500.0f;

          // Clamp to valid range
          rollCmd = std::max(-1.0f, std::min(1.0f, rollCmd));
          pitchCmd = std::max(-1.0f, std::min(1.0f, pitchCmd));
          throttleCmd = std::max(-1.0f, std::min(1.0f, throttleCmd));
        }

        // Create AircraftState with real timestamp and RC commands
        AircraftState state(
          static_cast<int>(blackboxAircraftStates.size()),
          speed,
          velocity_vector,
          earthToBody,
          position,
          pitchCmd,
          rollCmd,
          throttleCmd,
          relativeTimeUs
        );

        // Store vec vector for this state if available (match by INAV timestamp)
        if (timestampVecMap.find(inavMs) != timestampVecMap.end()) {
          stateVecMap[blackboxAircraftStates.size()] = timestampVecMap[inavMs];
        }

        blackboxAircraftStates.push_back(state);
        fullBlackboxAircraftStates.push_back(state);
      }
    }
  }

  file.close();

  // Save any unclosed span at end of file
  if (inSpan && !currentSpanData.vecs.empty()) {
    currentSpanData.endStateIdx = currentStateIdx > 0 ? currentStateIdx - 1 : 0;
    xiaoSpanData.push_back(currentSpanData);
  }

  // Store vec vectors globally for rendering
  craftToTargetVectors = currentVecPoints;

  // Count non-zero vec entries for debugging
  size_t nonZeroVecs = 0;
  for (const auto& v : craftToTargetVectors) {
    if (v.vector.norm() > 0.01f) nonZeroVecs++;
  }

  std::cout << "Parsed xiao log: " << blackboxAircraftStates.size() << " states, "
            << craftToTargetVectors.size() << " total vecs, "
            << nonZeroVecs << " non-zero vecs" << std::endl;

  // Show first few vecs with timestamps to debug duplicates
  if (!craftToTargetVectors.empty()) {
    std::cout << "\n=== VEC SAMPLE (first 10) ===" << std::endl;
    size_t showVecs = std::min(size_t(10), craftToTargetVectors.size());
    for (size_t i = 0; i < showVecs; i++) {
      const auto& v = craftToTargetVectors[i];
      std::cout << "  vec[" << i << "] = ["
                << v.vector[0] << ", "
                << v.vector[1] << ", "
                << v.vector[2] << "] @ " << v.timestampUs << "us, norm=" << v.vector.norm() << std::endl;
    }
    std::cout << "===================\n" << std::endl;
  }

  return !blackboxAircraftStates.empty();
}

// Load xiao log data
bool loadXiaoData() {
  return parseXiaoData(xiaoLogFile);
}

// Extract test spans from xiao log (spans where autoc=Y)
void extractXiaoTestSpans() {
  extern Renderer renderer;
  renderer.testSpans.clear();

  if (xiaoLines.empty() || fullBlackboxAircraftStates.empty()) {
    std::cerr << "No xiao data or aircraft states available for test span extraction" << std::endl;
    return;
  }

  // Regex patterns
  std::regex controlEnableRe(R"(GP Control: Switch enabled - test origin NED=\[([-0-9\.]+),\s*([-0-9\.]+),\s*([-0-9\.]+)\])");
  std::regex stateRe(R"(GP State:)");
  std::regex autocFlagRe(R"(autoc=(Y|N))");

  bool inControlSpan = false;
  bool inAutocSpan = false;
  TestSpan currentSpan;
  vec3 currentOrigin(0.0f, 0.0f, 0.0f);
  size_t stateIndex = 0;
  std::smatch matches;

  for (const std::string& line : xiaoLines) {
    // Check for control enable
    if (std::regex_search(line, matches, controlEnableRe)) {
      scalar origin_n = std::stof(matches[1].str());
      scalar origin_e = std::stof(matches[2].str());
      scalar origin_d = std::stof(matches[3].str());
      currentOrigin = vec3(origin_n, origin_e, origin_d);
      inControlSpan = true;
      continue;
    }

    // Skip if not in control span
    if (!inControlSpan) continue;

    // Check for GP State line
    if (std::regex_search(line, stateRe)) {
      if (stateIndex >= fullBlackboxAircraftStates.size()) break;

      // Check autoc flag
      std::smatch autocMatch;
      std::string autocFlag = "N";
      if (std::regex_search(line, autocMatch, autocFlagRe)) {
        autocFlag = autocMatch[1].str();
      }

      if (autocFlag == "Y" && !inAutocSpan) {
        // Start of autoc=Y span
        currentSpan.startIndex = stateIndex;
        currentSpan.startTime = stateIndex * 100;  // Fake timestamp
        currentSpan.origin = currentOrigin;
        inAutocSpan = true;
      } else if (autocFlag == "N" && inAutocSpan) {
        // End of autoc=Y span
        currentSpan.endIndex = stateIndex - 1;
        currentSpan.endTime = (stateIndex - 1) * 100;
        renderer.testSpans.push_back(currentSpan);
        inAutocSpan = false;
      }

      stateIndex++;
    }
  }

  // Handle span that continues to end of log
  if (inAutocSpan && stateIndex > 0) {
    currentSpan.endIndex = stateIndex - 1;
    currentSpan.endTime = (stateIndex - 1) * 100;
    renderer.testSpans.push_back(currentSpan);
  }

  // Validate spans
  std::vector<TestSpan> validSpans;
  for (const TestSpan& span : renderer.testSpans) {
    if (span.startIndex < fullBlackboxAircraftStates.size() &&
        span.endIndex < fullBlackboxAircraftStates.size() &&
        span.startIndex <= span.endIndex) {
      validSpans.push_back(span);
    } else {
      std::cerr << "Invalid span removed: indices " << span.startIndex << "-" << span.endIndex
                << " (max: " << fullBlackboxAircraftStates.size() - 1 << ")" << std::endl;
    }
  }
  renderer.testSpans = validSpans;

  // Copy vec data from xiaoSpanData to TestSpans
  for (size_t i = 0; i < renderer.testSpans.size() && i < xiaoSpanData.size(); i++) {
    renderer.testSpans[i].vecPoints = xiaoSpanData[i].vecs;
  }

  std::cout << "Extracted " << renderer.testSpans.size() << " valid test spans from " << stateIndex << " aircraft states" << std::endl;
  for (size_t i = 0; i < renderer.testSpans.size(); i++) {
    std::cout << "  Test " << (i+1) << ": indices " << renderer.testSpans[i].startIndex << "-" << renderer.testSpans[i].endIndex
              << " (" << (renderer.testSpans[i].endIndex - renderer.testSpans[i].startIndex + 1) << " states)"
              << " origin=[" << renderer.testSpans[i].origin[0] << ", "
              << renderer.testSpans[i].origin[1] << ", "
              << renderer.testSpans[i].origin[2] << "]"
              << " vecs=" << renderer.testSpans[i].vecPoints.size() << std::endl;
  }
}

void Renderer::jumpToNewestGeneration() {
  hideStopwatch();
  
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
  hideStopwatch();
  
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
  
  // Hide stopwatch when changing test data
  hideStopwatch();
  
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
  
  // Hide stopwatch when changing test data
  hideStopwatch();
  
  // Update the blackbox rendering to show only the selected test span
  updateBlackboxForCurrentTest();
  
  // Re-render using in-memory data (no S3 fetch)
  renderFullScene();
}

void Renderer::showAllFlight() {
  showingFullFlight = true;
  
  // Hide stopwatch when switching to full flight mode
  hideStopwatch();
  
  // Restore full blackbox data
  blackboxAircraftStates = fullBlackboxAircraftStates;
  
  // Update the blackbox rendering
  updateBlackboxForCurrentTest();
  
  // Re-render using in-memory data (no S3 fetch)
  renderFullScene();
}

void Renderer::createHighlightedFlightTapes(vec3 offset) {
  if (fullBlackboxAircraftStates.empty() || testSpans.empty()) {
    return;
  }
  
  // Remove the empty data we added earlier
  this->blackboxHighlightTapes->RemoveAllInputs();
  
  // Create dimmed version of full flight (25% brightness)
  std::vector<vec3> fullPoints = stateToVector(fullBlackboxAircraftStates);
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
        std::vector<vec3> spanPoints = stateToVector(spanStates);
        this->blackboxHighlightTapes->AddInputData(createTapeSet(offset, spanPoints, stateToOrientation(spanStates)));
      }
    }
  }
  
  // Ensure highlight actor is at full brightness
  blackboxHighlightActor->GetProperty()->SetOpacity(1.0);
  blackboxHighlightActor->GetBackfaceProperty()->SetOpacity(1.0);
}

void Renderer::createStopwatch() {
  // Create stopwatch geometry (circle with tick marks and hands)
  vtkNew<vtkPoints> points;
  vtkNew<vtkCellArray> lines;
  
  // Clock face circle (100px diameter)
  const scalar radius = static_cast<scalar>(50.0f);
  const int numPoints = 60;
  for (int i = 0; i < numPoints; i++) {
    scalar angle = static_cast<scalar>(2.0 * M_PI * i / numPoints);
    points->InsertNextPoint(radius * cos(angle), radius * sin(angle), 0);
  }
  
  // Add tick marks (10 ticks)
  for (int i = 0; i < 10; i++) {
    scalar angle = static_cast<scalar>(2.0 * M_PI * i / 10);
    scalar innerRadius = radius * static_cast<scalar>(0.85f);
    scalar outerRadius = radius * static_cast<scalar>(0.95f);
    
    // Inner point
    points->InsertNextPoint(innerRadius * cos(angle), innerRadius * sin(angle), 0);
    // Outer point
    points->InsertNextPoint(outerRadius * cos(angle), outerRadius * sin(angle), 0);
  }
  
  // Circle lines
  vtkNew<vtkPolyLine> circle;
  circle->GetPointIds()->SetNumberOfIds(numPoints + 1);
  for (int i = 0; i < numPoints; i++) {
    circle->GetPointIds()->SetId(i, i);
  }
  circle->GetPointIds()->SetId(numPoints, 0); // Close the circle
  lines->InsertNextCell(circle);
  
  // Tick mark lines
  for (int i = 0; i < 10; i++) {
    vtkNew<vtkLine> tick;
    tick->GetPointIds()->SetId(0, numPoints + i * 2);
    tick->GetPointIds()->SetId(1, numPoints + i * 2 + 1);
    lines->InsertNextCell(tick);
  }
  
  // Create polydata
  vtkNew<vtkPolyData> stopwatchData;
  stopwatchData->SetPoints(points);
  stopwatchData->SetLines(lines);
  
  // Create mapper and actor
  vtkNew<vtkPolyDataMapper2D> mapper;
  mapper->SetInputData(stopwatchData);
  
  stopwatchActor = vtkSmartPointer<vtkActor2D>::New();
  stopwatchActor->SetMapper(mapper);
  stopwatchActor->GetProperty()->SetColor(1.0, 1.0, 1.0); // White
  stopwatchActor->GetProperty()->SetOpacity(0.8);
  
  // Position will be set in updateStopwatch using screen coordinates
  stopwatchActor->SetPosition(0, 0);
  
  // Create time text actor
  stopwatchTimeActor = vtkSmartPointer<vtkTextActor>::New();
  stopwatchTimeActor->GetTextProperty()->SetFontSize(12);
  stopwatchTimeActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
  stopwatchTimeActor->GetTextProperty()->SetJustificationToCentered();
}

void Renderer::createControlsOverlay() {
  // Outline (stick box, crosshair, throttle outline)
  controlOutlineActor = vtkSmartPointer<vtkActor2D>::New();
  vtkNew<vtkPolyDataMapper2D> outlineMapper;
  outlineMapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  controlOutlineActor->SetMapper(outlineMapper);
  controlOutlineActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
  controlOutlineActor->GetProperty()->SetOpacity(0.8);

  // Moving stick indicator
  controlStickActor = vtkSmartPointer<vtkActor2D>::New();
  vtkNew<vtkPolyDataMapper2D> stickMapper;
  stickMapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  controlStickActor->SetMapper(stickMapper);
  controlStickActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
  controlStickActor->GetProperty()->SetOpacity(0.8);

  // Throttle fill
  throttleFillActor = vtkSmartPointer<vtkActor2D>::New();
  vtkNew<vtkPolyDataMapper2D> throttleMapper;
  throttleMapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  throttleFillActor->SetMapper(throttleMapper);
  throttleFillActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
  throttleFillActor->GetProperty()->SetOpacity(0.5);

  controlSourceActor = vtkSmartPointer<vtkTextActor>::New();
  controlSourceActor->GetTextProperty()->SetFontSize(12);
  controlSourceActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
  controlSourceActor->GetTextProperty()->SetJustificationToLeft();

  // Attitude indicator (sky/ground fill + outline/horizon)
  attitudeSkyActor = vtkSmartPointer<vtkActor2D>::New();
  vtkNew<vtkPolyDataMapper2D> skyMapper;
  skyMapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  attitudeSkyActor->SetMapper(skyMapper);
  attitudeSkyActor->GetProperty()->SetColor(0.2, 0.4, 0.9); // blue sky
  attitudeSkyActor->GetProperty()->SetOpacity(0.7);

  attitudeGroundActor = vtkSmartPointer<vtkActor2D>::New();
  vtkNew<vtkPolyDataMapper2D> groundMapper;
  groundMapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  attitudeGroundActor->SetMapper(groundMapper);
  attitudeGroundActor->GetProperty()->SetColor(0.8, 0.2, 0.2); // red ground
  attitudeGroundActor->GetProperty()->SetOpacity(0.7);

  attitudeOutlineActor = vtkSmartPointer<vtkActor2D>::New();
  vtkNew<vtkPolyDataMapper2D> attitudeOutlineMapper;
  attitudeOutlineMapper->SetInputData(vtkSmartPointer<vtkPolyData>::New());
  attitudeOutlineActor->SetMapper(attitudeOutlineMapper);
  attitudeOutlineActor->GetProperty()->SetColor(1.0, 1.0, 1.0);
  attitudeOutlineActor->GetProperty()->SetOpacity(0.8);

  velocityActor = vtkSmartPointer<vtkTextActor>::New();
  velocityActor->GetTextProperty()->SetFontSize(12);
  velocityActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
  velocityActor->GetTextProperty()->SetJustificationToLeft();
}

void Renderer::updateStopwatch(gp_scalar currentTime) {
  if (!stopwatchActor || !stopwatchTimeActor) return;
  
  stopwatchTime = currentTime;
  
  // Position stopwatch in upper right corner (2D screen coordinates)
  int* windowSize = renderWindow->GetSize();
  scalar centerX = static_cast<scalar>(windowSize[0] - 70); // 70px from right edge
  scalar centerY = static_cast<scalar>(windowSize[1] - 70); // 70px from top
  
  // Create updated stopwatch geometry with hands
  vtkNew<vtkPoints> points;
  vtkNew<vtkCellArray> lines;
  
  // Clock face circle
  const scalar radius = static_cast<scalar>(45.0f);
  const int numPoints = 60;
  for (int i = 0; i < numPoints; i++) {
    scalar angle = static_cast<scalar>(2.0 * M_PI * i / numPoints);
    points->InsertNextPoint(centerX + radius * cos(angle), centerY + radius * sin(angle), 0);
  }
  
  // Add tick marks (10 ticks)
  for (int i = 0; i < 10; i++) {
    scalar angle = static_cast<scalar>(2.0 * M_PI * i / 10 - M_PI/2); // Start from top
    scalar innerRadius = radius * static_cast<scalar>(0.85f);
    scalar outerRadius = radius * static_cast<scalar>(0.95f);
    
    points->InsertNextPoint(centerX + innerRadius * cos(angle), centerY + innerRadius * sin(angle), 0);
    points->InsertNextPoint(centerX + outerRadius * cos(angle), centerY + outerRadius * sin(angle), 0);
  }
  
  // Add hand center point
  points->InsertNextPoint(centerX, centerY, 0);
  int centerPointId = numPoints + 20;
  
  // Calculate hand angles (clockwise rotation)
  scalar secondAngle = static_cast<scalar>(-fmod(currentTime, static_cast<gp_scalar>(1.0f)) * 2.0 * M_PI + M_PI/2); // 1 rev/sec, clockwise from top
  scalar tenSecAngle = static_cast<scalar>(-fmod(currentTime, static_cast<gp_scalar>(10.0f)) * 2.0 * M_PI / 10.0 + M_PI/2); // 1 rev/10sec, clockwise
  
  // 1-second hand (longer, thinner)
  scalar secondRadius = radius * static_cast<scalar>(0.8f);
  points->InsertNextPoint(centerX + secondRadius * cos(secondAngle), centerY + secondRadius * sin(secondAngle), 0);
  int secondHandId = centerPointId + 1;
  
  // 10-second hand (shorter, thicker conceptually)
  scalar tenSecRadius = radius * static_cast<scalar>(0.6f);
  points->InsertNextPoint(centerX + tenSecRadius * cos(tenSecAngle), centerY + tenSecRadius * sin(tenSecAngle), 0);
  int tenSecHandId = centerPointId + 2;
  
  // Create circle lines
  vtkNew<vtkPolyLine> circle;
  circle->GetPointIds()->SetNumberOfIds(numPoints + 1);
  for (int i = 0; i < numPoints; i++) {
    circle->GetPointIds()->SetId(i, i);
  }
  circle->GetPointIds()->SetId(numPoints, 0);
  lines->InsertNextCell(circle);
  
  // Create tick mark lines
  for (int i = 0; i < 10; i++) {
    vtkNew<vtkLine> tick;
    tick->GetPointIds()->SetId(0, numPoints + i * 2);
    tick->GetPointIds()->SetId(1, numPoints + i * 2 + 1);
    lines->InsertNextCell(tick);
  }
  
  // Create hand lines
  vtkNew<vtkLine> secondHand;
  secondHand->GetPointIds()->SetId(0, centerPointId);
  secondHand->GetPointIds()->SetId(1, secondHandId);
  lines->InsertNextCell(secondHand);
  
  vtkNew<vtkLine> tenSecHand;
  tenSecHand->GetPointIds()->SetId(0, centerPointId);
  tenSecHand->GetPointIds()->SetId(1, tenSecHandId);
  lines->InsertNextCell(tenSecHand);
  
  // Update actor with new geometry
  vtkNew<vtkPolyData> stopwatchData;
  stopwatchData->SetPoints(points);
  stopwatchData->SetLines(lines);
  
  vtkPolyDataMapper2D* mapper = vtkPolyDataMapper2D::SafeDownCast(stopwatchActor->GetMapper());
  if (mapper) {
    mapper->SetInputData(stopwatchData);
  }
  
  // Update time text (centered horizontally below hands center)
  char timeStr[32];
  snprintf(timeStr, sizeof(timeStr), "%.2f", currentTime);
  stopwatchTimeActor->SetInput(timeStr);
  stopwatchTimeActor->SetPosition(centerX - 12, centerY - 25); // Centered horizontally, positioned below hands center
}

void Renderer::updateStopwatchPosition() {
  if (!stopwatchActor || !stopwatchTimeActor || !stopwatchVisible) {
    return;
  }
  
  // Recalculate position based on current window size
  int* windowSize = renderWindow->GetSize();
  scalar centerX = static_cast<scalar>(windowSize[0] - 70); // 70px from right edge
  scalar centerY = static_cast<scalar>(windowSize[1] - 70); // 70px from top
  
  // Update the stopwatch geometry with new center position
  // This requires recreating the geometry since the center position has changed
  if (stopwatchTime >= 0.0) { // Only update if we have valid stopwatch time
    updateStopwatch(stopwatchTime); // This will recreate the geometry with new position
  }

  // Keep controls aligned with the stopwatch box
  updateControlsPosition();
}

bool Renderer::getControlStateAtTime(gp_scalar currentSimTime, gp_scalar& pitch, gp_scalar& roll, gp_scalar& throttle, int arenaIndex, bool& usedBlackbox, const AircraftState*& chosenState) {
  const std::vector<AircraftState>* states = nullptr;
  bool useBlackboxTime = false;
  usedBlackbox = false;
  chosenState = nullptr;

  if (!blackboxAircraftStates.empty()) {
    states = &blackboxAircraftStates;
    useBlackboxTime = true;
    usedBlackbox = true;
  } else if (!evalResults.aircraftStateList.empty()) {
    if (arenaIndex < 0) arenaIndex = 0;
    if (arenaIndex >= static_cast<int>(evalResults.aircraftStateList.size())) {
      arenaIndex = static_cast<int>(evalResults.aircraftStateList.size()) - 1;
    }
    if (!evalResults.aircraftStateList[arenaIndex].empty()) {
      states = &evalResults.aircraftStateList[arenaIndex];
    }
  } else {
    return false;
  }

  if (!states || states->empty()) {
    return false;
  }

  scalar startTime = useBlackboxTime
    ? static_cast<scalar>((*states)[0].getSimTimeMsec()) / static_cast<scalar>(1000000.0f)
    : static_cast<scalar>((*states)[0].getSimTimeMsec()) / static_cast<scalar>(1000.0f);

  size_t chosenIndex = 0;
  for (size_t i = 0; i < states->size(); ++i) {
    scalar stateTime = useBlackboxTime
      ? static_cast<scalar>((*states)[i].getSimTimeMsec()) / static_cast<scalar>(1000000.0f)
      : static_cast<scalar>((*states)[i].getSimTimeMsec()) / static_cast<scalar>(1000.0f);
    scalar relativeTime = stateTime - startTime;
    if (relativeTime > currentSimTime) {
      break;
    }
    chosenIndex = i;
  }

  const AircraftState& state = (*states)[chosenIndex];
  pitch = state.getPitchCommand();
  roll = state.getRollCommand();
  throttle = state.getThrottleCommand();
  chosenState = &state;

  // Debug output for xiao mode RC commands
  static int debugCounter = 0;
  if (usedBlackbox && debugCounter++ % 30 == 0) {  // Print every 30th call
    std::cout << "RC at t=" << currentSimTime << "s: pitch=" << pitch
              << " roll=" << roll << " throttle=" << throttle << std::endl;
  }

  return true;
}

void Renderer::updateControlsOverlay(gp_scalar currentTime) {
  if (!renderWindow || !controlOutlineActor || !controlStickActor || !throttleFillActor || !controlsVisible) {
    return;
  }

  int selectedArena = focusMode ? focusArenaIndex : 0;
  bool usedBlackbox = false;
  const AircraftState* chosenState = nullptr;

  gp_scalar pitch = lastControlPitch;
  gp_scalar roll = lastControlRoll;
  gp_scalar throttle = lastControlThrottle;
  bool haveControls = getControlStateAtTime(currentTime, pitch, roll, throttle, selectedArena, usedBlackbox, chosenState);

  lastControlPitch = pitch;
  lastControlRoll = roll;
  lastControlThrottle = throttle;
  lastControlsTime = currentTime;

  int* windowSize = renderWindow->GetSize();
  scalar clockCenterX = static_cast<scalar>(windowSize[0] - 70); // Align with stopwatch anchor
  scalar clockCenterY = static_cast<scalar>(windowSize[1] - 70);

  scalar stickSize = static_cast<scalar>(70.0f);
  scalar stickHalf = stickSize / static_cast<scalar>(2.0f);
  scalar velocityAnchorX = clockCenterX - static_cast<scalar>(80.0f);
  scalar attitudeCenterX = velocityAnchorX - static_cast<scalar>(80.0f);
  scalar attitudeRadius = static_cast<scalar>(30.0f);
  scalar throttleLeft = attitudeCenterX - static_cast<scalar>(70.0f);
  scalar throttleWidth = static_cast<scalar>(12.0f);
  scalar throttleRight = throttleLeft + throttleWidth;
  scalar stickCenterX = throttleLeft - static_cast<scalar>(90.0f);
  scalar stickCenterY = clockCenterY;
  scalar stickLeft = stickCenterX - stickHalf;
  scalar stickRight = stickCenterX + stickHalf;
  scalar stickBottom = stickCenterY - stickHalf;
  scalar stickTop = stickCenterY + stickHalf;

  scalar throttleBottom = stickBottom;
  scalar throttleTop = stickTop;
  scalar attitudeCenterY = clockCenterY;

  vtkNew<vtkPoints> outlinePoints;
  vtkNew<vtkCellArray> outlineLines;
  auto addLine = [&](scalar x1, scalar y1, scalar x2, scalar y2) {
    vtkIdType id1 = outlinePoints->InsertNextPoint(x1, y1, 0);
    vtkIdType id2 = outlinePoints->InsertNextPoint(x2, y2, 0);
    vtkNew<vtkLine> line;
    line->GetPointIds()->SetId(0, id1);
    line->GetPointIds()->SetId(1, id2);
    outlineLines->InsertNextCell(line);
  };

  // Stick outline and crosshair
  addLine(stickLeft, stickBottom, stickRight, stickBottom);
  addLine(stickRight, stickBottom, stickRight, stickTop);
  addLine(stickRight, stickTop, stickLeft, stickTop);
  addLine(stickLeft, stickTop, stickLeft, stickBottom);
  addLine(stickCenterX, stickBottom, stickCenterX, stickTop);
  addLine(stickLeft, stickCenterY, stickRight, stickCenterY);

  // Throttle outline
  addLine(throttleLeft, throttleBottom, throttleRight, throttleBottom);
  addLine(throttleRight, throttleBottom, throttleRight, throttleTop);
  addLine(throttleRight, throttleTop, throttleLeft, throttleTop);
  addLine(throttleLeft, throttleTop, throttleLeft, throttleBottom);

  vtkNew<vtkPolyData> outlineData;
  outlineData->SetPoints(outlinePoints);
  outlineData->SetLines(outlineLines);
  vtkPolyDataMapper2D* outlineMapper = vtkPolyDataMapper2D::SafeDownCast(controlOutlineActor->GetMapper());
  if (outlineMapper) {
    outlineMapper->SetInputData(outlineData);
  }

  // Stick indicator (clamp to box). Positive pitch = stick back (downwards on radio face).
  roll = CLAMP_DEF(roll, -1.0f, 1.0f);
  pitch = CLAMP_DEF(pitch, -1.0f, 1.0f);
  throttle = CLAMP_DEF(throttle, -1.0f, 1.0f);

  scalar knobX = stickCenterX + roll * stickHalf;
  scalar knobY = stickCenterY - pitch * stickHalf;
  scalar knobRadius = static_cast<scalar>(6.0f);

  vtkNew<vtkPoints> stickPoints;
  vtkNew<vtkCellArray> stickLines;
  vtkNew<vtkPolyLine> knob;
  const int knobSegments = 20;
  knob->GetPointIds()->SetNumberOfIds(knobSegments + 1);
  for (int i = 0; i < knobSegments; ++i) {
    scalar angle = static_cast<scalar>(2.0 * M_PI * i / knobSegments);
    vtkIdType id = stickPoints->InsertNextPoint(knobX + knobRadius * cos(angle), knobY + knobRadius * sin(angle), 0);
    knob->GetPointIds()->SetId(i, id);
  }
  knob->GetPointIds()->SetId(knobSegments, 0);
  stickLines->InsertNextCell(knob);

  vtkNew<vtkPolyData> stickData;
  stickData->SetPoints(stickPoints);
  stickData->SetLines(stickLines);
  vtkPolyDataMapper2D* stickMapper = vtkPolyDataMapper2D::SafeDownCast(controlStickActor->GetMapper());
  if (stickMapper) {
    stickMapper->SetInputData(stickData);
  }

  // Throttle fill (bottom = idle)
  scalar throttleFraction = CLAMP_DEF((throttle + static_cast<scalar>(1.0f)) / static_cast<scalar>(2.0f), static_cast<scalar>(0.0f), static_cast<scalar>(1.0f));
  scalar fillTop = throttleBottom + (throttleTop - throttleBottom) * throttleFraction;

  vtkNew<vtkPoints> fillPoints;
  fillPoints->InsertNextPoint(throttleLeft, throttleBottom, 0);
  fillPoints->InsertNextPoint(throttleRight, throttleBottom, 0);
  fillPoints->InsertNextPoint(throttleRight, fillTop, 0);
  fillPoints->InsertNextPoint(throttleLeft, fillTop, 0);

  vtkNew<vtkCellArray> fillPolys;
  vtkNew<vtkPolygon> fillPoly;
  fillPoly->GetPointIds()->SetNumberOfIds(4);
  for (int i = 0; i < 4; ++i) {
    fillPoly->GetPointIds()->SetId(i, i);
  }
  fillPolys->InsertNextCell(fillPoly);

  vtkNew<vtkPolyData> fillData;
  fillData->SetPoints(fillPoints);
  fillData->SetPolys(fillPolys);
  vtkPolyDataMapper2D* throttleMapper = vtkPolyDataMapper2D::SafeDownCast(throttleFillActor->GetMapper());
  if (throttleMapper) {
    throttleMapper->SetInputData(fillData);
  }

  // Attitude indicator derived from quaternion
  gp_scalar attRoll = static_cast<gp_scalar>(0.0f);
  gp_scalar attPitch = static_cast<gp_scalar>(0.0f);
  if (chosenState) {
    // earth->body quaternion rotates body frame vectors to earth frame (NED)
    gp_quat q = chosenState->getOrientation();
    gp_vec3 forward = q * gp_vec3::UnitX();   // body +X -> earth frame
    gp_vec3 right = q * gp_vec3::UnitY();     // body +Y -> earth frame
    gp_vec3 up = q * -gp_vec3::UnitZ();       // body -Z (up) -> earth frame
    attPitch = std::atan2(-forward[2], std::sqrt(forward[0] * forward[0] + forward[1] * forward[1]));
    attRoll = std::atan2(right[2], -up[2]); // negate up-Z so level attitude yields 0 roll
  }
  scalar maxPitch = static_cast<scalar>(M_PI / 3.0); // clamp to +-60 deg for display
  // Positive pitch (nose up) should move horizon down to show more sky
  scalar pitchOffset = -CLAMP_DEF(attPitch / maxPitch, static_cast<scalar>(-1.0f), static_cast<scalar>(1.0f)) * attitudeRadius;
  scalar drawRoll = -attRoll; // horizon moves opposite aircraft roll
  scalar cosR = cos(drawRoll);
  scalar sinR = sin(drawRoll);

  auto rotateToLevel = [&](scalar x, scalar y) {
    scalar rx = x * cosR - y * sinR;
    scalar ry = x * sinR + y * cosR;
    return std::pair<scalar, scalar>(rx, ry);
  };
  auto rotateFromLevel = [&](scalar x, scalar y) {
    scalar rx = x * cosR + y * sinR;
    scalar ry = -x * sinR + y * cosR;
    return std::pair<scalar, scalar>(rx, ry);
  };

  // Precompute circle points
  const int attSegments = 80;
  std::vector<std::pair<scalar, scalar>> circlePoints;
  circlePoints.reserve(attSegments);
  for (int i = 0; i < attSegments; ++i) {
    scalar angle = static_cast<scalar>(2.0 * M_PI * i / attSegments);
    circlePoints.emplace_back(attitudeRadius * cos(angle), attitudeRadius * sin(angle));
  }

  auto clipHalf = [&](bool keepSky) -> vtkSmartPointer<vtkPolyData> {
    if (circlePoints.empty()) {
      return vtkSmartPointer<vtkPolyData>::New();
    }
    std::vector<std::pair<scalar, scalar>> input = circlePoints;
    std::vector<std::pair<scalar, scalar>> output;

    auto inside = [&](const std::pair<scalar, scalar>& p) {
      auto r = rotateToLevel(p.first, p.second);
      return keepSky ? (r.second >= pitchOffset) : (r.second <= pitchOffset);
    };

    auto intersect = [&](const std::pair<scalar, scalar>& p1, const std::pair<scalar, scalar>& p2) {
      auto r1 = rotateToLevel(p1.first, p1.second);
      auto r2 = rotateToLevel(p2.first, p2.second);
      if (r1.second == r2.second) return p1;
      scalar t = (pitchOffset - r1.second) / (r2.second - r1.second);
      scalar nx = p1.first + t * (p2.first - p1.first);
      scalar ny = p1.second + t * (p2.second - p1.second);
      return std::make_pair(nx, ny);
    };

    for (size_t i = 0; i < input.size(); ++i) {
      auto curr = input[i];
      auto next = input[(i + 1) % input.size()];
      bool currIn = inside(curr);
      bool nextIn = inside(next);

      if (currIn && nextIn) {
        output.push_back(next);
      } else if (currIn && !nextIn) {
        output.push_back(intersect(curr, next));
      } else if (!currIn && nextIn) {
        output.push_back(intersect(curr, next));
        output.push_back(next);
      }
    }

    vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New();
    if (output.size() < 3) {
      return data;
    }

    vtkNew<vtkPoints> pts;
    vtkNew<vtkPolygon> poly;
    poly->GetPointIds()->SetNumberOfIds(static_cast<vtkIdType>(output.size()));
    for (size_t i = 0; i < output.size(); ++i) {
      auto p = rotateFromLevel(output[i].first, output[i].second);
      vtkIdType id = pts->InsertNextPoint(attitudeCenterX + p.first, attitudeCenterY + p.second, 0);
      poly->GetPointIds()->SetId(static_cast<vtkIdType>(i), id);
    }
    vtkNew<vtkCellArray> polys;
    polys->InsertNextCell(poly);
    data->SetPoints(pts);
    data->SetPolys(polys);
    return data;
  };

  vtkPolyDataMapper2D* skyMap = vtkPolyDataMapper2D::SafeDownCast(attitudeSkyActor->GetMapper());
  vtkPolyDataMapper2D* groundMap = vtkPolyDataMapper2D::SafeDownCast(attitudeGroundActor->GetMapper());
  if (skyMap) skyMap->SetInputData(clipHalf(true));
  if (groundMap) groundMap->SetInputData(clipHalf(false));

  // Outline only (horizon implied by sky/ground split)
  vtkNew<vtkPoints> attPoints;
  vtkNew<vtkCellArray> attLines;
  vtkNew<vtkPolyLine> attCircle;
  attCircle->GetPointIds()->SetNumberOfIds(attSegments + 1);
  for (int i = 0; i < attSegments; ++i) {
    vtkIdType id = attPoints->InsertNextPoint(attitudeCenterX + circlePoints[i].first,
      attitudeCenterY + circlePoints[i].second, 0);
    attCircle->GetPointIds()->SetId(i, id);
  }
  attCircle->GetPointIds()->SetId(attSegments, 0);
  attLines->InsertNextCell(attCircle);

  vtkNew<vtkPolyData> attOutlineData;
  attOutlineData->SetPoints(attPoints);
  attOutlineData->SetLines(attLines);
  vtkPolyDataMapper2D* attOutlineMapper = vtkPolyDataMapper2D::SafeDownCast(attitudeOutlineActor->GetMapper());
  if (attOutlineMapper) {
    attOutlineMapper->SetInputData(attOutlineData);
  }

  // Fixed aircraft reference symbol at center (small crosshair)
  {
    vtkNew<vtkPoints> refPoints;
    vtkNew<vtkCellArray> refLines;
    scalar wingSpan = attitudeRadius * static_cast<scalar>(0.6f);
    scalar tailSpan = attitudeRadius * static_cast<scalar>(0.3f);
    vtkIdType w1 = refPoints->InsertNextPoint(attitudeCenterX - wingSpan, attitudeCenterY, 0);
    vtkIdType w2 = refPoints->InsertNextPoint(attitudeCenterX + wingSpan, attitudeCenterY, 0);
    vtkIdType t1 = refPoints->InsertNextPoint(attitudeCenterX, attitudeCenterY - tailSpan, 0);
    vtkIdType t2 = refPoints->InsertNextPoint(attitudeCenterX, attitudeCenterY + tailSpan, 0);
    vtkNew<vtkLine> wings;
    wings->GetPointIds()->SetId(0, w1);
    wings->GetPointIds()->SetId(1, w2);
    vtkNew<vtkLine> tail;
    tail->GetPointIds()->SetId(0, t1);
    tail->GetPointIds()->SetId(1, t2);
    refLines->InsertNextCell(wings);
    refLines->InsertNextCell(tail);
    vtkNew<vtkPolyData> refData;
    refData->SetPoints(refPoints);
    refData->SetLines(refLines);
    if (attOutlineMapper) {
      // Overlay onto same outline actor for simplicity
      vtkNew<vtkAppendPolyData> append;
      append->AddInputData(attOutlineData);
      append->AddInputData(refData);
      append->Update();
      attOutlineMapper->SetInputData(append->GetOutput());
    }
  }

  // Update arena label
  if (controlSourceActor) {
    int fontSize = std::max(8, std::min(16, static_cast<int>(windowSize[1] * 0.018)));
    controlSourceActor->GetTextProperty()->SetFontSize(fontSize);
    std::ostringstream label;
    if (!haveControls) {
      label << "No control data";
    } else if (usedBlackbox) {
      if (inXiaoMode) {
        label << "Xiao";
      } else {
        label << "Blackbox";
      }
    } else {
      ScenarioMetadata meta{};
      if (!evalResults.scenarioList.empty() && selectedArena >= 0 && selectedArena < static_cast<int>(evalResults.scenarioList.size())) {
        meta = evalResults.scenarioList[selectedArena];
      } else {
        meta = evalResults.scenario;
      }
      int pathIdx = (meta.pathVariantIndex >= 0) ? meta.pathVariantIndex : selectedArena;
      int windIdx = meta.windVariantIndex;
      label << "Path " << pathIdx << " / Wind " << windIdx;
    }
    controlSourceActor->SetInput(label.str().c_str());
    controlSourceActor->SetPosition(stickLeft, stickTop + static_cast<scalar>(12.0f));
  }

  // Velocity indicator (approximate forward speed in m/s)
  if (velocityActor) {
    int fontSize = std::max(8, std::min(16, static_cast<int>(windowSize[1] * 0.018)));
    velocityActor->GetTextProperty()->SetFontSize(fontSize);
    gp_scalar speed = 0.0f;
    if (haveControls && chosenState) {
      gp_vec3 vel = chosenState->getVelocity();
      speed = vel.norm();
    }
    std::ostringstream velStr;
    velStr << "Vel " << std::fixed << std::setprecision(1) << speed << " m/s";
    velocityActor->SetInput(velStr.str().c_str());
    velocityActor->SetPosition(attitudeCenterX - attitudeRadius,
      attitudeCenterY + attitudeRadius + static_cast<scalar>(10.0f));
  }
}

void Renderer::updateControlsPosition() {
  if (!controlsVisible) {
    return;
  }
  updateControlsOverlay(lastControlsTime);
}

void Renderer::toggleFocusMode() {
  if (evalResults.pathList.empty()) {
    return;
  }

  focusMode = true;

  if (focusArenaIndex < 0 || focusArenaIndex >= static_cast<int>(evalResults.pathList.size())) {
    focusArenaIndex = 0;
  }
  setFocusArena(focusArenaIndex);
  renderWindow->Render();
}

void Renderer::adjustFocusArena(int delta) {
  if (evalResults.pathList.empty()) {
    return;
  }
  focusMode = true;
  int total = static_cast<int>(evalResults.pathList.size());
  int newIndex = focusArenaIndex + delta;
  if (total > 0) {
    newIndex = (newIndex % total + total) % total;
    setFocusArena(newIndex);
    renderWindow->Render();
  }
}

void Renderer::focusMoveLeft() {
  if (evalResults.pathList.empty()) {
    return;
  }
  int total = static_cast<int>(evalResults.pathList.size());
  if (total <= 1) {
    return;
  }
  focusMode = true;
  int columns = static_cast<int>(std::ceil(std::sqrt(total)));
  int rows = static_cast<int>(std::ceil(static_cast<scalar>(total) / static_cast<scalar>(columns)));
  int row = focusArenaIndex / columns;
  int col = focusArenaIndex % columns;
  int startRow = row;
  do {
    row = (row + rows - 1) % rows;
    int candidate = row * columns + col;
    if (candidate < total) {
      focusArenaIndex = candidate;
      setFocusArena(focusArenaIndex);
      renderWindow->Render();
      return;
    }
  } while (row != startRow);
  setFocusArena(focusArenaIndex);
  renderWindow->Render();
}

void Renderer::focusMoveRight() {
  if (evalResults.pathList.empty()) {
    return;
  }
  int total = static_cast<int>(evalResults.pathList.size());
  if (total <= 1) {
    return;
  }
  focusMode = true;
  int columns = static_cast<int>(std::ceil(std::sqrt(total)));
  int rows = static_cast<int>(std::ceil(static_cast<scalar>(total) / static_cast<scalar>(columns)));
  int row = focusArenaIndex / columns;
  int col = focusArenaIndex % columns;
  int startRow = row;
  do {
    row = (row + 1) % rows;
    int candidate = row * columns + col;
    if (candidate < total) {
      focusArenaIndex = candidate;
      setFocusArena(focusArenaIndex);
      renderWindow->Render();
      return;
    }
  } while (row != startRow);
  setFocusArena(focusArenaIndex);
  renderWindow->Render();
}

void Renderer::focusMoveUp() {
  if (evalResults.pathList.empty()) {
    return;
  }
  int total = static_cast<int>(evalResults.pathList.size());
  if (total <= 1) {
    return;
  }
  focusMode = true;
  int columns = static_cast<int>(std::ceil(std::sqrt(total)));
  int row = focusArenaIndex / columns;
  int col = focusArenaIndex % columns;
  int startCol = col;
  do {
    col = (col + 1) % columns;
    int candidate = row * columns + col;
    if (candidate < total) {
      focusArenaIndex = candidate;
      setFocusArena(focusArenaIndex);
      renderWindow->Render();
      return;
    }
  } while (col != startCol);
  setFocusArena(focusArenaIndex);
  renderWindow->Render();
}

void Renderer::focusMoveDown() {
  if (evalResults.pathList.empty()) {
    return;
  }
  int total = static_cast<int>(evalResults.pathList.size());
  if (total <= 1) {
    return;
  }
  focusMode = true;
  int columns = static_cast<int>(std::ceil(std::sqrt(total)));
  int row = focusArenaIndex / columns;
  int col = focusArenaIndex % columns;
  int startCol = col;
  do {
    col = (col + columns - 1) % columns;
    int candidate = row * columns + col;
    if (candidate < total) {
      focusArenaIndex = candidate;
      setFocusArena(focusArenaIndex);
      renderWindow->Render();
      return;
    }
  } while (col != startCol);
  setFocusArena(focusArenaIndex);
  renderWindow->Render();
}

void Renderer::setFocusArena(int arenaIdx) {
  if (evalResults.pathList.empty()) {
    return;
  }

  int total = static_cast<int>(evalResults.pathList.size());
  if (total == 0) {
    return;
  }

  arenaIdx = (arenaIdx % total + total) % total;
  focusArenaIndex = arenaIdx;

  vtkRenderer* activeRenderer = renderWindow->GetRenderers()->GetFirstRenderer();
  if (!activeRenderer) {
    return;
  }

  vec3 offset = renderingOffset(arenaIdx);
  scalar focusZ = 0.0f;
  if (arenaIdx < static_cast<int>(evalResults.pathList.size()) && !evalResults.pathList[arenaIdx].empty()) {
    focusZ = evalResults.pathList[arenaIdx].front().start[2];
  }

  scalar camX = offset[0] - static_cast<scalar>(70.0f);
  scalar camY = offset[1] - static_cast<scalar>(10.0f);
  scalar camZ = offset[2] - static_cast<scalar>(100.0f);

  focusCameraPosition = {camX, camY, camZ};
  focusCameraFocalPoint = {offset[0], offset[1], offset[2] - static_cast<scalar>(10.0f)};
  focusCameraViewUp = {0.0f, 0.0f, -1.0f};

  vtkCamera* camera = activeRenderer->GetActiveCamera();
  if (camera) {
    camera->SetPosition(focusCameraPosition[0], focusCameraPosition[1], focusCameraPosition[2]);
    camera->SetFocalPoint(focusCameraFocalPoint[0], focusCameraFocalPoint[1], focusCameraFocalPoint[2]);
    camera->SetViewUp(focusCameraViewUp[0], focusCameraViewUp[1], focusCameraViewUp[2]);
    activeRenderer->ResetCameraClippingRange();
  }
}

void Renderer::hideStopwatch() {
  if (stopwatchActor && stopwatchTimeActor) {
    vtkRenderer* renderer = vtkRenderer::SafeDownCast(renderWindow->GetRenderers()->GetFirstRenderer());
    if (renderer) {
      renderer->RemoveActor2D(stopwatchActor);
      renderer->RemoveActor2D(stopwatchTimeActor);
      if (controlOutlineActor && controlStickActor && throttleFillActor) {
        renderer->RemoveActor2D(controlOutlineActor);
        renderer->RemoveActor2D(controlStickActor);
        renderer->RemoveActor2D(throttleFillActor);
        if (controlSourceActor) {
          renderer->RemoveActor2D(controlSourceActor);
        }
      }
      if (attitudeSkyActor && attitudeGroundActor && attitudeOutlineActor) {
        renderer->RemoveActor2D(attitudeSkyActor);
        renderer->RemoveActor2D(attitudeGroundActor);
        renderer->RemoveActor2D(attitudeOutlineActor);
      }
      if (velocityActor) {
        renderer->RemoveActor2D(velocityActor);
      }
    }
    stopwatchVisible = false;
    controlsVisible = false;
  }
}

void Renderer::togglePlaybackAnimation() {
  if (isPlaybackActive) {
    // Stop animation
    isPlaybackActive = false;
    isPlaybackPaused = false;
    totalPausedTime = std::chrono::duration<gp_scalar>::zero();
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
    totalPausedTime = std::chrono::duration<gp_scalar>::zero();
    animationStartTime = std::chrono::steady_clock::now();
    stopwatchVisible = true;
    controlsVisible = true;
    
    // Show stopwatch actors (need to get the renderer from the render window)
    vtkRenderer* renderer = vtkRenderer::SafeDownCast(renderWindow->GetRenderers()->GetFirstRenderer());
  if (renderer) {
    renderer->AddActor2D(stopwatchActor);
    renderer->AddActor2D(stopwatchTimeActor);
    renderer->AddActor2D(controlOutlineActor);
    renderer->AddActor2D(controlStickActor);
    renderer->AddActor2D(throttleFillActor);
    renderer->AddActor2D(controlSourceActor);
    renderer->AddActor2D(attitudeSkyActor);
    renderer->AddActor2D(attitudeGroundActor);
    renderer->AddActor2D(attitudeOutlineActor);
    renderer->AddActor2D(velocityActor);
  }
    
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
  scalar elapsed = static_cast<scalar>(std::chrono::duration_cast<std::chrono::milliseconds>(effectiveElapsed).count()) / static_cast<scalar>(1000.0f);
  
  // Path always drives the timing - find longest path duration across all arenas
  scalar primaryDuration = 0.0f;
  int maxArenas = evalResults.pathList.size();
  for (int i = 0; i < maxArenas; i++) {
    if (!evalResults.pathList[i].empty()) {
      // Use path timestamp directly (convert from milliseconds to seconds)
      scalar pathDuration = static_cast<scalar>(evalResults.pathList[i].back().simTimeMsec) / static_cast<scalar>(1000.0f);
      primaryDuration = std::max(primaryDuration, pathDuration);
    }
  }
  
  // Use real-time playback - elapsed seconds = simulation seconds
  scalar currentSimTime = elapsed;
  
  // Update stopwatch
  if (stopwatchVisible) {
    updateStopwatch(currentSimTime);
  }
  
  // Update control HUD alongside stopwatch
  if (controlsVisible) {
    updateControlsOverlay(currentSimTime);
  }
  
  // Check if animation is complete
  if (primaryDuration > 0.0 && currentSimTime >= primaryDuration) {
    currentSimTime = primaryDuration;
    isPlaybackActive = false;
    animationTimerId = 0;
    // Keep stopwatch visible showing final time
    if (stopwatchVisible) {
      updateStopwatch(currentSimTime);
    }
    std::cout << "Playback animation completed (duration: " << primaryDuration << "s)" << std::endl;
  }
  
  // Clear existing data
  this->paths->RemoveAllInputs();
  this->actuals->RemoveAllInputs();
  this->segmentGaps->RemoveAllInputs();
  this->blackboxTapes->RemoveAllInputs();
  this->blackboxHighlightTapes->RemoveAllInputs();
  this->xiaoVecArrows->RemoveAllInputs();

  // Always ensure these have at least empty data to prevent VTK pipeline errors
  vtkNew<vtkPolyData> emptyHighlightData;
  vtkNew<vtkPoints> emptyHighlightPoints;
  emptyHighlightData->SetPoints(emptyHighlightPoints);
  this->blackboxHighlightTapes->AddInputData(emptyHighlightData);

  vtkNew<vtkPolyData> emptyVecData;
  vtkNew<vtkPoints> emptyVecPoints;
  emptyVecData->SetPoints(emptyVecPoints);
  this->xiaoVecArrows->AddInputData(emptyVecData);

  
  // Render with synchronized time progress
  int renderArenas = (!blackboxAircraftStates.empty()) ? 1 : static_cast<int>(evalResults.pathList.size());
  for (int i = 0; i < renderArenas; i++) {
    vec3 offset = renderingOffset(i);
    
    std::vector<vec3> p = pathToVector(evalResults.pathList[i]);
    std::vector<vec3> a = stateToVector(evalResults.aircraftStateList[i]);
    
    // Time-based filtering: show data up to currentSimTime
    
    // Filter path points by simulation timestamp
    std::vector<vec3> visiblePathVector;
    if (!evalResults.pathList[i].empty()) {
      std::vector<Path> visiblePath;
      for (const auto& pathPoint : evalResults.pathList[i]) {
        scalar pathPointTime = static_cast<scalar>(pathPoint.simTimeMsec) / static_cast<scalar>(1000.0f); // Convert to seconds
        if (pathPointTime <= currentSimTime) {
          visiblePath.push_back(pathPoint);
        }
      }
      if (!visiblePath.empty()) {
        visiblePathVector = pathToVector(visiblePath);
      }
    }
    // Always add path data (empty if no visible path) to prevent VTK warnings
    this->paths->AddInputData(createPointSet(offset, visiblePathVector));
    
    // Filter aircraft states by timestamp
    std::vector<AircraftState> visibleStates;
    std::vector<vec3> visibleStateVector;
    if (!evalResults.aircraftStateList[i].empty()) {
      for (const auto& state : evalResults.aircraftStateList[i]) {
        scalar stateTime = static_cast<scalar>(state.getSimTimeMsec()) / static_cast<scalar>(1000.0f);
        if (stateTime <= currentSimTime) {
          visibleStates.push_back(state);
        }
      }
      if (!visibleStates.empty()) {
        visibleStateVector = stateToVector(visibleStates);
      }
    }
    // Always add aircraft data (empty if no visible states) to prevent VTK warnings
    this->actuals->AddInputData(createTapeSet(offset, visibleStateVector, 
      visibleStates.empty() ? std::vector<vec3>() : stateToOrientation(visibleStates)));
    
    // Add segment gaps only if we have visible states and the full path exists
    if (!visibleStates.empty() && !evalResults.pathList[i].empty()) {
      // Use full path for segment connections (aircraft states reference original path indices)
      std::vector<vec3> fullPathVector = pathToVector(evalResults.pathList[i]);
      
      // Further filter visible states to only include those with valid path references
      std::vector<AircraftState> validStates;
      for (const auto& state : visibleStates) {
        if (state.getThisPathIndex() < fullPathVector.size()) {
          validStates.push_back(state);
        }
      }
      
      if (!validStates.empty()) {
        this->segmentGaps->AddInputData(createSegmentSet(offset, validStates, fullPathVector));
      } else {
        // Add empty segment gaps to prevent VTK warnings
        std::vector<AircraftState> emptyStates;
        std::vector<vec3> emptyPath;
        this->segmentGaps->AddInputData(createSegmentSet(offset, emptyStates, emptyPath));
      }
    } else {
      // Add empty segment gaps to prevent VTK warnings
      std::vector<AircraftState> emptyStates;
      std::vector<vec3> emptyPath;
      this->segmentGaps->AddInputData(createSegmentSet(offset, emptyStates, emptyPath));
    }
    
    // Add blackbox data to first arena only (with animation)
    if (i == 0 && !blackboxAircraftStates.empty()) {
      vec3 blackboxOffset = offset;
      std::vector<vec3> a_bb = stateToVector(blackboxAircraftStates);

      // Calculate blackbox progress to sync with path timing (used for all blackbox/xiao animation)
      scalar blackboxProgress = static_cast<scalar>(1.0f); // Default to show all
      if (!blackboxAircraftStates.empty()) {
        // The blackbox test data starts at some offset but should sync with path time 0
        // Simply use the blackbox duration and sync with currentSimTime
        scalar blackboxStartTime = static_cast<scalar>(blackboxAircraftStates.front().getSimTimeMsec());
        scalar blackboxEndTime = static_cast<scalar>(blackboxAircraftStates.back().getSimTimeMsec());
        scalar blackboxDuration = (blackboxEndTime - blackboxStartTime) / static_cast<scalar>(1000000.0f); // Convert microseconds to seconds

        if (blackboxDuration > static_cast<scalar>(0.0f)) {
          // Blackbox animates synchronized with path: both start at currentSimTime=0
          blackboxProgress = std::min(static_cast<scalar>(1.0f), currentSimTime / blackboxDuration);
        }
      }

      if (a_bb.size() >= 2) {

        if ((inDecodeMode || inXiaoMode) && !testSpans.empty() && !showingFullFlight) {
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates), blackboxProgress));
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        } else if ((inDecodeMode || inXiaoMode) && showingFullFlight && !testSpans.empty()) {
          // For full flight mode, filter full flight blackbox data by time
          scalar fullStartTime = static_cast<scalar>(fullBlackboxAircraftStates.front().getSimTimeMsec()) / static_cast<scalar>(1000000.0f); // microseconds to seconds
          std::vector<AircraftState> visibleFullStates;
          
          for (const auto& state : fullBlackboxAircraftStates) {
            scalar stateTime = static_cast<scalar>(state.getSimTimeMsec()) / static_cast<scalar>(1000000.0f); // microseconds to seconds
            scalar relativeTime = stateTime - fullStartTime; // normalize to start at 0
            if (relativeTime <= currentSimTime) {
              visibleFullStates.push_back(state);
            }
          }
          
          if (visibleFullStates.size() >= 2) {
            std::vector<vec3> visibleFullVector = stateToVector(visibleFullStates);
            this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, visibleFullVector, stateToOrientation(visibleFullStates)));
            blackboxActor->GetProperty()->SetOpacity(0.25);
            blackboxActor->GetBackfaceProperty()->SetOpacity(0.25);
            
            // Animate highlighted test spans with time-based filtering
            for (const TestSpan& span : testSpans) {
              size_t startIdx = std::min(span.startIndex, fullBlackboxAircraftStates.size());
              size_t endIdx = std::min(span.endIndex + 1, fullBlackboxAircraftStates.size());
              
              if (startIdx < endIdx && startIdx < fullBlackboxAircraftStates.size()) {
                std::vector<AircraftState> spanStates;
                for (size_t j = startIdx; j < endIdx; j++) {
                  spanStates.push_back(fullBlackboxAircraftStates[j]);
                }
                
                // Filter span states by time
                std::vector<AircraftState> visibleSpanStates;
                for (const auto& state : spanStates) {
                  scalar stateTime = static_cast<scalar>(state.getSimTimeMsec()) / static_cast<scalar>(1000000.0f); // microseconds to seconds
                  scalar relativeTime = stateTime - fullStartTime; // use full flight start time
                  if (relativeTime <= currentSimTime) {
                    visibleSpanStates.push_back(state);
                  }
                }
                
                if (visibleSpanStates.size() >= 2) {
                  std::vector<vec3> visibleSpanVector = stateToVector(visibleSpanStates);
                  this->blackboxHighlightTapes->AddInputData(createTapeSet(blackboxOffset, visibleSpanVector, stateToOrientation(visibleSpanStates)));
                }
              }
            }
          }
        } else {
          // Standard blackbox mode with time-based filtering
          std::vector<AircraftState> visibleBlackboxStates;
          scalar blackboxStartTime = static_cast<scalar>(blackboxAircraftStates.front().getSimTimeMsec()) / static_cast<scalar>(1000000.0f); // microseconds to seconds
          
          for (const auto& state : blackboxAircraftStates) {
            scalar stateTime = static_cast<scalar>(state.getSimTimeMsec()) / static_cast<scalar>(1000000.0f); // microseconds to seconds
            scalar relativeTime = stateTime - blackboxStartTime; // normalize to start at 0
            if (relativeTime <= currentSimTime) {
              visibleBlackboxStates.push_back(state);
            }
          }
          
          if (visibleBlackboxStates.size() >= 2) {
            std::vector<vec3> visibleBlackboxVector = stateToVector(visibleBlackboxStates);
            this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, visibleBlackboxVector, stateToOrientation(visibleBlackboxStates)));
            blackboxActor->GetProperty()->SetOpacity(1.0);
            blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
          }
        }
      }

      // Render vec arrows for xiao mode (with time-based animation)
      if (inXiaoMode && !craftToTargetVectors.empty() && !blackboxAircraftStates.empty()) {
        // Filter vecs by timestamp, same as blackbox states and rabbits
        scalar vecStartTime = static_cast<scalar>(blackboxAircraftStates.front().getSimTimeMsec()) / static_cast<scalar>(1000000.0f);

        for (const auto& timestampedVec : craftToTargetVectors) {
          scalar vecTime = static_cast<scalar>(timestampedVec.timestampUs) / static_cast<scalar>(1000000.0f);
          scalar relativeTime = vecTime - vecStartTime;

          // Only show vecs that have occurred by currentSimTime
          if (relativeTime <= currentSimTime) {
            vec3 vecDir = timestampedVec.vector;
            if (vecDir.norm() > 0.01f) {
              // Find closest aircraft state by time to get position
              scalar closestTimeDiff = std::numeric_limits<scalar>::max();
              vec3 craftPos = blackboxAircraftStates[0].getPosition();

              for (const auto& state : blackboxAircraftStates) {
                scalar stateTime = static_cast<scalar>(state.getSimTimeMsec()) / static_cast<scalar>(1000000.0f);
                scalar stateRelativeTime = stateTime - vecStartTime;
                scalar timeDiff = std::abs(stateRelativeTime - relativeTime);
                if (timeDiff < closestTimeDiff) {
                  closestTimeDiff = timeDiff;
                  craftPos = state.getPosition();
                }
              }

              craftPos = craftPos + blackboxOffset;
              scalar vecLength = vecDir.norm();

              vtkNew<vtkArrowSource> arrowSource;
              arrowSource->SetTipResolution(8);
              arrowSource->SetShaftResolution(8);
              arrowSource->SetShaftRadius(0.1);    // Match blue line thickness (0.1 radius)
              arrowSource->SetTipRadius(0.2);      // Proportional tip thickness
              arrowSource->SetTipLength(0.12);     // Stubbier tip
              arrowSource->Update();

              vtkNew<vtkTransform> transform;
              transform->Translate(craftPos[0], craftPos[1], craftPos[2]);

              vec3 vecNorm = vecDir / vecLength;
              vec3 defaultDir(1.0, 0.0, 0.0);
              vec3 rotAxis = defaultDir.cross(vecNorm);
              scalar rotAxisLen = rotAxis.norm();
              if (rotAxisLen > 0.001f) {
                rotAxis = rotAxis / rotAxisLen;
                scalar dotProd = std::max(-1.0f, std::min(1.0f, defaultDir.dot(vecNorm)));
                scalar angle = std::acos(dotProd) * 180.0f / M_PI;
                transform->RotateWXYZ(angle, rotAxis[0], rotAxis[1], rotAxis[2]);
              }
              // Scale only in X direction (length) to keep shaft/tip thickness absolute
              transform->Scale(vecLength, 1.0, 1.0);

              vtkNew<vtkTransformPolyDataFilter> transformFilter;
              transformFilter->SetInputConnection(arrowSource->GetOutputPort());
              transformFilter->SetTransform(transform);
              transformFilter->Update();

              this->xiaoVecArrows->AddInputData(transformFilter->GetOutput());
            }
          }
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

  // Update xiao-specific actors
  if (inXiaoMode) {
    if (this->xiaoVecArrows->GetNumberOfInputConnections(0) > 0) {
      this->xiaoVecArrows->Update();
    }
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
  this->xiaoVecArrows->RemoveAllInputs();

  // Always ensure these have at least empty data to prevent VTK pipeline errors
  vtkNew<vtkPolyData> emptyHighlightData;
  vtkNew<vtkPoints> emptyHighlightPoints;
  emptyHighlightData->SetPoints(emptyHighlightPoints);
  this->blackboxHighlightTapes->AddInputData(emptyHighlightData);

  vtkNew<vtkPolyData> emptyVecData;
  vtkNew<vtkPoints> emptyVecPoints;
  emptyVecData->SetPoints(emptyVecPoints);
  this->xiaoVecArrows->AddInputData(emptyVecData);


  // Render with full progress (1.0) using existing data
  int maxArenas = evalResults.pathList.size();
  for (int i = 0; i < maxArenas; i++) {
    vec3 offset = renderingOffset(i);
    
    std::vector<vec3> p = pathToVector(evalResults.pathList[i]);
    std::vector<vec3> a = stateToVector(evalResults.aircraftStateList[i]);
    
    if (!p.empty()) {
      this->paths->AddInputData(createPointSet(offset, p)); // Full progress (no timeProgress param)
    }
    if (!a.empty()) {
      this->actuals->AddInputData(createTapeSet(offset, a, stateToOrientation(evalResults.aircraftStateList[i]))); // Full progress
    }
    if (!a.empty() && !p.empty()) {
      this->segmentGaps->AddInputData(createSegmentSet(offset, evalResults.aircraftStateList[i], p)); // Full progress
    }
    
    // Add blackbox/xiao data to first arena only (full progress)
    if (i == 0 && !blackboxAircraftStates.empty()) {
      vec3 blackboxOffset = offset;
      std::vector<vec3> a_bb = stateToVector(blackboxAircraftStates);

      if (a_bb.size() >= 2) {
        if ((inDecodeMode || inXiaoMode) && !testSpans.empty() && !showingFullFlight) {
          // Single test span mode - full brightness
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates)));
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        } else if ((inDecodeMode || inXiaoMode) && showingFullFlight && !testSpans.empty()) {
          // All flight mode - dimmed background with highlighted test spans
          createHighlightedFlightTapes(blackboxOffset);
        } else {
          // No test spans or other modes - normal rendering
          this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates)));
          blackboxActor->GetProperty()->SetOpacity(1.0);
          blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
        }
      }

      // Render vec arrows for xiao mode
      if (inXiaoMode && !craftToTargetVectors.empty() && !blackboxAircraftStates.empty()) {
        // Create arrows from craft positions pointing in vec direction
        size_t arrowCount = 0;

        for (const auto& timestampedVec : craftToTargetVectors) {
          vec3 vecDir = timestampedVec.vector;
          // Only render if vec is non-zero
          if (vecDir.norm() > 0.01f) {
            // Find closest aircraft state by timestamp to get position
            scalar vecTime = static_cast<scalar>(timestampedVec.timestampUs);
            scalar closestTimeDiff = std::numeric_limits<scalar>::max();
            vec3 craftPos = blackboxAircraftStates[0].getPosition();

            for (const auto& state : blackboxAircraftStates) {
              scalar stateTime = static_cast<scalar>(state.getSimTimeMsec());
              scalar timeDiff = std::abs(stateTime - vecTime);
              if (timeDiff < closestTimeDiff) {
                closestTimeDiff = timeDiff;
                craftPos = state.getPosition();
              }
            }

            craftPos = craftPos + blackboxOffset;
            scalar vecLength = vecDir.norm();

            // Create arrow pointing along vecDir
            vtkNew<vtkArrowSource> arrowSource;
            arrowSource->SetTipResolution(8);
            arrowSource->SetShaftResolution(8);
            arrowSource->SetShaftRadius(0.1);    // Match blue line thickness (0.1 radius)
            arrowSource->SetTipRadius(0.2);      // Proportional tip thickness
            arrowSource->SetTipLength(0.12);     // Stubbier tip
            arrowSource->Update();

            // Transform arrow to point from craft to target
            vtkNew<vtkTransform> transform;
            transform->Translate(craftPos[0], craftPos[1], craftPos[2]);

            // Calculate rotation to align arrow with vec direction
            vec3 vecNorm = vecDir / vecLength;
            vec3 defaultDir(1.0, 0.0, 0.0);  // vtkArrowSource points along +X by default

            // Compute rotation axis and angle
            vec3 rotAxis = defaultDir.cross(vecNorm);
            scalar rotAxisLen = rotAxis.norm();
            if (rotAxisLen > 0.001f) {
              rotAxis = rotAxis / rotAxisLen;
              scalar dotProd = std::max(-1.0f, std::min(1.0f, defaultDir.dot(vecNorm)));
              scalar angle = std::acos(dotProd) * 180.0f / M_PI;
              transform->RotateWXYZ(angle, rotAxis[0], rotAxis[1], rotAxis[2]);
            }

            // Scale only in X direction (length) to keep shaft/tip thickness absolute
            transform->Scale(vecLength, 1.0, 1.0);

            vtkNew<vtkTransformPolyDataFilter> transformFilter;
            transformFilter->SetInputConnection(arrowSource->GetOutputPort());
            transformFilter->SetTransform(transform);
            transformFilter->Update();

            this->xiaoVecArrows->AddInputData(transformFilter->GetOutput());
            arrowCount++;
          }
        }
        if (arrowCount > 0) {
          std::cout << "Rendered " << arrowCount << " vec arrows" << std::endl;
        }
      }

    }
  }

  // Handle xiao mode when there are no arenas (maxArenas == 0)
  if (inXiaoMode && maxArenas == 0 && !blackboxAircraftStates.empty()) {
    vec3 blackboxOffset(0.0f, 0.0f, 0.0f);  // No arena offset
    std::vector<vec3> a_bb = stateToVector(blackboxAircraftStates);

    // Render blackbox tape
    if (a_bb.size() >= 2) {
      this->blackboxTapes->AddInputData(createTapeSet(blackboxOffset, a_bb, stateToOrientation(blackboxAircraftStates)));
      blackboxActor->GetProperty()->SetOpacity(1.0);
      blackboxActor->GetBackfaceProperty()->SetOpacity(1.0);
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
  if (inXiaoMode && this->xiaoVecArrows->GetNumberOfInputConnections(0) > 0) {
    this->xiaoVecArrows->Update();
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
  std::cout << "  -x, --xiaofile FILE      Specify xiao log file to overlay\n";
  std::cout << "                           Mutually exclusive with -d option\n";
  std::cout << "  -i, --config FILE        Use specified config file (default: autoc.ini)\n";
  std::cout << "  -h, --help               Show this help message\n";
  std::cout << "\n";
  std::cout << "Examples:\n";
  std::cout << "  " << progName << "                                    # Render all arenas without blackbox data\n";
  std::cout << "  " << progName << " -d 'blackbox_decode file.bbl'     # Pipe through decoder command\n";
  std::cout << "  " << progName << " -d -                              # Read CSV from stdin\n";
  std::cout << "  cat data.csv | " << progName << " -d -               # Read CSV from stdin (alternative)\n";
  std::cout << "  " << progName << " -x flight.txt                     # Overlay xiao log data\n";
}

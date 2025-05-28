#include <regex>

#include "renderer.h"

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/ListObjectsV2Request.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/client/ClientConfiguration.h>

EvalResults evalResults;
std::string computedKeyName = "";
Renderer renderer;

std::shared_ptr<Aws::S3::S3Client> getS3Client() {
  // real S3 or local minio?
  Aws::Client::ClientConfiguration clientConfig;
  Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::RequestDependent;
  if (strcmp("default", "minio" /*extraCfg.s3Profile*/) != 0) {
    clientConfig.endpointOverride = "http://localhost:9000"; // MinIO server address
    clientConfig.scheme = Aws::Http::Scheme::HTTP; // Use HTTP instead of HTTPS
    clientConfig.verifySSL = false; // Disable SSL verification for local testing

    policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Never;
  }

  auto credentialsProvider = Aws::MakeShared<Aws::Auth::ProfileConfigFileAWSCredentialsProvider>("CredentialsProvider", "minio" /*extraCfg.s3Profile*/);

  return Aws::MakeShared<Aws::S3::S3Client>("S3Client",
    credentialsProvider,
    clientConfig,
    Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Always,
    false
  );
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
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  for (const auto& point : points) {
    Eigen::Vector3d rPoint = point + offset;
    vtp->InsertNextPoint(rPoint[0], rPoint[1], rPoint[2]);
  }

  vtkSmartPointer<vtkPolyLine> lines = vtkSmartPointer<vtkPolyLine>::New();
  lines->GetPointIds()->SetNumberOfIds(points.size());
  for (int i = 0; i < points.size(); ++i) {
    lines->GetPointIds()->SetId(i, i);
  }

  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  cells->InsertNextCell(lines);

  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(vtp);
  polyData->SetLines(cells);

  // vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  // glyphFilter->SetInputData(polyData);
  // glyphFilter->Update();

  // vtkSmartPointer<vtkPolyData> pointPolyData = vtkSmartPointer<vtkPolyData>::New();
  // pointPolyData->ShallowCopy(glyphFilter->GetOutput());

  return polyData;
}

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(Eigen::Vector3d offset, const std::vector<AircraftState> state, const std::vector<Eigen::Vector3d> end) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < state.size(); i++) {
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

  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);
  polyData->SetLines(lines);

  return polyData;
}

/*
 ** actual data is rendered as a tape with a top and bottom
 */
vtkSmartPointer<vtkPolyData> Renderer::createTapeSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> points,
  const std::vector<Eigen::Vector3d> normals) {
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkDoubleArray> orientations = vtkSmartPointer<vtkDoubleArray>::New();
  orientations->SetNumberOfComponents(3);
  orientations->SetName("Orientations");

  vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
  polyLine->GetPointIds()->SetNumberOfIds(points.size());

  for (size_t i = 0; i < points.size(); ++i) {
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

  // Create a ribbon filter
  vtkSmartPointer<vtkRibbonFilter> ribbonFilter = vtkSmartPointer<vtkRibbonFilter>::New();
  ribbonFilter->SetInputData(polyData);
  ribbonFilter->SetWidth(0.5);
  ribbonFilter->Update();

  return ribbonFilter->GetOutput();
}

/*
 * re-render the generation number asked for
 */
bool Renderer::updateGenerationDisplay(int newGen) {
  // now do initial fetch
  Aws::S3::Model::GetObjectRequest request;
  request.SetBucket("autoc-storage"); // TODO extraCfg
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

  for (int i = 0; i < evalResults.pathList.size(); i++) {
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
  }

  this->planeData->Update();
  this->paths->Update();
  this->actuals->Update();
  this->segmentGaps->Update();

  // Update the window title
  std::string title = keyName + " - " + std::to_string(10000 - newGen);
  renderWindow->SetWindowName(title.c_str());

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


void Renderer::initialize() {
  // Create a renderer and render window interactor
  vtkNew<vtkRenderer> renderer;
  renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(1080, 900);

  renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

  // Set the interactor style to extended trackball camera
  vtkNew<CustomInteractorStyle> interactorStyle;
  renderWindowInteractor->SetInteractorStyle(interactorStyle);
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Create command objects for the callbacks
  vtkNew<NextModelCommand> nextModelCommand;
  nextModelCommand->SetRenderer(this);
  vtkNew<PreviousModelCommand> previousModelCommand;
  previousModelCommand->SetRenderer(this);

  // Add observers for the custom events
  interactorStyle->AddObserver(NextModelEvent, nextModelCommand);
  interactorStyle->AddObserver(PreviousModelEvent, previousModelCommand);

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

  // Temporary until first update
  vtkNew<vtkPolyData> emptyPolyData;
  vtkNew<vtkPoints> emptyPoints;
  emptyPolyData->SetPoints(emptyPoints);

  // Add the empty polydata to elements
  paths->AddInputData(emptyPolyData);
  actuals->AddInputData(emptyPolyData);
  segmentGaps->AddInputData(emptyPolyData);
  planeData->AddInputData(emptyPolyData);

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
  property->SetColor(1.0, 0.5, 0.0);  // Orange for front face (top)
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
    points.push_back(s.getOrientation() * -Eigen::Vector3d::UnitZ());
  }
  return points;
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
  // which build are we going to look at
  std::string keyName = "";

  if (argc > 1) {
    computedKeyName = argv[1];
  }
  else {
    // later, we'll search for latest build
  }
  // AWS setup
  Aws::SDKOptions options;
  Aws::InitAPI(options);

  std::shared_ptr<Aws::S3::S3Client> s3_client = getS3Client();

  // should we look up the latest run?
  if (computedKeyName.empty()) {
    Aws::S3::Model::ListObjectsV2Request listFolders;
    listFolders.SetBucket("autoc-storage"); // TODO extraCfg
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
  listItem.SetBucket("autoc-storage"); // TODO extraCfg
  listItem.SetPrefix(computedKeyName + "gen");
  bool isTruncated = false;
  do {
    auto outcome = s3_client->ListObjectsV2(listItem);
    if (outcome.IsSuccess()) {
      // Objects are already in reverse lexicographical order
      for (const auto& object : outcome.GetResult().GetContents()) {
        keyName = object.GetKey();
        break;
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

  // kick it off
  renderer.renderWindowInteractor->Start();

  return 0;
}

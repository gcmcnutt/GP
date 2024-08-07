#include "renderer.h"

void PrintPolyDataInfo(vtkPolyData* polyData)
{
    vtkPoints* points = polyData->GetPoints();
    vtkCellArray* cells = polyData->GetPolys();
    vtkIdList* idList = vtkIdList::New();

    if (points == NULL) {
        printf("No points in the polydata\n");
    } else {
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
  int sideLength = std::ceil(std::sqrt(extraCfg.simNumPathsPerGen));
  
  int row = i / sideLength;
  int col = i % sideLength;

  // for now put them in a line
  double xOffset = col * (FIELD_SIZE + FIELD_GAP);
  double yOffset = row * (FIELD_SIZE + FIELD_GAP);

  return Eigen::Vector3d(xOffset, yOffset, 0.0);
}

// VTK timer event callback
void Renderer::Execute(vtkObject* caller, unsigned long eventId, void* vtkNotUsed(callData)) {
  std::lock_guard<std::recursive_mutex> lock(dataMutex);
  if (vtkCommand::TimerEvent == eventId) {
    if (!newDataAvailable) {
      return;
    }

    // Update mappers
    vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkTubeFilter> tubeFilter1 = vtkSmartPointer<vtkTubeFilter>::New();
    tubeFilter1->SetInputConnection(paths->GetOutputPort());
    tubeFilter1->SetRadius(0.2);
    tubeFilter1->SetNumberOfSides(10);
    mapper1->SetInputConnection(tubeFilter1->GetOutputPort());

    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkTubeFilter> tubeFilter2 = vtkSmartPointer<vtkTubeFilter>::New();
    tubeFilter2->SetInputConnection(actuals->GetOutputPort());
    tubeFilter2->SetRadius(0.3);
    tubeFilter2->SetNumberOfSides(10);
    mapper2->SetInputConnection(tubeFilter2->GetOutputPort());

    vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkTubeFilter> tubeFilter3 = vtkSmartPointer<vtkTubeFilter>::New();
    tubeFilter3->SetInputConnection(segmentGaps->GetOutputPort());
    tubeFilter3->SetRadius(0.1);
    tubeFilter3->SetNumberOfSides(10);
    mapper3->SetInputConnection(tubeFilter3->GetOutputPort());

    // Update actors
    actor1->SetMapper(mapper1);
    actor2->SetMapper(mapper2);
    actor3->SetMapper(mapper3);

    // render
    vtkRenderWindowInteractor* iren = static_cast<vtkRenderWindowInteractor*>(caller);
    iren->GetRenderWindow()->Render();

    newDataAvailable = false;
  }
  else if (vtkCommand::ExitEvent == eventId) {
    exitFlag = true;
  }
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

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(Eigen::Vector3d offset, const std::vector<Eigen::Vector3d> start, const std::vector<Eigen::Vector3d> end) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < start.size(); i++) {
    Eigen::Vector3d rStart = start[i] + offset;
    Eigen::Vector3d rEnd = end[i] + offset;
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

void Renderer::addPathElementList(std::vector<Path> plan, std::vector<Path> actual) {
  std::lock_guard<std::recursive_mutex> lock(dataMutex);
  actualList.push_back(actual);
  pathList.push_back(plan);
}

void Renderer::RenderInBackground(vtkSmartPointer<vtkRenderWindow> renderWindow) {
  // Create a renderer and render window interactor
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(renderer);
  renderWindow->SetSize(1080, 900);

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();

    // Set the interactor style to trackball camera
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Configure the camera
  vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetPosition(-50, 0, -2);        // behind the action
  camera->SetFocalPoint(0, 0, -10);       // TODO center this on arenas
  camera->SetViewUp(0, 0, -1);           // Set the view up vector
  camera->SetViewAngle(60);             // Set the field of view (FOV) in degrees

  // Apply the camera settings to the renderer
  renderer->SetActiveCamera(camera);
  renderer->ResetCameraClippingRange(); // Adjust clipping range based on the scene

  // Create the axis actor
  vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
  
  // Create the orientation marker widget
  vtkSmartPointer<vtkOrientationMarkerWidget> orientationMarker = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
  orientationMarker->SetOrientationMarker(axes);
  orientationMarker->SetInteractor(renderWindowInteractor);
  orientationMarker->SetViewport(0.0, 0.0, 0.2, 0.2); // Position and size in the window
  orientationMarker->SetEnabled(1);
  orientationMarker->InteractiveOn();

  vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();

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

  // create the checkerboards
  planeData = vtkSmartPointer<vtkAppendPolyData>::New();
  for (int j = 0; j < extraCfg.simNumPathsPerGen; j++) {
    Eigen::Vector3d offset = renderingOffset(j);

    // Create a plane source at z = 0
    vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();

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
    vtkSmartPointer<vtkUnsignedCharArray> cellData = vtkSmartPointer<vtkUnsignedCharArray>::New();
    cellData->SetNumberOfComponents(4);
    cellData->SetNumberOfTuples(planeSource->GetOutput()->GetNumberOfCells());

    // checkerboard
    for (int i = 0; i < planeSource->GetOutput()->GetNumberOfCells(); i++) {
      if (i % 2 ^ (i / 10) % 2) {
        double rgb[4] = {255.0, 255.0, 255.0, 100.0};
        cellData->InsertTuple(i, rgb);
      } else {      
        double rgb[4] = {0.0, 0.0, 0.0, 100.0};
        cellData->InsertTuple(i, rgb);
      }
    }
    planeSource->GetOutput()->GetCellData()->SetScalars(cellData);
    planeSource->Update();
    planeData->AddInputData(planeSource->GetOutput());
  }

  // Enable anti-aliasing (multi-sampling)
  renderWindow->SetMultiSamples(4); // Use 4x MSAA

  // Enable depth peeling for proper transparency rendering
  renderWindow->SetAlphaBitPlanes(1);
  renderWindow->SetMultiSamples(0);
  renderer->SetUseDepthPeeling(1);
  renderer->SetMaximumNumberOfPeels(100);  // Maximum number of depth peels
  renderer->SetOcclusionRatio(0.1);        // Occlusion ratio

  // Create an actor for the plane
  planeActor = vtkSmartPointer<vtkActor>::New();vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  planeMapper->SetInputConnection(planeData->GetOutputPort());
  planeActor->SetMapper(planeMapper);

  renderer->AddActor(planeActor);
  renderer->AddActor(actor1);
  renderer->AddActor(actor2);
  renderer->AddActor(actor3);

  // Add the timer callback and close callback
  renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, this);
  renderWindowInteractor->AddObserver(vtkCommand::ExitEvent, this);
  renderWindowInteractor->CreateRepeatingTimer(100); // 100 ms interval

  // Initialize the rendering and interaction
  renderWindow->Render();
  renderWindowInteractor->Start();
}

std::vector<Eigen::Vector3d> Renderer::pathToVector(std::vector<Path> path) {
  std::vector<Eigen::Vector3d> points;
  for (const auto& p : path) {
    points.push_back(p.start);
  }
  return points;
}

void Renderer::update() {
  std::lock_guard<std::recursive_mutex> lock(dataMutex);
  this->paths = vtkSmartPointer<vtkAppendPolyData>::New();
  this->actuals = vtkSmartPointer<vtkAppendPolyData>::New();
  this->segmentGaps = vtkSmartPointer<vtkAppendPolyData>::New();

  for (int i = 0; i < actualList.size(); i++) {
    Eigen::Vector3d offset = renderingOffset(i);

    // paths
    std::vector<Eigen::Vector3d> p = pathToVector(pathList[i]);
    this->paths->AddInputData(createPointSet(offset, p));
    this->paths->Update();

    // actuals
    std::vector<Eigen::Vector3d> a = pathToVector(actualList[i]);
    this->actuals->AddInputData(createPointSet(offset, a));
    this->actuals->Update();

    // segment gaps
    this->segmentGaps->AddInputData(createSegmentSet(offset, a, p));
    this->segmentGaps->Update();
  }
  actualList.clear();
  pathList.clear();

  newDataAvailable = true;
}

void Renderer::start() {
  // Create a VTK render window
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

  std::thread renderThread([this, renderWindow]() { RenderInBackground(renderWindow); });
  renderThread.detach(); // Detach the thread to run independently
}

bool Renderer::isRunning() {
  std::lock_guard<std::recursive_mutex> lock(dataMutex);
  return !exitFlag;
}
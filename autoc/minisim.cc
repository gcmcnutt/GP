/* test sim for aircraft */
#include <iostream>
#include <cmath>
#include <thread>
#include <mutex>

#include "minisim.h"
#include "pathgen.h"

using namespace std;

AircraftState::AircraftState() {
      this->dRelVel = 0;
      this->dPhi = 0;
      this->dTheta = 0;
      this->dPsi = 0;
      this->X = 0;
      this->Y = 0;
      this->Z = 0;
      this->R_X = 0;
      this->R_Y = 0;
      this->R_Z = 0;
} 

AircraftState::AircraftState(double dRelVel, double dPhi, double dTheta, double dPsi, double X, double Y, double Z, double R_X, double R_Y, double R_Z) {
      this->dRelVel = dRelVel;
      this->dPhi = dPhi;
      this->dTheta = dTheta;
      this->dPsi = dPsi;
      this->X = X;
      this->Y = Y;
      this->Z = Z;
      this->R_X = R_X;
      this->R_Y = R_Y;
      this->R_Z = R_Z;
}


Aircraft::Aircraft(AircraftState *state) {
      this->state = state;
      this->pitchCommand = 0;
      this->rollCommand = 0;
      this->throttleCommand = 0;
}

void Aircraft::setState(AircraftState *state) {
      this->state = state;
}

AircraftState *Aircraft::getState() {
      return state;
}

double Aircraft::setPitchCommand(double pitchCommand) {
  this->pitchCommand = pitchCommand;
  return pitchCommand;
}

double Aircraft::getPitchCommand() {
  return pitchCommand;
}

double Aircraft::setRollCommand(double rollCommand) {
  this->rollCommand = rollCommand;
  return rollCommand;
}

double Aircraft::getRollCommand() {
  return rollCommand;
}

double Aircraft::setThrottleCommand(double throttleCommand) {
  this->throttleCommand = throttleCommand;
  return throttleCommand;
}

double Aircraft::getThrottleCommand() {
  return throttleCommand;
}

void Aircraft::advanceState(double dt) {
  // get velocity
  double dVel = state->dRelVel;

  // get current roll state
  double rollCurrent = state->dPhi;

  // get roll command: negative is roll left, positive is roll right
  double rollCommand = max(min(getRollCommand(), 1.0), -1.0); // XXX constant

  // compute new roll orientation
  double rollNew = remainder(rollCurrent + rollCommand * dt * MAX_ROLL_RATE_RADSEC, M_PI * 2);

  // compute ground left/right force from current roll
  double dForceLR = sin(rollNew) * dt * MAX_YAW_RATE_RADSEC;

  // get current heading
  double headingCurrent = state->dPsi;

  // update heading based on roll
  double newHeading = remainder(headingCurrent + dForceLR, M_PI * 2);

  // get position
  Point3D position = {state->X, state->Y, state->Z};

  // update XY position based on heading, velocity, and dt
  position.x += dVel * std::cos(newHeading) * dt; // XXX need to work on z+ is actually down
  position.y += dVel * std::sin(newHeading) * dt;

  // update state as a result
  state->X = position.x;
  state->Y = position.y;
  state->dPhi = rollNew;
  state->dPsi = newHeading;
}

void Aircraft::toString(char *output) {
  sprintf(output, "AircraftState: %f %f %f %f %f %f %f %f %f %f %f  Command: %f %f %f\n", state->dRelVel, state->dPhi, state->dTheta,
    state->dPsi, state->dPhi, state->X, state->Y, state->Z, state->R_X, state->R_Y, state->R_Z,
    pitchCommand, rollCommand, throttleCommand);
}




// VTK timer event callback
void Renderer::Execute(vtkObject* caller, unsigned long eventId, void* vtkNotUsed(callData)) {
  if (vtkCommand::TimerEvent == eventId) {
    std::lock_guard<std::mutex> lock(dataMutex);

    if (!newDataAvailable) {
      return;
    }

    // Update mappers
    vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper1->SetInputData(path);

    vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper2->SetInputData(actual);

    vtkSmartPointer<vtkPolyDataMapper> mapper3 = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper3->SetInputData(segmentGap);

    // Create a mapper for the plane
    vtkSmartPointer<vtkPolyDataMapper> planeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    planeMapper->SetInputConnection(planeSource->GetOutputPort());

    // Update actors
    planeActor->SetMapper(planeMapper);
    actor1->SetMapper(mapper1);
    actor2->SetMapper(mapper2);
    actor3->SetMapper(mapper3);

    // render
    vtkRenderWindowInteractor* iren = static_cast<vtkRenderWindowInteractor*>(caller);
    iren->GetRenderWindow()->Render();

    newDataAvailable = false;
  }
}


vtkSmartPointer<vtkPolyData> Renderer::createPointSet(const std::vector<Point3D> points) {
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  for (const auto& point : points) {
      vtp->InsertNextPoint(point.x, point.y, point.z);
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

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(const std::vector<Point3D> start, const std::vector<Point3D> end) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < start.size(); i++) {
    points->InsertNextPoint(start[i].x, start[i].y, start[i].z);
    points->InsertNextPoint(end[i].x, end[i].y, end[i].z);
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
  camera->SetPosition(0, -50, 2);        // behind the action
  camera->SetFocalPoint(0, 0, 10);       // Start of initial height
  camera->SetViewUp(0, 0, 1);           // Set the view up vector
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

  // Create a plane source at z = 0
  planeSource = vtkSmartPointer<vtkPlaneSource>::New();

  double width = 100.0;
  double height = 100.0;
  int resolution = 10;
  planeSource->SetOrigin(-width / 2.0, -height / 2.0, 0.0);
  planeSource->SetPoint1(width / 2.0, -height / 2.0, 0.0);
  planeSource->SetPoint2(-width / 2.0, height / 2.0, 0.0);
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

  // Enable anti-aliasing (multi-sampling)
  renderWindow->SetMultiSamples(4); // Use 4x MSAA

  // Enable depth peeling for proper transparency rendering
  renderer->SetUseDepthPeeling(1);
  renderer->SetMaximumNumberOfPeels(100);  // Maximum number of depth peels
  renderer->SetOcclusionRatio(0.1);        // Occlusion ratio

  // Create an actor for the plane
  planeActor = vtkSmartPointer<vtkActor>::New();

  renderer->AddActor(planeActor);
  renderer->AddActor(actor1);
  renderer->AddActor(actor2);
  renderer->AddActor(actor3);

  // Add the timer callback
  renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, this);
  renderWindowInteractor->CreateRepeatingTimer(100); // 100 ms interval

  // Initialize the rendering and interaction
  renderWindow->Render();
  renderWindowInteractor->Start();
}

void Renderer::update(std::vector<Point3D> path, std::vector<Point3D> actual) {
  std::lock_guard<std::mutex> lock(dataMutex);
  this->path = createPointSet(path);
  this->actual = createPointSet(actual);
  this->segmentGap = createSegmentSet(actual, path);

  newDataAvailable = true;
}

void Renderer::start() {
  // Create a VTK render window
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

  std::thread renderThread([this, renderWindow]() { RenderInBackground(renderWindow); });
  renderThread.detach(); // Detach the thread to run independently

  // TODO
  // renderThread.join();
}




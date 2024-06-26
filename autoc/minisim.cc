/* test sim for aircraft */
#include <iostream>
#include <cmath>
#include <thread>
#include <mutex>
#include <array>

#include "minisim.h"
#include "pathgen.h"

using namespace std;

Aircraft::Aircraft(double dRelVel, Eigen::Quaterniond aircraft_orientation, Eigen::Vector3d position, double R_X, double R_Y, double R_Z) {
      this->dRelVel = dRelVel;
      this->aircraft_orientation = aircraft_orientation;
      this->position = position;
      this->R_X = R_X;
      this->R_Y = R_Y;
      this->R_Z = R_Z;
      
      this->pitchCommand = 0;
      this->rollCommand = 0;
      this->throttleCommand = 0;
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
  // get current roll state, compute left/right force (positive roll is right)
  double rollCommand = std::clamp(getRollCommand(), -1.0, 1.0);
  double delta_roll = remainder(rollCommand * dt * MAX_ROLL_RATE_RADSEC, M_PI);

  // get current pitch state, compute up/down force (positive pitch is up)
  double pitchCommand = std::clamp(getPitchCommand(), -1.0, 1.0);
  double delta_pitch = remainder(pitchCommand * dt * MAX_PITCH_RATE_RADSEC, M_PI);

  // adjust velocity as a function of throttle (-1:1)
  double throttle = std::clamp(getThrottleCommand(), -1.0, +1.0);
  dRelVel = SIM_INITIAL_VELOCITY + (throttle * SIM_THROTTLE_SCALE);

  // Convert pitch and roll updates to quaternions (in the body frame)
  Eigen::Quaterniond delta_pitch_quat(Eigen::AngleAxisd(delta_pitch, Eigen::Vector3d::UnitY()));
  Eigen::Quaterniond delta_roll_quat(Eigen::AngleAxisd(delta_roll, Eigen::Vector3d::UnitX()));

  // Apply the pitch and roll adjustments to the aircraft's orientation
  aircraft_orientation = delta_pitch_quat * aircraft_orientation;
  aircraft_orientation = delta_roll_quat * aircraft_orientation;

  // Normalize the resulting quaternion
  aircraft_orientation.normalize();
 
  // Define the initial velocity vector in the body frame
  Eigen::Vector3d velocity_body(dRelVel * dt, 0, 0);

  // Rotate the velocity vector using the updated quaternion
  Eigen::Vector3d velocity_world = aircraft_orientation * velocity_body;

  // adjust position
  position += velocity_world;
}

void Aircraft::toString(char *output) {
  sprintf(output, "AircraftState: %f %f %f %f %f %f %f %f %f %f %f  Command: %f %f %f\n", dRelVel, 
    aircraft_orientation.w(), aircraft_orientation.x(), aircraft_orientation.y(), aircraft_orientation.z(),
    position[0], position[1], position[2], R_X, R_Y, R_Z,
    pitchCommand, rollCommand, throttleCommand);
}




// VTK timer event callback
void Renderer::Execute(vtkObject* caller, unsigned long eventId, void* vtkNotUsed(callData)) {
  std::lock_guard<std::mutex> lock(dataMutex);
  if (vtkCommand::TimerEvent == eventId) {

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
  else if (vtkCommand::ExitEvent == eventId) {
    exitFlag = true;
  }
}


vtkSmartPointer<vtkPolyData> Renderer::createPointSet(const std::vector<Eigen::Vector3d> points) {
  vtkSmartPointer<vtkPoints> vtp = vtkSmartPointer<vtkPoints>::New();
  for (const auto& point : points) {
      vtp->InsertNextPoint(point[0], point[1], point[2]);
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

vtkSmartPointer<vtkPolyData> Renderer::createSegmentSet(const std::vector<Eigen::Vector3d> start, const std::vector<Eigen::Vector3d> end) {
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < start.size(); i++) {
    points->InsertNextPoint(start[i][0], start[i][1], start[i][2]);
    points->InsertNextPoint(end[i][0], end[i][1], end[i][2]);
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
  camera->SetPosition(-50, 0, -2);        // behind the action
  camera->SetFocalPoint(0, 0, SIM_INITIAL_ALTITUDE);       // Start of initial height
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

  // Add the timer callback and close callback
  renderWindowInteractor->AddObserver(vtkCommand::TimerEvent, this);
  renderWindowInteractor->AddObserver(vtkCommand::ExitEvent, this);
  renderWindowInteractor->CreateRepeatingTimer(100); // 100 ms interval

  // Initialize the rendering and interaction
  renderWindow->Render();
  renderWindowInteractor->Start();
}

void Renderer::update(std::vector<Eigen::Vector3d> path, std::vector<Eigen::Vector3d> actual) {
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

bool Renderer::isRunning() {
  std::lock_guard<std::mutex> lock(dataMutex);
  return !exitFlag;
}

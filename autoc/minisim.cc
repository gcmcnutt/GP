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

  // get heading CW from North
  double heading = state->dPsi;

  // get roll command: negative is roll left, positive is roll right (-1:1)
  double roll = max(min(getRollCommand(), MAX_SERVO_DEFLECTION), -MAX_SERVO_DEFLECTION) / MAX_SERVO_DEFLECTION;

  // get position
  Point3D position = {state->X, state->Y, state->Z};

  // update heading based on roll
  heading = remainder(heading + roll * dt * MAX_DELTA_ANGLE_RADSEC, M_PI * 2);

  // update XY position based on heading, velocity, and dt
  position.x += dVel * std::cos(heading) * dt;
  position.y += dVel * std::sin(heading) * dt;

  // update state as a result
  state->X = position.x;
  state->Y = position.y;
  state->dPsi = heading;
}

void Aircraft::toString(char *output) {
  sprintf(output, "AircraftState: %f %f %f %f %f %f %f %f %f %f  Command: %f %f %f\n", state->dRelVel, state->dPhi, state->dTheta,
    state->dPsi, state->X, state->Y, state->Z, state->R_X, state->R_Y, state->R_Z,
    pitchCommand, rollCommand, throttleCommand);
}


vtkSmartPointer<vtkPolyData> createPointSet(const std::vector<Point3D> points) {
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

void Renderer::RenderInBackground(vtkSmartPointer<vtkRenderWindow> renderWindow)
{
    // Create a renderer and render window interactor
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(800, 600);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Initialize the rendering and interaction
    renderWindow->Render();
    renderWindowInteractor->Initialize();

    vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();


    // Rendering loop
    while (true)
    {
        std::unique_lock<std::mutex> lock(dataMutex);
        dataCondition.wait(lock, [this] { return newDataAvailable; });

        // Update the rendering with new data
        renderer->RemoveAllViewProps();

        // Create mappers
        vtkSmartPointer<vtkPolyDataMapper> mapper1 = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper1->SetInputData(path);

        vtkSmartPointer<vtkPolyDataMapper> mapper2 = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper2->SetInputData(actual);

        // Create actors
        vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
        actor1->SetMapper(mapper1);
        actor1->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
        actor1->GetProperty()->SetPointSize(4);

        vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
        actor2->SetMapper(mapper2);
        actor2->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
        actor2->GetProperty()->SetPointSize(4);

        // Create a renderer
        renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor1);
        renderer->AddActor(actor2);
        renderer->SetBackground(colors->GetColor3d("Black").GetData());

        // Create a render window
        // renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        // renderWindow->SetSize(800, 600);

        // Create an interactor
        renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        renderWindowInteractor->SetRenderWindow(renderWindow);

        // Start the interaction
        renderWindow->Render();

        newDataAvailable = false;
    }
}

void Renderer::update(std::vector<Point3D> path, std::vector<Point3D> actual) {
  {
    std::lock_guard<std::mutex> lock(dataMutex);
    this->path = createPointSet(path);
    this->actual = createPointSet(actual);
    newDataAvailable = true;
  }
  dataCondition.notify_one();
}

void Renderer::start() {
    // Create a VTK render window
    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

    std::thread renderThread([this, renderWindow]() { RenderInBackground(renderWindow); });
    renderThread.detach(); // Detach the thread to run independently

    // TODO
    // renderThread.join();
}




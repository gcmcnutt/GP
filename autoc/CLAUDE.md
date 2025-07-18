# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

AutoC is a genetic programming system for autonomous aircraft control that uses the GPC++ kernel library. It evolves control programs for aircraft to follow flight paths using genetic programming techniques. The system includes simulation, evaluation, visualization, and blackbox data analysis capabilities.

## Build System

The autoc project uses CMake and requires several dependencies:

```bash
# Dependencies (Ubuntu 22.04)
sudo apt-get install -y libboost-all-dev libeigen3-dev libvtk9-dev xvfb g++ cmake gdb qtbase5-dev

# For headless operation (visualization)
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
```

**Key dependencies:**
- VTK 9.x for 3D visualization
- Eigen3 for linear algebra
- Boost (thread, system, serialization, log)
- AWS SDK for C++ (S3 storage)

## Architecture

**Core Components:**

1. **Genetic Programming Engine** (`autoc.cc/h`)
   - Implements GP operators: ADD, SUB, MUL, DIV, IF, EQ, GT, SIN, COS
   - Aircraft control functions: SETPITCH, SETROLL, SETTHROTTLE
   - Sensor functions: GETDPHI, GETDTHETA, GETDTARGET, GETDHOME, GETVEL
   - Fitness evaluation based on distance, alignment, and control smoothness

2. **Aircraft Simulation** (`minisim.cc/h`, `aircraft_state.h`)
   - `AircraftState`: Position, orientation (quaternion), velocity, control commands
   - Physics simulation with roll/pitch rate limits, throttle control
   - NED coordinate system (North-East-Down)

3. **Path Generation** (`pathgen.cc/h`)
   - Generates flight paths for GP evaluation
   - Methods: random, classic, line, computedPaths
   - Configurable via `autoc.ini`

4. **3D Visualization** (`renderer.cc/h`)
   - VTK-based rendering of flight paths, aircraft trajectories, and blackbox data
   - Multiple visualization layers: target paths (red), actual paths (blue/yellow ribbon), blackbox data (green/orange ribbon)
   - Interactive navigation with keyboard controls (N/P for next/previous generation)

5. **Blackbox Data Analysis**
   - Parses INAV blackbox CSV data (position in cm, attitude in centidegrees)
   - Converts to `AircraftState` objects for visualization
   - Handles coordinate transforms: NED to VTK display coordinates

**Data Flow:**
1. Path generation creates target waypoints
2. GP individuals control simulated aircraft to follow paths
3. Fitness evaluation scores performance (distance, alignment, control smoothness)
4. Results stored in S3 for visualization
5. Renderer displays evolution progress and flight trajectories

## Configuration

**Runtime configuration** (`autoc.ini`):
- `PopulationSize`: GP population size (default: 500)
- `NumberOfGenerations`: Evolution generations (default: 100)
- `SimNumPathsPerGeneration`: Number of paths per evaluation (default: 6)
- `EvalThreads`: Parallel evaluation threads (default: 12)
- `PathGeneratorMethod`: Path generation strategy
- `S3Bucket`/`S3Profile`: AWS S3 storage configuration

## Key Data Structures

**AircraftState**: Core aircraft representation
- Position: `Eigen::Vector3d` (NED coordinates)
- Orientation: `Eigen::Quaterniond` (aircraft attitude)
- Control commands: pitch, roll, throttle (-1 to 1)
- Time: simulation timestamp

**Path**: Target waypoint representation
- Start position and orientation
- Distance and angular measurements from start

**Coordinate Systems:**
- **NED**: North-East-Down (aircraft/simulation standard)
- **VTK**: Right-handed with +Z up (visualization)
- **Blackbox**: INAV format (position in cm, attitude in centidegrees)

## Visualization Controls

**Renderer keyboard controls:**
- `N`: Next generation
- `P`: Previous generation
- Mouse: Standard trackball camera controls

**Visualization layers:**
- **Red lines**: Target flight paths
- **Blue/Yellow ribbon**: Simulated aircraft trajectories (top/bottom faces)
- **Green/Orange ribbon**: Blackbox data trajectories (top/bottom faces)
- **Checkerboard plane**: Ground reference at Z=0

## Blackbox Data Integration

The renderer can visualize INAV blackbox data by:
1. Reading CSV with columns: `navPos[0]`, `navPos[1]`, `navPos[2]`, `attitude[0]`, `attitude[1]`, `attitude[2]`
2. Converting positions from centimeters to meters
3. Converting attitudes from centidegrees to radians
4. Creating `AircraftState` objects with proper quaternion orientation
5. Using same rendering pipeline as simulated data

Example usage:
```bash
cat blackbox.csv | ./renderer -k keyname
```

## AWS Integration

The system uses S3 for storing evolution results:
- Generation data stored as serialized objects
- Key format: `autoc-{timestamp}/gen{number}.dmp`
- Configurable bucket and profile in `autoc.ini`
- Supports local MinIO for development
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GPC++ is a C++ genetic programming kernel class library originally developed by Adam Fraser (1993-1994) and Thomas Weinbrenner (1996-1997). The library implements genetic programming techniques including automatically defined functions (ADFs), tournament selection, crossover, mutation, and population management.

This repository contains both the core GP library and various example applications, with `autoc/` being a modern aircraft control evolution system that includes a sophisticated bytecode interpreter system.

## Build System

The project uses a hybrid build system:

**Core Library (Make-based):**
```bash
# Build everything (library + traditional examples)
make

# Build only the core library
cd src && make

# Clean object files
make clean

# Clean everything including library and executables
make superclean

# Install library and headers (to /usr/local by default)
make install
```

**AutoC Project (CMake-based):**
```bash
# Dependencies (Ubuntu 22.04)
sudo apt-get install -y libboost-all-dev libeigen3-dev libvtk9-dev xvfb g++ cmake gdb qtbase5-dev

# Build autoc system (requires core library built first)
# IMPORTANT: Build from ~/GP/build directory, not ~/GP/autoc/build
cd ~/GP
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ../autoc
make

# For headless operation (visualization)
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
```

**Key Build Commands:**
```bash
# Full rebuild cycle
make clean && rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Debug ../autoc && cd .. && make

# Run evolution
./build/autoc

# Extract GP to bytecode
./build/gpextractor -b -o gp_program.dat

# Run evaluation mode
./build/autoc  # (with EvaluateMode=1 in autoc.ini)

# Visualize results
./build/renderer -k keyname
```

## Architecture

**Two-Tier Structure:**

**Core Library (`src/`):**
- Built as static library `lib/libgp.a`
- All modules compile to object files that are archived together
- Main headers: `include/gp.h`, `include/gpconfig.h`
- Pure GP functionality without domain-specific code

**Key Components:**
- Population management (`pop.cc`, `contain.cc`)
- Genetic operations (`cross.cc`, `mutate.cc`, `select.cc`)
- Tree nodes and genes (`node.cc`, `gene.cc`)
- Evaluation framework (`eval.cc`)
- Configuration system (`config.cc`)
- Load/save functionality (`loadsave.cc`)

**Traditional Example Applications:**
- `symbreg/`: Symbolic regression (evolves X⁴+X³+X²+X)
- `lawn/`: Lawn mower problem
- `ant/`: Santa Fe trail following
- `skeleton/`: Template for new problems

**Modern AutoC System (`autoc/`):**
- Aircraft control evolution using genetic programming
- Dual-mode architecture: GP tree evaluation + bytecode interpretation
- 3D visualization with VTK
- Real-time simulation and blackbox data analysis
- AWS S3 integration for result storage
- RPC-based distributed simulation

## AutoC Bytecode System

The AutoC system implements a sophisticated dual-mode evaluation architecture:

**Training Mode (GP Tree Evaluation):**
- Traditional recursive tree traversal during evolution
- Implemented in `autoc-eval.cc` with `MyGene::evaluate()`
- Direct function calls for each node type
- Serializes GP trees for transport to simulation nodes

**Evaluation Mode (Bytecode Interpretation):**
- Stack-based bytecode interpreter for deployment
- Implemented in `gp_bytecode.cc` with `GPBytecodeInterpreter::evaluate()`
- Compact binary representation for efficient transport
- Identical aircraft control semantics to GP mode

**Key Components:**
- `gpextractor`: Converts S3-stored GP trees to bytecode files
- `gp_bytecode.h/cc`: Stack-based interpreter with 256-element stack
- `minisim.cc`: Simulator that auto-detects GP tree vs bytecode data
- `autoc.cc`: Main evolution engine with configurable evaluation modes

**Critical PROGN Instruction:**
The PROGN bytecode instruction preserves GP tree semantics for operations with side effects:
```cpp
// Executes both children but returns only second child's value
// Essential for control sequences like PROGN(SETROLL(x), SETPITCH(y))
```

**RPC Communication:**
- Uses Boost.Asio TCP sockets between autoc and minisim
- Boost text serialization for GP trees
- Boost binary serialization for bytecode interpreters
- Automatic mode detection based on data format

## Key Configuration Files

**Core Library** (`Makefile.ini`):
- Compiler: g++ (configurable)
- Library: Static archive (`.a`) by default, can be changed to dynamic (`.so`)
- Flags: `-g -Wno-write-strings` by default

**AutoC System** (`autoc/autoc.ini`):
- `PopulationSize`: GP population size (default: 500)
- `NumberOfGenerations`: Evolution generations (default: 50)
- `SimNumPathsPerGeneration`: Number of paths per evaluation (default: 6)
- `EvalThreads`: Parallel evaluation threads (default: 1)
- `PathGeneratorMethod`: Path generation strategy (random, classic, line, computedPaths)
- `EvaluateMode`: 0=normal GP evolution, 1=bytecode verification
- `BytecodeFile`: Path to bytecode file for evaluation mode
- `S3Bucket`/`S3Profile`: AWS S3 storage configuration

## Development Dependencies

**Core Library:** Standard C++ compiler (no templates/exceptions required)

**AutoC System:** 
- VTK 9.x for 3D visualization
- Eigen3 for linear algebra
- Boost (thread, system, serialization, log)
- AWS SDK for C++ (S3 storage)
- Qt5 for GUI components

## AutoC Data Structures

**AircraftState**: Core aircraft representation
- Position: `Eigen::Vector3d` (NED coordinates)
- Orientation: `Eigen::Quaterniond` (aircraft attitude)
- Control commands: pitch, roll, throttle (-1 to 1)
- Time: simulation timestamp

**Path**: Target waypoint representation
- Start position and orientation
- Distance and angular measurements from start

**GPBytecode**: Bytecode instruction format
- 9-byte structure: opcode, argc, constant value
- Stack-based execution model
- Supports all GP operators and aircraft control functions

## GP Operators and Sensors

**Mathematical Operations:**
- `ADD, SUB, MUL, DIV` - Basic arithmetic
- `SIN, COS` - Trigonometric functions
- `IF, EQ, GT` - Conditional and comparison operations
- `CLAMP, ATAN2, ABS, SQRT, MIN, MAX` - Mathematical helper functions

**Aircraft Control:**
- `SETPITCH, SETROLL, SETTHROTTLE` - Set control commands (-1 to 1)
- `GETPITCH, GETROLL, GETTHROTTLE` - Get current control commands

**Navigation Sensors:**
- `GETDPHI(steps)` - Roll angle to target at path step offset
- `GETDTHETA(steps)` - Pitch angle to target at path step offset  
- `GETDTARGET(steps)` - Distance-based throttle estimate to target
- `GETDHOME` - Distance to home/origin point
- `GETVEL` - Current aircraft speed magnitude

**Attitude Sensors:**
- `GETROLL_RAD` - Current roll angle in radians
- `GETPITCH_RAD` - Current pitch angle in radians
- `GETALPHA` - Angle of attack (alpha)
- `GETBETA` - Sideslip angle (beta)

**Velocity Sensors:**
- `GETVELX` - North velocity component (NED)
- `GETVELY` - East velocity component (NED)  
- `GETVELZ` - Down velocity component (NED)

**Constants:**
- `PI, ZERO, ONE, TWO` - Mathematical constants
- `PROGN` - Sequential execution operator (returns second argument)

**Platform Compatibility:**
All mathematical functions use platform-specific macros (e.g., `CLAMP_DEF`) to support both full C++ builds and Arduino/embedded platforms.

## Coordinate Systems

- **NED**: North-East-Down (aircraft/simulation standard)
- **VTK**: Right-handed with +Z up (visualization)
- **Blackbox**: INAV format (position in cm, attitude in centidegrees)

## Common Workflows

**GP Evolution:**
1. Configure `autoc.ini` with EvaluateMode=0
2. Run `./build/autoc` to evolve controllers
3. Results stored in S3 with format: `autoc-{timestamp}/gen{number}.dmp`

**Bytecode Deployment:**
1. Extract best GP: `./build/gpextractor -b -o gp_program.dat`
2. Configure `autoc.ini` with EvaluateMode=1, BytecodeFile=gp_program.dat
3. Run `./build/autoc` to verify identical behavior

**Visualization:**
1. Start virtual display: `Xvfb :99 -screen 0 1024x768x24 &`
2. Export display: `export DISPLAY=:99`
3. Run renderer: `./build/renderer -k keyname`
4. Controls: N/P for next/previous generation, mouse for camera

**Debugging Mode Differences:**
- Both modes should produce identical fitness values
- Store pre-evaluation control values for comparison
- Use debug output to trace SETROLL/GETROLL operations if needed

## File Dependencies

All source files depend on the main headers. The dependency chain flows from the core library outward to the examples, which link against `libgp.a`. The autoc project has its own dependency tree and uses the core library as a foundation for aircraft control evolution.

The bytecode system creates an additional dependency path: GP trees (S3) → gpextractor → bytecode files → evaluation mode, enabling deployment independent of the full GP system.
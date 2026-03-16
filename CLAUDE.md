# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GPC++ is a C++ genetic programming kernel class library originally developed by Adam Fraser (1993-1994) and Thomas Weinbrenner (1996-1997). The library implements genetic programming techniques including automatically defined functions (ADFs), tournament selection, crossover, mutation, and population management.

This repository contains both the core GP library and various example applications, with `autoc/` being a modern aircraft control evolution system that includes a sophisticated bytecode interpreter system.

## Documentation Structure

| Location | Description |
|----------|-------------|
| `CLAUDE.md` | This file - main project guidance |
| `specs/BACKLOG.md` | Project backlog and TODO items |
| `specs/` | Feature specs (speckit workflow) |
| `docs/` | Reference documentation (loaded on-demand) |
| `autoc/specs/` | Design specifications (ZZZ- prefix = archived/done) |
| `.specify/memory/constitution.md` | Project constitution and principles |

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

# Build autoc system (standalone, no libgp.a dependency)
cd ~/GP/autoc
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

# For headless operation (visualization)
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
```

**Build Scripts (in `autoc/`):**
```bash
# Clean debug build - full rebuild with cmake reconfigure (for testing/debugging)
cd ~/GP/autoc && bash rebuild.sh

# Clean optimized build - full rebuild with -O3 (for performance runs)
cd ~/GP/autoc && bash rebuild-perf.sh

# Incremental build - just recompile changed files (fastest)
cd ~/GP/autoc/build && make

# Extract best NN weights from S3 archive
./build/nnextractor -k keyname -o nn_weights.dat -i autoc.ini

# Generate embedded C++ from NN weight file
./build/nn2cpp -i nn_weights.dat -o ~/xiao-gp/generated/nn_program_generated.cpp
```

**Key Build Commands:**
```bash
# Full rebuild cycle (equivalent to rebuild.sh)
cd ~/GP/autoc && rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Debug .. && make

# Run evolution (from ~/GP/autoc)
./build/autoc

# Run evaluation mode (with EvaluateMode=1 and NNWeightFile set in autoc.ini)
./build/autoc

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

**NN Evaluation (`autoc/`, spec 013-neuroevolution):**
- NNGenome: topology + weights + fitness + generation metadata
- Forward pass: feedforward with tanh activation, fast_tanh LUT (512 entries)
- Input: 22 normalized sensor values (temporal error history, quaternion attitude, velocity, alpha/beta, commands)
- Output: 3 control commands (pitch, roll, throttle) via tanh → [-1, 1]
- Serialization: custom binary "NN01" format (not Boost), detected by magic bytes
- S3 prefix: `nn-{timestamp}/` vs `autoc-{timestamp}/` for GP

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
- `ControllerType`: "GP" or "NN" (default: GP)
- `NNMutationSigma`: initial mutation sigma (default: 0.1)
- `NNCrossoverAlpha`: BLX-alpha blend factor (-1 = uniform random)
- `NNWeightFile`: weight file path for eval mode (default: nn_weights.dat)
- `NNInitMethod`: xavier or uniform (default: xavier)

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
- `GETDTARGET(steps)` - *(deprecated)* Composite throttle estimate: `CLAMP((distance-10)/speed, -1, 1)`. Replaced by GETDIST primitives. Opcode retained for backward compatibility.
- `GETDHOME` - Distance to home/origin point
- `GETVEL` - Current aircraft speed magnitude

**Distance Sensors** (see specs/012-distance-temporal-nodes):
- `GETDIST` - Raw Euclidean distance to rabbit (meters, nullary)
- `GETDIST_PREV(n)` - Buffered distance at history index n (meters, unary)
- `GETDIST_RATE` - Rate of distance change (m/s, nullary, clamped [-10,10])

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

**Temporal State** (derivative-style control):
- `GETDPHI_PREV(n)` - Previous dPhi error at history index n
- `GETDTHETA_PREV(n)` - Previous dTheta error at history index n
- `GETDPHI_RATE` - Rate of change of dPhi (rad/s)
- `GETDTHETA_RATE` - Rate of change of dTheta (rad/s)

**Platform Compatibility:**
All mathematical functions use platform-specific macros (e.g., `CLAMP_DEF`) to support both full C++ builds and Arduino/embedded platforms.

**Portable Evaluator** (`gp_evaluator_portable.cc/h`):
- Single codebase for desktop and embedded evaluation
- LUT-based trigonometry (512-element SIN table with interpolation)
- `GP_BUILD` define for full desktop builds
- `GP_TEST` define exposes internal functions for testing

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

**NN Evolution:**
1. Configure `autoc.ini` with ControllerType=NN and NN params
2. Run `./build/autoc` to evolve NN controllers
3. Results stored in S3 with format: `nn-{timestamp}/gen{number}.dmp`

**NN Deployment to Embedded:**
1. Extract weights: `./build/nnextractor -k keyname -o nn_weights.dat`
2. Generate code: `./build/nn2cpp -i nn_weights.dat -o ~/xiao-gp/generated/nn_program_generated.cpp`
3. Build xiao-gp with NN env: `pio run -e xiaoblesense_nn`

**Debugging Mode Differences:**
- Both modes should produce identical fitness values
- Store pre-evaluation control values for comparison
- Use debug output to trace SETROLL/GETROLL operations if needed

## File Dependencies

All source files depend on the main headers. The dependency chain flows from the core library outward to the examples, which link against `libgp.a`. The autoc project has its own dependency tree and uses the core library as a foundation for aircraft control evolution.

The bytecode system creates an additional dependency path: GP trees (S3) → gpextractor → bytecode files → evaluation mode, enabling deployment independent of the full GP system.

## Key Design Specs

| Spec | Description |
|------|-------------|
| `docs/COORDINATE_CONVENTIONS.md` | NED frame, quaternion conventions, Euler extraction |
| `autoc/specs/LAYERED_CONTROLLER.md` | Safety/strategy layer architecture |
| `autoc/specs/ZZZ-FASTMATH.md` | LUT-based trig, platform-portable math (completed) |
| `autoc/specs/ZZZ-TEMPORAL_STATE.md` | GETDPHI_PREV/RATE nodes (completed) |

## Active Technologies
- C++17 (CMake 3.10+, g++) + Eigen3 (quaternions, vectors), GoogleTest 1.14.0 (testing) (001-gp-eval-tests)
- N/A (in-memory test state only) (001-gp-eval-tests)
- C++17 (g++, CMake 3.10+) + Eigen3 (vectors, quaternions), Boost (serialization, logging), GoogleTest 1.14.0 (002-path-interpolation)
- N/A (in-memory state only) (002-path-interpolation)
- C++17 (g++, CMake 3.10+) + Eigen3 (vectors, quaternions), Boost (serialization, logging, threads) (003-variations-redux)
- N/A (in-memory state, S3 for evolution artifacts) (003-variations-redux)
- C++17 (g++, CMake 3.10+) + Eigen3 (vectors, quaternions), Boost (serialization, logging, threads) (005-entry-fitness-ramp)
- N/A (in-memory state, S3 for evolution artifacts) (005-entry-fitness-ramp)
- C++17 (g++, CMake 3.10+) + Eigen3 (vectors), Boost (serialization, logging), GoogleTest 1.14.0 (012-distance-temporal-nodes)
- N/A (in-memory ring buffers, S3 for evolution artifacts) (012-distance-temporal-nodes)
- C++17 (g++, CMake 3.10+) + Eigen3 (vectors, quaternions), Boost (serialization, logging, threads), GoogleTest 1.14.0 (013-neuroevolution)
- S3 for evolution archives, local files for diagnostics (data.dat, data.stc) (013-neuroevolution)
- C++17 (g++, CMake 3.10+) + Eigen3 (linear algebra), AWS SDK (S3), inih (config parser, vendored). Boost removed during refactoring. (014-nn-training-signal)

## Recent Changes
- 001-gp-eval-tests: Added C++17 (CMake 3.10+, g++) + Eigen3 (quaternions, vectors), GoogleTest 1.14.0 (testing)

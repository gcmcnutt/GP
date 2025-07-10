# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GPC++ is a C++ genetic programming kernel class library originally developed by Adam Fraser (1993-1994) and Thomas Weinbrenner (1996-1997). The library implements genetic programming techniques including automatically defined functions (ADFs), tournament selection, crossover, mutation, and population management.

This repository contains both the core GP library and various example applications, with `autoc/` being a modern aircraft control evolution system.

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
# Build autoc system (requires core library built first)
cd autoc
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

# For headless operation
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
```

The build process:
1. Creates `lib/` directory if it doesn't exist
2. Builds the core GP library (`libgp.a`) in `src/`
3. Traditional examples use Make and link against `libgp.a`
4. AutoC uses CMake and requires additional dependencies (VTK, Eigen3, Boost, AWS SDK)

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
- 3D visualization with VTK
- Real-time simulation and blackbox data analysis
- AWS S3 integration for result storage
- Uses CMake build system with modern C++ dependencies

Each example reads a `.ini` configuration file at runtime to set GP parameters without recompilation.

## Key Configuration

**Core Library** (`Makefile.ini`):
- Compiler: g++ (configurable)
- Library: Static archive (`.a`) by default, can be changed to dynamic (`.so`)
- Flags: `-g -Wno-write-strings` by default

**AutoC System** (`autoc/autoc.ini`):
- Population size, generations, crossover probability
- Simulation parameters (paths per generation, evaluation threads)
- AWS S3 storage configuration
- Path generation methods

## Development Dependencies

**Core Library:** Standard C++ compiler (no templates/exceptions required)

**AutoC System:** 
- VTK 9.x for 3D visualization
- Eigen3 for linear algebra
- Boost (thread, system, serialization, log)
- AWS SDK for C++ (S3 storage)
- Qt5 for GUI components

## File Dependencies

All source files depend on the main headers. The dependency chain flows from the core library outward to the examples, which link against `libgp.a`. The autoc project has its own dependency tree and uses the core library as a foundation for aircraft control evolution.
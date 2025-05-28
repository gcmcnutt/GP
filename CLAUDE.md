# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GPC++ is a C++ genetic programming kernel class library originally developed by Adam Fraser (1993-1994) and Thomas Weinbrenner (1996-1997). The library implements genetic programming techniques including automatically defined functions (ADFs), tournament selection, crossover, mutation, and population management.

## Build System

The project uses a traditional Make-based build system:

```bash
# Build everything (library + examples)
make

# Build only the core library
cd src && make

# Clean object files
make clean

# Clean everything including library and executables  
make superclean

# Strip debug information
make strip

# Install library and headers (to /usr/local by default)
make install
```

The build process:
1. Creates `lib/` directory if it doesn't exist
2. Builds the core GP library (`libgp.a`) in `src/`
3. Builds example applications in subdirectories

Configuration is in `Makefile.ini`:
- Compiler: g++ (configurable)
- Library: Static archive (`.a`) by default, can be changed to dynamic (`.so`)
- Flags: `-g -Wno-write-strings` by default

## Architecture

**Core Library (`src/`):**
- Built as static library `lib/libgp.a`
- All modules compile to object files that are archived together
- Main headers: `include/gp.h`, `include/gpconfig.h`

**Key Components:**
- Population management (`pop.cc`, `contain.cc`)
- Genetic operations (`cross.cc`, `mutate.cc`, `select.cc`)
- Tree nodes and genes (`node.cc`, `gene.cc`)
- Evaluation framework (`eval.cc`)
- Configuration system (`config.cc`)
- Load/save functionality (`loadsave.cc`)

**Example Applications:**
- `symbreg/`: Symbolic regression (evolves X⁴+X³+X²+X)
- `lawn/`: Lawn mower problem
- `ant/`: Santa Fe trail following
- `skeleton/`: Template for new problems
- `autoc/`: Additional example (uses CMake)

Each example reads a `.ini` configuration file at runtime to set GP parameters without recompilation.

## File Dependencies

All source files depend on the main headers. The dependency chain flows from the core library outward to the examples, which link against `libgp.a`.
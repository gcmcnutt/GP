# AutoC - Autonomous Aircraft Control via Genetic Programming

AutoC is a genetic programming system for autonomous aircraft control that uses the GPC++ kernel library. It evolves control programs for aircraft to follow flight paths using genetic programming techniques.

## Quick Start

See the main project [CLAUDE.md](../CLAUDE.md) for comprehensive build instructions and architecture overview.

### Build

```bash
# Dependencies (Ubuntu 22.04)
sudo apt-get install -y libboost-all-dev libeigen3-dev libvtk9-dev xvfb g++ cmake gdb qtbase5-dev

# Build from scratch (debug)
cd ~/GP/autoc && bash rebuild.sh

# Incremental build
cd ~/GP && make

# Run tests
cd ~/GP/build && ctest --output-on-failure
```

### Headless Operation

```bash
Xvfb :99 -screen 0 1024x768x24 &
export DISPLAY=:99
```

## Documentation

| Location | Description |
|----------|-------------|
| [../CLAUDE.md](../CLAUDE.md) | Main project guidance, build system, architecture |
| [../specs/BACKLOG.md](../specs/BACKLOG.md) | Project backlog and TODO items |
| [specs/](specs/) | Design specifications for autoc features |
| [../specs/](../specs/) | Feature specs (speckit workflow) |
| [../.specify/](../.specify/) | Project constitution and templates |

### Design Specs (autoc/specs/)

| Spec | Status | Description |
|------|--------|-------------|
| COORDINATE_CONVENTIONS.md | Active | NED/quaternion conventions |
| LAYERED_CONTROLLER.md | Active | Safety/strategy layer design |
| FASTMATH.md | Active | LUT-based trig functions |
| RAMP_LANDSCAPE.md | Active | Evaluation landscape analysis |
| PROFILING.md | Reference | Performance profiling notes |
| ZZZ-*.md | Archived | Completed or superseded specs |

## Key Components

- **autoc.cc** - Main GP evolution engine
- **minisim.cc** - Aircraft physics simulation
- **renderer.cc** - VTK-based 3D visualization
- **gpextractor** - GP tree to bytecode converter
- **gp_evaluator_portable.cc** - Portable evaluator (desktop + embedded)
- **gp_bytecode.cc** - Stack-based bytecode interpreter

## Configuration

Runtime configuration via `autoc.ini`:

```ini
PopulationSize=500
NumberOfGenerations=100
SimNumPathsPerGeneration=6
EvalThreads=12
PathGeneratorMethod=random
```

## AWS Integration

Results stored in S3:
- Key format: `autoc-{timestamp}/gen{number}.dmp`
- Configure bucket/profile in `autoc.ini`

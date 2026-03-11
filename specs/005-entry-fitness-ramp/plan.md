# Implementation Plan: Entry Position Variations with Intercept-Budget Fitness Scaling

**Branch**: `005-entry-fitness-ramp` | **Date**: 2026-03-10 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/005-entry-fitness-ramp/spec.md`

## Summary

Add entry position offsets to scenario variations and implement per-step intercept-budget fitness scaling. The scaling function estimates time-to-intercept from initial displacement/heading, then ramps the distance and attitude penalty from a floor (~0.1) to full (1.0) over that budget. This provides GP with gradient signal during the intercept phase rather than crushing all individuals with large distance errors.

## Technical Context

**Language/Version**: C++17 (g++, CMake 3.10+)
**Primary Dependencies**: Eigen3 (vectors, quaternions), Boost (serialization, logging, threads)
**Storage**: N/A (in-memory state, S3 for evolution artifacts)
**Testing**: GoogleTest 1.14.0 (unit tests in `autoc/tests/`)
**Target Platform**: Linux (DGX workstation), headless evolution
**Project Type**: Library + CLI application (autoc evolution engine)
**Performance Goals**: No measurable regression in per-generation evaluation time
**Constraints**: Must be backward-compatible (sigma=0 produces identical fitness to current system)
**Scale/Scope**: ~49 scenarios per generation, fitness computed per-scenario with per-step scaling

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Testing-First | PASS | Unit tests for intercept budget estimation and scaling function. Integration validation via evolution runs. |
| II. Build Stability | PASS | Changes in autoc only (plus crrcsim for position offset application). GP `make`, crrcsim `make`, xiao-gp `pio run` all verified. |
| III. Dual-Mode Parity | PASS | Fitness scaling is in the fitness computation loop (autoc.cc), not in GP evaluation or bytecode interpretation. Both modes use the same fitness function. No new GP operators added. |

## Project Structure

### Documentation (this feature)

```text
specs/005-entry-fitness-ramp/
├── plan.md              # This file
├── spec.md              # Feature specification (clarified)
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
autoc/
├── autoc.h              # MODIFY: Add INTERCEPT_SCALE_FLOOR/CEILING, ENTRY_SAFE_* constants
├── autoc.cc             # MODIFY: Add intercept budget estimation, apply scaling in fitness loop
├── variation_generator.h # MODIFY: Add entryNorthOffset, entryEastOffset, entryAltOffset; cylindrical generation
├── config_manager.cc    # MODIFY: Add EntryPositionSigma config parameter
├── minisim.h            # MODIFY: Add position offset fields to ScenarioMetadata (version bump to 6)
├── minisim.cc           # MODIFY: Apply position offset to initial aircraft position
└── tests/
    └── gp_evaluator_tests.cc  # MODIFY: Add intercept budget and scaling function tests

~/crsim/crrcsim-0.9.13/
├── src/global.h         # MODIFY: Add entryNorthOffset, entryEastOffset
├── src/global.cpp       # MODIFY: Initialize position offset globals
├── src/crrc_main.cpp    # MODIFY: Apply position offsets to posX, posY in initAirplaneState
└── src/mod_inputdev/inputdev_autoc/
    └── inputdev_autoc.cpp  # MODIFY: Read position offsets from ScenarioMetadata
```

**Structure Decision**: No new files needed. All changes modify existing files following established patterns (variation_generator → ScenarioMetadata → simulator application).

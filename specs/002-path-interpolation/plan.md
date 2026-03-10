# Implementation Plan: Path Interpolation & Evaluator Improvements

**Branch**: `002-path-interpolation` | **Date**: 2026-03-06 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/002-path-interpolation/spec.md`

## Summary

Fix the path indexing jitter bug by replacing discrete waypoint lookup (`getPathIndex()`) with time-based position interpolation. Verify that FITNESS_SIMPLIFY and RAMP_LANDSCAPE are working correctly.

**Primary change**: New `getInterpolatedTargetPosition()` function using binary search + linear lerp between waypoints, eliminating discrete jumps that cause erratic sensor values on real-time systems.

## Technical Context

**Language/Version**: C++17 (g++, CMake 3.10+)
**Primary Dependencies**: Eigen3 (vectors, quaternions), Boost (serialization, logging), GoogleTest 1.14.0
**Storage**: N/A (in-memory state only)
**Testing**: GoogleTest (`autoc/tests/gp_evaluator_tests.cc`)
**Target Platform**: Linux (desktop), embedded (xiao-gp via portable evaluator)
**Project Type**: Library (GP evaluator with portable core)
**Performance Goals**: <1% sensor value change from 10ms timing jitter
**Constraints**: Must work on both desktop (full Eigen) and embedded (portable math)
**Scale/Scope**: Single aircraft evaluation, paths with ~100-1000 waypoints

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Evidence |
|-----------|--------|----------|
| I. Testing-First | ✅ PASS | Spec defines 4 interpolation tests: InterpolationMidpoint, InterpolationContinuity, InterpolationJitterRobust, EarlyTimestampOvershoot fix |
| II. Build Stability | ✅ PASS | Success criteria #4: `cd ~/GP && make` succeeds; existing 77 tests must pass |
| III. Dual-Mode Parity | ✅ PASS | Changes in `gp_evaluator_portable.cc` are shared by both GP tree and bytecode modes |

**Gate result**: PASS - proceed to Phase 0

## Project Structure

### Documentation (this feature)

```text
specs/002-path-interpolation/
├── spec.md              # Feature specification (clarified)
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
autoc/
├── aircraft_state.h          # MODIFY: Remove getPathIndex(), add getMaxTimeMsec()
├── gp_evaluator_portable.h   # MODIFY: Add getInterpolatedTargetPosition() declaration
├── gp_evaluator_portable.cc  # MODIFY: Implement interpolation, update executeGetD*()
├── autoc-eval.cc             # VERIFY: Uses portable evaluator (no changes needed)
├── autoc.cc                  # VERIFY: FITNESS_SIMPLIFY, RAMP_LANDSCAPE implementation
├── config_manager.cc         # VERIFY: VariationRampStep config parsing
└── tests/
    └── gp_evaluator_tests.cc # MODIFY: Add interpolation tests, update existing tests
```

**Structure Decision**: Single project with portable evaluator core. All changes are in `autoc/` subdirectory. No new directories needed.

## Constitution Check (Post-Design)

*Re-evaluation after Phase 1 design completion*

| Principle | Status | Evidence |
|-----------|--------|----------|
| I. Testing-First | ✅ PASS | Tests defined in data-model.md validation rules, spec verification section |
| II. Build Stability | ✅ PASS | No new dependencies; existing build commands unchanged |
| III. Dual-Mode Parity | ✅ PASS | `getInterpolatedTargetPosition()` in portable evaluator shared by both modes |

**Gate result**: PASS - design is constitution-compliant

## Complexity Tracking

No constitution violations requiring justification.

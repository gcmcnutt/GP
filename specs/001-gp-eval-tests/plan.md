# Implementation Plan: GP Evaluator Regression Tests

**Branch**: `001-gp-eval-tests` | **Date**: 2026-03-06 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/001-gp-eval-tests/spec.md`

## Summary

Comprehensive regression test suite for the GP evaluator achieving 100% coverage of all 42 operator nodes in `allNodes[]`. Tests are clean slate validation using Eigen quaternion operations as ground truth, covering both tree evaluator and bytecode interpreter paths with bit-accurate comparisons. Emphasis on 3D quaternion transform testing for all orientation-dependent nodes.

## Technical Context

**Language/Version**: C++17 (CMake 3.10+, g++)
**Primary Dependencies**: Eigen3 (quaternions, vectors), GoogleTest 1.14.0 (testing)
**Storage**: N/A (in-memory test state only)
**Testing**: GoogleTest via CMake FetchContent, integrated with `ctest`
**Target Platform**: x86_64 (primary), ARM (cross-validation)
**Project Type**: Test suite extension (existing infrastructure)
**Performance Goals**: Test suite completes in <5 seconds
**Constraints**: FP32 (gp_scalar) precision throughout, bit-accurate comparisons
**Scale/Scope**: 42 nodes × 2 evaluation modes × multiple quaternion orientations = ~200+ test cases

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Evidence/Notes |
|-----------|--------|----------------|
| **I. Testing-First** | ✅ PASS | Feature IS testing - fully aligned. Extends existing test infrastructure. |
| **II. Build Stability** | ✅ PASS | Tests integrate with existing `make test` / `ctest`. No new build targets required. |
| **III. Dual-Mode Parity** | ✅ PASS | Tests explicitly cover BOTH tree evaluator AND bytecode interpreter per FR-001. |

**Build verification commands** (per constitution):
- `cd ~/GP && make` - includes test compilation
- `ctest` or `make test` - runs test suite

**Gate Result**: PASS - No violations. Proceed to Phase 0.

## Project Structure

### Documentation (this feature)

```text
specs/001-gp-eval-tests/
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output
├── quickstart.md        # Phase 1 output
└── tasks.md             # Phase 2 output (via /speckit.tasks)
```

### Source Code (repository root)

```text
autoc/
├── tests/
│   └── gp_evaluator_tests.cc    # Extended test file (primary target)
├── gp_evaluator_portable.h       # Evaluator interface under test
├── gp_evaluator_portable.cc      # Evaluator implementation under test
├── gp_bytecode.h                 # Bytecode interpreter interface
├── gp_bytecode.cc                # Bytecode interpreter implementation
├── autoc-eval.cc                 # allNodes[] definition (42 nodes)
├── aircraft_state.h              # AircraftState class
└── gp_types.h                    # gp_scalar, gp_vec3, gp_quat typedefs
```

**Structure Decision**: Single test file extension pattern. The existing `gp_evaluator_tests.cc` will be expanded with new test fixtures and helper functions. No new directories needed.

## Complexity Tracking

> No Constitution Check violations - this section is empty.

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| (none)    | -          | -                                   |

## Post-Design Constitution Re-Check

*Gate: Verified after Phase 1 design completion.*

| Principle | Status | Post-Design Evidence |
|-----------|--------|---------------------|
| **I. Testing-First** | ✅ PASS | Design includes test patterns, coverage tracking, and quickstart guide |
| **II. Build Stability** | ✅ PASS | Uses existing CMake/ctest integration, no new build targets |
| **III. Dual-Mode Parity** | ✅ PASS | BytecodeEval test suite mirrors TreeEval for all nodes |

**No new violations introduced.** Design maintains full constitution compliance.

## Phase Artifacts Summary

| Phase | Artifact | Status | Description |
|-------|----------|--------|-------------|
| 0 | [research.md](research.md) | ✅ Complete | 5 research questions resolved |
| 1 | [data-model.md](data-model.md) | ✅ Complete | 6 test framework entities |
| 1 | contracts/ | ⏭️ Skipped | N/A - test suite has no external interfaces |
| 1 | [quickstart.md](quickstart.md) | ✅ Complete | Build, run, debug instructions |
| 2 | [tasks.md](tasks.md) | ✅ Complete | 90 tasks across 9 phases |

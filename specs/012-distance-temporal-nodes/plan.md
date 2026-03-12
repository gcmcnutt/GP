# Implementation Plan: Distance Temporal Sensor Nodes

**Branch**: `012-distance-temporal-nodes` | **Date**: 2026-03-12 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/012-distance-temporal-nodes/spec.md`

## Summary

Add three new GP sensor nodes (GETDIST, GETDIST_PREV, GETDIST_RATE) providing raw distance-to-rabbit in meters with temporal history and derivative, following the proven GETDPHI_PREV/RATE pattern. Deprecate the composite GETDTARGET node from active training sets. Implementation touches 15 files across opcode enums, evaluators, bytecode tools, config, docs, and tests.

## Technical Context

**Language/Version**: C++17 (g++, CMake 3.10+)
**Primary Dependencies**: Eigen3 (vectors), Boost (serialization, logging), GoogleTest 1.14.0
**Storage**: N/A (in-memory ring buffers, S3 for evolution artifacts)
**Testing**: GoogleTest — existing suite at `autoc/tests/gp_evaluator_tests.cc` (101+ tests)
**Target Platform**: Linux (desktop training), Arduino/XIAO BLE (embedded deployment via portable evaluator)
**Project Type**: Library + simulation system
**Performance Goals**: No measurable overhead — distance already computed each tick for fitness; buffering adds one float write per tick
**Constraints**: Ring buffer must fit embedded memory (<100 bytes); opcodes must append to enum end for bytecode backward compatibility
**Scale/Scope**: 15 files modified, ~200 lines of new code, ~50 lines of test code

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| **I. Testing-First** | PASS | Unit tests for all 3 nodes across GP tree + bytecode evaluators (FR-013) |
| **II. Build Stability** | PASS | Incremental — new opcodes appended to enum end, no breaking changes. Build verified via `cd ~/GP && make` |
| **III. Dual-Mode Parity** | PASS | All nodes implemented in GP tree evaluator AND bytecode interpreter AND portable evaluator (FR-007). Bytecode parity verified via EvaluateMode=1 (SC-005) |

No violations. No complexity tracking needed.

## Project Structure

### Documentation (this feature)

```text
specs/012-distance-temporal-nodes/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: no unknowns, documents pattern analysis
├── data-model.md        # Phase 1: ring buffer and opcode data model
└── checklists/
    └── requirements.md  # Spec quality checklist
```

### Source Code (modified files)

```text
autoc/
├── autoc.h                          # Opcode enum: +GETDIST, GETDIST_PREV, GETDIST_RATE
├── autoc.ini                        # TrainingNodes: +GETDIST nodes, -GETDTARGET (deprecated)
├── autoc-eval.ini                   # Mirror autoc.ini changes
├── autoc-eval.cc                    # allNodes[] registration table
├── aircraft_state.h                 # distHistory_[] buffer, recordErrorHistory(), clearHistory(), getHistoricalDist()
├── minisim.cc                       # Record distance to history buffer before GP eval
├── gp_evaluator_portable.h          # Opcode enum mirror, function declarations
├── gp_evaluator_portable.cc         # evaluateGPOperator() + evaluateBytecodePortable() switch cases, execute functions
├── gp_bytecode.cc                   # (delegates to portable — verify no direct dispatch)
├── gpextractor.cc                   # GP-to-bytecode: child processing for GETDIST_PREV (unary)
├── bytecode2cpp.cc                  # Stack analysis, operator names, code generation
└── tests/
    └── gp_evaluator_tests.cc        # Unit tests for all 3 nodes

CLAUDE.md                            # GP Operators documentation update
```

**Structure Decision**: No new files created. All changes are additions to existing files, following the exact pattern established by GETDPHI_PREV/GETDPHI_RATE in the same files.

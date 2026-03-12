# Tasks: Distance Temporal Sensor Nodes

**Input**: Design documents from `/specs/012-distance-temporal-nodes/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md

**Tests**: Tests are REQUIRED per constitution principle I (Testing-First) and FR-013.

**Organization**: Tasks grouped by user story. US1 (raw distance nodes) is the MVP. US2 (GETDTARGET deprecation) and US3 (cross-platform consistency) build on US1.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story (US1, US2, US3)

---

## Phase 1: Setup

**Purpose**: No new project setup needed — all changes are to existing files.

(No tasks — existing project, incremental changes only.)

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Opcode enum and history buffer infrastructure that ALL user stories depend on.

- [ ] T001 Add GETDIST, GETDIST_PREV, GETDIST_RATE to Operators enum before _END marker in autoc/autoc.h
- [ ] T002 Mirror enum additions in autoc/gp_evaluator_portable.h — same order, same values as autoc.h
- [ ] T003 Add distHistory_[HISTORY_SIZE] array to AircraftState in autoc/aircraft_state.h. Extend recordErrorHistory() to accept and store distance parameter. Extend clearHistory() to zero distHistory_[]. Add getHistoricalDist(int n) accessor following getHistoricalDPhi() pattern.
- [ ] T004 Update minisim.cc to compute distance-to-rabbit (Euclidean norm of craftPos minus interpolated rabbit position from PathProvider, same distance used in fitness calculation) and pass it to recordErrorHistory() before GP evaluation, at the same call site as existing dPhi/dTheta recording (~line 240)

**Checkpoint**: Enum defined, history buffer records distance each tick. No evaluation yet.

---

## Phase 3: User Story 1 — Raw Distance Sensor with Derivative Control (Priority: P1)

**Goal**: GP trees can use GETDIST, GETDIST_PREV(n), and GETDIST_RATE to sense and react to distance.

**Independent Test**: Unit tests pass for all 3 nodes; short evolution run uses distance nodes in GP trees.

### Tests for User Story 1

- [ ] T005 [P] [US1] Write unit tests for GETDIST in autoc/tests/gp_evaluator_tests.cc: verify returns Euclidean distance in meters for known aircraft/rabbit positions
- [ ] T006 [P] [US1] Write unit tests for GETDIST_PREV in autoc/tests/gp_evaluator_tests.cc: verify history lookup for indices 0-9, buffer wraparound, first-tick (return 0.0), index clamping (n > historyCount, negative n, fractional n)
- [ ] T007 [P] [US1] Write unit tests for GETDIST_RATE in autoc/tests/gp_evaluator_tests.cc: verify rate = (dist[0]-dist[1])/dt for known sequences, closing (negative) vs opening (positive) sign convention, first-tick returns 0.0, clamping to [-10,10] m/s, dt default when timestamps equal

### Implementation for User Story 1

- [ ] T008 [P] [US1] Declare executeGetDist(), executeGetDistPrev(), executeGetDistRate() in autoc/gp_evaluator_portable.h
- [ ] T009 [US1] Implement executeGetDist() in autoc/gp_evaluator_portable.cc: compute Euclidean distance from aircraftState position to rabbit's interpolated position (same distance used in fitness calc)
- [ ] T010 [US1] Implement executeGetDistPrev() in autoc/gp_evaluator_portable.cc: look up distHistory_ via getHistoricalDist(int(arg)), following executeGetDPhiPrev() pattern
- [ ] T011 [US1] Implement executeGetDistRate() in autoc/gp_evaluator_portable.cc: compute (dist[0]-dist[1])/dt using timestamp-based dt (same formula as executeGetDPhiRate()), clamp to [-10,10] m/s, return 0.0 if historyCount < 2
- [ ] T012 [US1] Add switch cases for GETDIST, GETDIST_PREV, GETDIST_RATE in evaluateGPOperator() in autoc/gp_evaluator_portable.cc — GETDIST and GETDIST_RATE as zero-arg terminals, GETDIST_PREV as unary
- [ ] T013 [US1] Add switch cases for GETDIST, GETDIST_PREV, GETDIST_RATE in evaluateBytecodePortable() in autoc/gp_evaluator_portable.cc — GETDIST_PREV as unary (pop 1, push 1), GETDIST and GETDIST_RATE as nullary (push 1)
- [ ] T014 [US1] Add entries for GETDIST (0 args), GETDIST_PREV (1 arg), GETDIST_RATE (0 args) to allNodes[] table in autoc/autoc-eval.cc with correct name strings
- [ ] T015 [US1] Run tests: `cd ~/GP && make && ./build/autoc_tests` — verify all existing tests pass plus new GETDIST/GETDIST_PREV/GETDIST_RATE tests

**Checkpoint**: All 3 distance nodes functional in GP tree and bytecode evaluation. Tests pass.

---

## Phase 4: User Story 2 — GETDTARGET Deprecation (Priority: P1)

**Goal**: GETDTARGET removed from active training configuration, replaced by GETDIST primitives. Backward compatibility preserved.

**Independent Test**: autoc.ini parsed without GETDTARGET in node set; existing GETDTARGET bytecode still evaluates correctly.

### Implementation for User Story 2

- [ ] T016 [US2] Update TrainingNodes in autoc/autoc.ini: remove GETDTARGET, add GETDIST, GETDIST_PREV, GETDIST_RATE. Add comment: `# GETDTARGET deprecated — composite signal (distance-10)/speed replaced by cleaner GETDIST primitives`
- [ ] T017 [US2] Update TrainingNodes in autoc/autoc-eval.ini: same changes as autoc.ini
- [ ] T018 [US2] Update node set reference comment block in autoc/autoc.ini: add GETDIST/GETDIST_PREV/GETDIST_RATE to DISTANCE SENSORS category, move GETDTARGET to DEPRECATED section with explanation
- [ ] T019 [US2] Update node set reference comment block in autoc/autoc-eval.ini: same changes as autoc.ini
- [ ] T020 [US2] Write backward-compat test in autoc/tests/gp_evaluator_tests.cc: verify GETDTARGET opcode still evaluates correctly (returns CLAMP((distance-10)/relVel, -1, 1)) even though removed from TrainingNodes

**Checkpoint**: Config files updated. GETDTARGET evaluates correctly but is not in active training set.

---

## Phase 5: User Story 3 — Cross-Platform Consistency (Priority: P2)

**Goal**: All evaluation paths handle the new nodes identically. Bytecode tools generate correct code.

**Independent Test**: Bytecode round-trip produces identical outputs to GP tree evaluation for all 3 nodes.

### Tests for User Story 3

- [ ] T021 [P] [US3] Write bytecode round-trip tests in autoc/tests/gp_evaluator_tests.cc: for each of GETDIST, GETDIST_PREV, GETDIST_RATE, verify bytecode evaluation matches GP tree evaluation for known input sequences

### Implementation for User Story 3

- [ ] T022 [US3] Add GETDIST_PREV case to generateBytecode() in autoc/gpextractor.cc: handle as unary (1 child to process), following GETDPHI_PREV pattern
- [ ] T023 [US3] Add GETDIST, GETDIST_RATE cases to generateBytecode() in autoc/gpextractor.cc: handle as nullary (no children), following GETDPHI_RATE pattern
- [ ] T024 [US3] Update analyzeStackDepth() in autoc/bytecode2cpp.cc: add GETDIST_PREV as unary (net 0 stack), GETDIST and GETDIST_RATE as nullary (push 1)
- [ ] T025 [US3] Update getOperatorName() in autoc/bytecode2cpp.cc: add string mappings for all 3 new opcodes
- [ ] T026 [US3] Update generateInstruction() in autoc/bytecode2cpp.cc: add cases for GETDIST (nullary), GETDIST_PREV (unary), GETDIST_RATE (nullary) following existing patterns
- [ ] T027 [US3] Run full test suite and verify bytecode round-trip tests pass: `cd ~/GP && make && ./build/autoc_tests`

**Checkpoint**: All bytecode tools handle new nodes. Cross-platform parity verified.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Documentation and build verification

- [ ] T028 [P] Update CLAUDE.md GP Operators section: add GETDIST, GETDIST_PREV, GETDIST_RATE descriptions under new "Distance Sensors" category. Add deprecation note to GETDTARGET entry in "Navigation Sensors".
- [ ] T029 [P] Update CLAUDE.md Active Technologies section for 012-distance-temporal-nodes
- [x] T030 [P] Prefix eval-mode output files: auto-prefix data.dat/data.stc with "eval-" when EvaluateMode=1 in autoc/autoc.cc (backlog item: Configurable Output File Prefixes)
- [ ] T031 Verify build stability: `cd ~/GP && make` succeeds, `./build/autoc_tests` passes all tests (existing + new)
- [ ] T032 Run short evolution (5 generations) with updated TrainingNodes to verify GETDIST/GETDIST_PREV/GETDIST_RATE appear in evolved GP trees — MANUAL validation

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 2 (Foundational)**: No dependencies — start immediately
- **Phase 3 (US1)**: Depends on Phase 2 (enum + history buffer)
- **Phase 4 (US2)**: Depends on Phase 3 (nodes must be functional before config switch)
- **Phase 5 (US3)**: Depends on Phase 3 (nodes must exist for bytecode tools)
- **Phase 6 (Polish)**: Depends on Phases 3-5

### User Story Dependencies

- **US1 (P1)**: Depends on Foundational only — core implementation
- **US2 (P1)**: Depends on US1 — can't deprecate GETDTARGET until replacements work
- **US3 (P2)**: Depends on US1 — can't test bytecode parity without working nodes. Can run in parallel with US2.

### Within Each Phase

- T001, T002 can run in parallel (different files)
- T005, T006, T007 can run in parallel (same file but independent test functions)
- T008 can run in parallel with T005-T007 (header vs test file)
- T022-T026 can run in parallel (different files: gpextractor.cc, bytecode2cpp.cc)
- T028, T029 can run in parallel (same file but independent sections)

---

## Parallel Example: User Story 1

```bash
# After Phase 2 completes, launch tests and header in parallel:
Task T005: "Unit test GETDIST in autoc/tests/gp_evaluator_tests.cc"
Task T006: "Unit test GETDIST_PREV in autoc/tests/gp_evaluator_tests.cc"
Task T007: "Unit test GETDIST_RATE in autoc/tests/gp_evaluator_tests.cc"
Task T008: "Declare functions in autoc/gp_evaluator_portable.h"

# Then implement sequentially (same file):
Task T009-T013: "Implement + dispatch in autoc/gp_evaluator_portable.cc"
Task T014: "Register in autoc/autoc-eval.cc"
Task T015: "Build and run tests"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 2: Enum + history buffer (T001-T004)
2. Complete Phase 3: GETDIST/GETDIST_PREV/GETDIST_RATE nodes + tests (T005-T015)
3. **STOP and VALIDATE**: All tests pass, nodes work in GP tree and bytecode
4. Can run evolution immediately with manually edited TrainingNodes

### Incremental Delivery

1. Phase 2 → Foundational infrastructure ready
2. Phase 3 (US1) → Distance nodes work → Test independently
3. Phase 4 (US2) → Config updated, GETDTARGET deprecated → Ready for training runs
4. Phase 5 (US3) → Bytecode tools updated → Deployment-ready
5. Phase 6 → Docs and final verification

---

## Notes

- All tasks follow the existing GETDPHI_PREV/GETDPHI_RATE pattern — reference those implementations directly
- The portable evaluator (gp_evaluator_portable.cc) is the single source of truth for evaluation logic — both GP tree and bytecode paths delegate to it
- Opcode enum ordering is critical: append before _END, never reorder existing entries
- recordErrorHistory() signature change (adding distance parameter) affects both aircraft_state.h and all call sites in minisim.cc

# Tasks: Path Interpolation & Evaluator Improvements

**Input**: Design documents from `/specs/002-path-interpolation/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md

**Tests**: Tests are REQUIRED per Constitution Principle I (Testing-First).

**Organization**: Tasks organized by feature priority from spec.md scope table.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (US1, US2, US3)
- Include exact file paths in descriptions

## Path Conventions

All changes in `autoc/` subdirectory per plan.md structure decision.

---

## Phase 1: Setup

**Purpose**: No setup required - existing project structure, no new dependencies

- [ ] T001 Verify build passes before starting: `cd ~/GP && make`

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Header changes that other tasks depend on

**⚠️ CRITICAL**: US1 implementation depends on these completions

- [ ] T002 Add `getMaxTimeMsec()` method to PathProvider in `autoc/aircraft_state.h`
- [ ] T003 Add `MAX_OFFSET_STEPS` constant (10) in `autoc/aircraft_state.h`
- [ ] T004 Add `getInterpolatedTargetPosition()` declaration in `autoc/gp_evaluator_portable.h`

**Checkpoint**: Header declarations ready - implementation can now begin

---

## Phase 3: User Story 1 - Path Interpolation (Priority: P1) 🎯 MVP

**Goal**: Replace discrete waypoint lookup with time-based interpolation, eliminating sensor jumps from timing jitter

**Independent Test**: Run `autoc_tests --gtest_filter=NavigationOps.Interpolation*` - all pass with <1% sensor change from 10ms jitter

### Tests for User Story 1 ⚠️

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T005 [P] [US1] Add TEST(NavigationOps, InterpolationMidpoint) in `autoc/tests/gp_evaluator_tests.cc` - verify lerp at 50% between waypoints
- [ ] T006 [P] [US1] Add TEST(NavigationOps, InterpolationContinuity) in `autoc/tests/gp_evaluator_tests.cc` - verify no jumps at waypoint boundaries
- [ ] T007 [P] [US1] Add TEST(NavigationOps, InterpolationJitterRobust) in `autoc/tests/gp_evaluator_tests.cc` - verify 10ms jitter causes <1% change
- [ ] T008 [P] [US1] Add TEST(NavigationOps, InterpolationBoundaryClamp) in `autoc/tests/gp_evaluator_tests.cc` - verify ±10 step clamp (±1 second)
- [ ] T009 [P] [US1] Add TEST(NavigationOps, InterpolationNaNHandling) in `autoc/tests/gp_evaluator_tests.cc` - verify NaN input returns current rabbit position
- [ ] T010 [US1] Verify new tests fail (no implementation yet): `cd ~/GP && make && ./build/autoc_tests --gtest_filter=NavigationOps.Interpolation*`

### Implementation for User Story 1

- [ ] T011 [US1] Implement `getInterpolatedTargetPosition()` with binary search + linear lerp in `autoc/gp_evaluator_portable.cc`
- [ ] T012 [US1] Update `executeGetDPhi()` to use interpolated position in `autoc/gp_evaluator_portable.cc`
- [ ] T013 [US1] Update `executeGetDTheta()` to use interpolated position in `autoc/gp_evaluator_portable.cc`
- [ ] T014 [US1] Update `executeGetDTarget()` to use interpolated position in `autoc/gp_evaluator_portable.cc`
- [ ] T015 [US1] Update existing TEST(NavigationOps, EarlyTimestampOvershoot) to verify fix in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T016 [US1] Remove `getPathIndex()` from `autoc/aircraft_state.h`
- [ ] T017 [US1] Verify all interpolation tests pass: `cd ~/GP && make && ./build/autoc_tests --gtest_filter=NavigationOps.*`
- [ ] T018 [US1] Verify all 77+ existing tests still pass: `cd ~/GP && make && ./build/autoc_tests`

**Checkpoint**: Path interpolation complete - sensors now continuous, jitter-robust

---

## Phase 4: User Story 2 - Verify FITNESS_SIMPLIFY (Priority: P2)

**Goal**: Confirm 2-objective fitness (distance + attitude) with path-relative scaling is working correctly

**Independent Test**: Code inspection confirms DISTANCE_POWER=1.2 used, attitude_scale computed per-path; unit test for edge cases passes

### Tests for User Story 2 ⚠️

- [ ] T019 [P] [US2] Add TEST(Fitness, AttitudeScaleComputation) in `autoc/tests/gp_evaluator_tests.cc` - verify `attitude_scale = path_distance / max(path_turn_rad, 1.0f)` for various path geometries
- [ ] T020 [P] [US2] Add TEST(Fitness, AttitudeScaleStraightLine) in `autoc/tests/gp_evaluator_tests.cc` - verify straight-line path (turn_rad→0) uses clamp of 1.0f

### Implementation for User Story 2

- [ ] T021 [US2] Code inspection: Confirm `DISTANCE_POWER = 1.2` defined in `autoc/autoc.h:19`
- [ ] T022 [US2] Code inspection: Confirm `attitude_scale` computed per-path in `autoc/autoc.cc` fitness calculation
- [ ] T023 [US2] Code inspection: Verify removed metrics (direction alignment, energy deviation) are NOT computed in `autoc/autoc.cc`
- [ ] T024 [US2] Verify fitness tests pass: `cd ~/GP && make && ./build/autoc_tests --gtest_filter=Fitness.*`
- [ ] T025 [US2] Run short evolution (~10 gens, small pop) and verify fitness values reasonable (not NaN, not exploding)

**Checkpoint**: FITNESS_SIMPLIFY verified working

---

## Phase 5: User Story 3 - Verify RAMP_LANDSCAPE (Priority: P2)

**Goal**: Confirm gradual variation scaling ramps from 0→1 over training generations

**Independent Test**: Unit test for scale computation passes; evolution run shows "VariationScale" increasing in logs

### Tests for User Story 3 ⚠️

- [ ] T026 [P] [US3] Add TEST(Ramp, VariationScaleDisabled) in `autoc/tests/gp_evaluator_tests.cc` - verify rampStep=0 returns scale=1.0
- [ ] T027 [P] [US3] Add TEST(Ramp, VariationScaleGenZero) in `autoc/tests/gp_evaluator_tests.cc` - verify gen=0 returns scale=0.0
- [ ] T028 [P] [US3] Add TEST(Ramp, VariationScaleMidTraining) in `autoc/tests/gp_evaluator_tests.cc` - verify gen=50/100 returns scale≈0.5
- [ ] T029 [P] [US3] Add TEST(Ramp, VariationScaleEndTraining) in `autoc/tests/gp_evaluator_tests.cc` - verify gen=99/100 returns scale≈0.95

### Implementation for User Story 3

- [ ] T030 [US3] Code inspection: Confirm `computeVariationScale()` implementation matches spec in `autoc/autoc.cc:116-125`
- [ ] T031 [US3] Code inspection: Confirm scale applied to all variation sigmas in `autoc/autoc.cc:populateVariationOffsets()`
- [ ] T032 [US3] Code inspection: Confirm `VariationRampStep` config loaded in `autoc/config_manager.cc`
- [ ] T033 [US3] Verify ramp tests pass: `cd ~/GP && make && ./build/autoc_tests --gtest_filter=Ramp.*`
- [ ] T034 [US3] Run short evolution (~20 gens) with VariationRampStep=5 and verify log shows "VariationScale" increasing
- [ ] T035 [US3] Test with VariationRampStep=0 and verify full sigmas used from start (no ramp)

**Checkpoint**: RAMP_LANDSCAPE verified working

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Final verification and documentation

- [ ] T036 Verify full test suite passes: `cd ~/GP && make && ./build/autoc_tests`
- [ ] T037 Verify build on clean checkout: `cd ~/GP/autoc && bash rebuild.sh`
- [ ] T038 Update spec.md status from "Clarified" to "Implemented" in `specs/002-path-interpolation/spec.md`
- [ ] T039 [P] Mark feature complete in BACKLOG.md - move "Path Interpolation" item to completed section in `specs/BACKLOG.md`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - verify build works
- **Foundational (Phase 2)**: Depends on Setup - adds header declarations
- **US1 Path Interpolation (Phase 3)**: Depends on Foundational - BLOCKS US2/US3
- **US2 FITNESS_SIMPLIFY (Phase 4)**: Can start after US1 (clean sensor signals needed)
- **US3 RAMP_LANDSCAPE (Phase 5)**: Can start after US1 (or parallel with US2)
- **Polish (Phase 6)**: Depends on US1, US2, US3 completion

### User Story Dependencies

- **US1 (P1)**: Must complete first - provides jitter-robust foundation
- **US2 (P2)**: Depends on US1 - clean fitness signals needed for verification
- **US3 (P2)**: Depends on US1 - can run parallel with US2 if desired

### Within Each User Story

- Tests MUST be written and FAIL before implementation
- Implementation follows test order
- Verify tests pass after implementation
- Story complete before moving to next

### Parallel Opportunities

**Phase 2 (Foundational)**: T002, T003, T004 touch different areas of headers - could be sequential in single file edit session

**Phase 3 (US1 Tests)**: T005-T009 are all new test functions - can write in parallel batches

**Phase 4 (US2 Tests)**: T019, T020 are independent test functions - can write in parallel

**Phase 5 (US3 Tests)**: T026-T029 are independent test functions - can write in parallel

---

## Parallel Example: User Story 1 Tests

```bash
# Launch all US1 tests together (write to same file but different TEST blocks):
Task: "T005 [P] [US1] Add TEST(NavigationOps, InterpolationMidpoint)"
Task: "T006 [P] [US1] Add TEST(NavigationOps, InterpolationContinuity)"
Task: "T007 [P] [US1] Add TEST(NavigationOps, InterpolationJitterRobust)"
Task: "T008 [P] [US1] Add TEST(NavigationOps, InterpolationBoundaryClamp)"
Task: "T009 [P] [US1] Add TEST(NavigationOps, InterpolationNaNHandling)"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (verify build)
2. Complete Phase 2: Foundational (header declarations)
3. Complete Phase 3: User Story 1 (path interpolation)
4. **STOP and VALIDATE**: Run all tests, verify <1% jitter sensitivity
5. Can deploy/merge MVP at this point

### Incremental Delivery

1. Complete Setup + Foundational → Foundation ready
2. Add US1 → Test independently → **MVP complete!**
3. Add US2 → Verify FITNESS_SIMPLIFY → Document findings
4. Add US3 → Verify RAMP_LANDSCAPE → Document findings
5. Each story adds verification without breaking previous work

### Suggested Order

```
T001 → T002 → T003 → T004           # Foundation
T005-T009 (parallel) → T010         # US1 tests (fail first)
T011 → T012 → T013 → T014 → T015 → T016 → T017 → T018  # US1 implementation
T019-T020 (parallel) → T021 → T022 → T023 → T024 → T025  # US2
T026-T029 (parallel) → T030 → T031 → T032 → T033 → T034 → T035  # US3
T036 → T037 → T038 → T039           # Polish
```

---

## Notes

- [P] tasks = different TEST blocks, can batch write
- [Story] label maps task to specific user story for traceability
- US1 is MVP - can stop there if needed
- US2 and US3 are verification only (code already implemented)
- All test commands assume build from `~/GP` directory
- Constitution requires tests pass before merge

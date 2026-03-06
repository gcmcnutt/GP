# Tasks: GP Evaluator Regression Tests

**Input**: Design documents from `/specs/001-gp-eval-tests/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, quickstart.md

**Note**: This feature IS testing - all tasks implement test code. Tests are the deliverable.

**Organization**: Tasks grouped by user story. Each story covers a category of GP operator nodes.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different test suites, no dependencies)
- **[Story]**: Which user story this task belongs to (US1-US6)
- All paths relative to repository root

## Path Conventions

```text
autoc/
├── tests/
│   └── gp_evaluator_tests.cc    # Primary test file (extended)
├── gp_evaluator_portable.h       # May need GP_TEST exposure
└── gp_evaluator_portable.cc      # LUT functions to expose
```

---

## Phase 1: Setup (Test Infrastructure Foundation)

**Purpose**: Establish test helpers and infrastructure before implementing individual node tests

- [ ] T001 Add `extern` declarations for `allNodes[]` and `allNodesCount` in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T002 [P] Add `#include <set>` for coverage tracking in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T003 [P] Expose `fastSin`, `fastCos`, `fastAtan2` via `GP_TEST` ifdef in `autoc/gp_evaluator_portable.h`

---

## Phase 2: Foundational (Test Helpers - BLOCKING)

**Purpose**: Core test utilities that ALL user story tests depend on

**⚠️ CRITICAL**: No node tests can begin until these helpers are complete

- [ ] T004 Implement `QuatHelper` namespace with `fromEuler()` and `fromAxisAngle()` in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T005 Add `QuatHelper` flight attitude presets: `level()`, `pitchedUp()`, `pitchedDown()`, `bankedRight()`, `bankedLeft()`, `yawed()` in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T006 Add `QuatHelper::climbingTurn()` for combined rotations and `deg()` helper in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T007 Implement `NodeCoverageTracker` class with `markTested()`, `isTested()`, `getUntestedOps()` in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T008 Add `TEST_NODE(op)` macro for coverage registration in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T009 [P] Implement `ExpectedValue::sin()`, `cos()`, `atan2()` using exposed LUT functions in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T010 [P] Implement `ExpectedValue::dPhi()` and `dTheta()` using Eigen operations in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T011 [P] Implement `ExpectedValue::alpha()`, `beta()`, `rollRad()`, `pitchRad()` using Eigen in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T012 Extend `makeState()` overloads to accept position, orientation, velocity in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T013 Extend `makePath()` to accept `simTimeMsec` parameter in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T014 Add `TestPathProvider::addPath()`, `setCurrentIndex()`, `clear()` helper methods in `autoc/tests/gp_evaluator_tests.cc`

**Checkpoint**: Foundation ready - node tests can now begin

---

## Phase 3: User Story 1 - Test Mathematical Operations (Priority: P1) 🎯 MVP

**Goal**: Verify all 12 math operator nodes with bit-accurate comparisons

**Nodes**: ADD, SUB, MUL, DIV, SIN, COS, ATAN2, CLAMP, ABS, SQRT, MIN, MAX

**Independent Test**: `./autoc_tests --gtest_filter="MathOps.*"` - requires no path/aircraft setup

### Tree Evaluator Tests

- [ ] T015 [P] [US1] Implement `TEST(MathOps, AddSubMulDiv)` for basic arithmetic in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T016 [P] [US1] Implement `TEST(MathOps, DivisionByZero)` verifying DIV(x,0)=0 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T017 [P] [US1] Implement `TEST(MathOps, SinCosBitAccurate)` using LUT expected values in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T018 [P] [US1] Implement `TEST(MathOps, SinCosAllQuadrants)` covering 0, π/2, π, 3π/2 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T019 [P] [US1] Implement `TEST(MathOps, Atan2AllQuadrants)` for (+,+), (+,-), (-,+), (-,-) in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T020 [P] [US1] Implement `TEST(MathOps, ClampMinMax)` for CLAMP, MIN, MAX in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T021 [P] [US1] Implement `TEST(MathOps, AbsSqrt)` including SQRT(-1)=0 edge case in `autoc/tests/gp_evaluator_tests.cc`

### Bytecode Interpreter Tests

- [ ] T022 [US1] Implement `TEST(BytecodeMath, AddSubMulDiv)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T023 [US1] Implement `TEST(BytecodeMath, SinCosBitAccurate)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T024 [US1] Implement `TEST(BytecodeMath, Atan2AllQuadrants)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T025 [US1] Implement `TEST(BytecodeMath, ClampMinMaxAbsSqrt)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`

**Checkpoint**: Math operations fully tested - 12 nodes covered

---

## Phase 4: User Story 2 - Test Control Side Effects (Priority: P1)

**Goal**: Verify all 6 control nodes with SET/GET cycle and clamping

**Nodes**: SETPITCH, SETROLL, SETTHROTTLE, GETPITCH, GETROLL, GETTHROTTLE

**Independent Test**: `./autoc_tests --gtest_filter="ControlOps.*"` - requires only AircraftState

### Tree Evaluator Tests

- [ ] T026 [P] [US2] Implement `TEST(ControlOps, SetPitchClampHigh)` verifying >1.0 clamps to 1.0 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T027 [P] [US2] Implement `TEST(ControlOps, SetRollClampLow)` verifying <-1.0 clamps to -1.0 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T028 [P] [US2] Implement `TEST(ControlOps, SetThrottleClamp)` verifying both bounds in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T029 [US2] Implement `TEST(ControlOps, SetGetCycle)` verifying SET followed by GET returns correct value in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T030 [US2] Implement `TEST(ControlOps, GetBeforeSet)` verifying initial GET values in `autoc/tests/gp_evaluator_tests.cc`

### Bytecode Interpreter Tests

- [ ] T031 [US2] Implement `TEST(BytecodeControl, SetPitchRollThrottle)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T032 [US2] Implement `TEST(BytecodeControl, GetPitchRollThrottle)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`

**Checkpoint**: Control operations fully tested - 6 nodes covered

---

## Phase 5: User Story 3 - Test Path Following Nodes (Priority: P1)

**Goal**: Verify path-following nodes with time offsets and 3D quaternion orientations

**Nodes**: GETDPHI, GETDTHETA, GETDTARGET

**Independent Test**: `./autoc_tests --gtest_filter="NavigationOps.*"` - requires TestPathProvider setup

### Tree Evaluator Tests - Identity Orientation Baseline

- [ ] T033 [P] [US3] Implement `TEST(NavigationOps, GetDPhiDThetaZeroOffset)` with target directly ahead in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T034 [P] [US3] Implement `TEST(NavigationOps, GetDThetaTargetAbove)` verifying positive result in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T035 [P] [US3] Implement `TEST(NavigationOps, GetDThetaTargetBelow)` verifying negative result in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T036 [P] [US3] Implement `TEST(NavigationOps, GetDPhiTargetRight)` verifying positive result in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T037 [P] [US3] Implement `TEST(NavigationOps, GetDPhiTargetLeft)` verifying negative result in `autoc/tests/gp_evaluator_tests.cc`

### Tree Evaluator Tests - Time Offset Indexing

- [ ] T038 [US3] Implement `TEST(NavigationOps, TimeOffsetPositive)` testing GETDPHI(2) targets waypoint+2 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T039 [US3] Implement `TEST(NavigationOps, TimeOffsetClampEnd)` testing offset beyond path length in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T040 [US3] Implement `TEST(NavigationOps, TimeOffsetClampStart)` testing negative offset clamps to 0 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T041 [US3] Implement `TEST(NavigationOps, GetDTarget)` verifying throttle estimate calculation in `autoc/tests/gp_evaluator_tests.cc`

### Tree Evaluator Tests - 3D Quaternion Orientations

- [ ] T042 [US3] Implement `TEST(NavigationOps, GetDPhiBankedRight)` using QuatHelper::bankedRight(deg(45)) in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T043 [US3] Implement `TEST(NavigationOps, GetDThetaPitchedUp)` using QuatHelper::pitchedUp(deg(30)) in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T044 [US3] Implement `TEST(NavigationOps, GetDPhiClimbingTurn)` using combined pitch+roll quaternion in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T045 [US3] Implement `TEST(NavigationOps, TargetAllOctants)` testing ±X, ±Y, ±Z target combinations in `autoc/tests/gp_evaluator_tests.cc`

### Bytecode Interpreter Tests

- [ ] T046 [US3] Implement `TEST(BytecodeNavigation, GetDPhiDTheta)` mirroring identity baseline tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T047 [US3] Implement `TEST(BytecodeNavigation, GetDTarget)` mirroring tree test in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T048 [US3] Implement `TEST(BytecodeNavigation, QuaternionOrientations)` mirroring 3D quaternion tests in `autoc/tests/gp_evaluator_tests.cc`

**Checkpoint**: Path following fully tested - 3 nodes covered with quaternion variations

---

## Phase 6: User Story 4 - Test Sensor Nodes (Priority: P1)

**Goal**: Verify sensor nodes with comprehensive 3D quaternion orientations

**Nodes**: GETVEL, GETVELX, GETVELY, GETVELZ, GETALPHA, GETBETA, GETROLL_RAD, GETPITCH_RAD, GETDHOME

**Independent Test**: `./autoc_tests --gtest_filter="SensorOps.*"` - requires AircraftState with velocity/orientation

### Tree Evaluator Tests - Velocity Sensors

- [ ] T049 [P] [US4] Implement `TEST(SensorOps, GetVelMagnitude)` verifying speed magnitude in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T050 [P] [US4] Implement `TEST(SensorOps, GetVelXYZ)` verifying NED velocity components in `autoc/tests/gp_evaluator_tests.cc`

### Tree Evaluator Tests - Alpha/Beta with Quaternions

- [ ] T051 [P] [US4] Implement `TEST(SensorOps, AlphaBetaLevel)` with identity orientation returning 0 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T052 [US4] Implement `TEST(SensorOps, AlphaPitchedUp)` using QuatHelper::pitchedUp(deg(30)) in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T053 [US4] Implement `TEST(SensorOps, BetaYawed)` using QuatHelper::yawed(deg(45)) in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T054 [US4] Implement `TEST(SensorOps, AlphaBetaRolled90)` using QuatHelper::bankedRight(deg(90)) in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T055 [US4] Implement `TEST(SensorOps, AlphaBetaCombined)` using combined pitch+roll+yaw in `autoc/tests/gp_evaluator_tests.cc`

### Tree Evaluator Tests - Roll/Pitch from Quaternion

- [ ] T056 [P] [US4] Implement `TEST(SensorOps, GetRollRadPurePitch)` verifying pure pitch orientations in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T057 [P] [US4] Implement `TEST(SensorOps, GetPitchRadPureRoll)` verifying pure roll orientations in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T058 [US4] Implement `TEST(SensorOps, GetRollPitchGimbalLock)` testing pitch=±90° edge case in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T059 [US4] Implement `TEST(SensorOps, GetRollPitchAllQuadrants)` covering 0°, 30°, 60°, 90°, 180° in `autoc/tests/gp_evaluator_tests.cc`

### Tree Evaluator Tests - Distance to Home

- [ ] T060 [US4] Implement `TEST(SensorOps, GetDHome)` verifying distance calculation in `autoc/tests/gp_evaluator_tests.cc`

### Bytecode Interpreter Tests

- [ ] T061 [US4] Implement `TEST(BytecodeSensors, GetVelXYZ)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T062 [US4] Implement `TEST(BytecodeSensors, AlphaBetaQuaternions)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T063 [US4] Implement `TEST(BytecodeSensors, GetRollPitchRad)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T064 [US4] Implement `TEST(BytecodeSensors, GetDHome)` mirroring tree test in `autoc/tests/gp_evaluator_tests.cc`

**Checkpoint**: Sensor operations fully tested - 9 nodes covered with quaternion variations

---

## Phase 7: User Story 5 - Test Temporal History Nodes (Priority: P2)

**Goal**: Verify temporal nodes with history buffer management

**Nodes**: GETDPHI_PREV, GETDTHETA_PREV, GETDPHI_RATE, GETDTHETA_RATE

**Independent Test**: `./autoc_tests --gtest_filter="TemporalOps.*"` - requires AircraftState with history

### Tree Evaluator Tests

- [ ] T065 [P] [US5] Implement `TEST(TemporalOps, GetDPhiPrevIndex0)` verifying most recent history in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T066 [P] [US5] Implement `TEST(TemporalOps, GetDPhiPrevIndex2)` verifying older history access in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T067 [P] [US5] Implement `TEST(TemporalOps, GetDThetaPrev)` mirroring dPhi tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T068 [US5] Implement `TEST(TemporalOps, EmptyHistoryReturnsZero)` verifying empty buffer behavior in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T069 [US5] Implement `TEST(TemporalOps, HistoryIndexOutOfBounds)` verifying clamping in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T070 [US5] Implement `TEST(TemporalOps, GetDPhiRate)` verifying rate calculation in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T071 [US5] Implement `TEST(TemporalOps, GetDThetaRate)` mirroring dPhi rate test in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T072 [US5] Implement `TEST(TemporalOps, RateWithZeroDt)` verifying default dt behavior in `autoc/tests/gp_evaluator_tests.cc`

### Bytecode Interpreter Tests

- [ ] T073 [US5] Implement `TEST(BytecodeTemporal, GetPrevNodes)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T074 [US5] Implement `TEST(BytecodeTemporal, GetRateNodes)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`

**Checkpoint**: Temporal operations fully tested - 4 nodes covered

---

## Phase 8: User Story 6 - Test Logic and Constant Nodes (Priority: P3)

**Goal**: Verify logic and constant nodes for completeness

**Nodes**: IF, EQ, GT, PROGN, ZERO, ONE, TWO, OP_PI

**Independent Test**: `./autoc_tests --gtest_filter="LogicOps.*"` - requires minimal setup

### Tree Evaluator Tests

- [ ] T075 [P] [US6] Implement `TEST(LogicOps, IfTruthyBranch)` verifying IF(1,a,b)=a in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T076 [P] [US6] Implement `TEST(LogicOps, IfFalsyBranch)` verifying IF(0,a,b)=b in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T077 [P] [US6] Implement `TEST(LogicOps, EqEqual)` verifying EQ(5,5)=1 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T078 [P] [US6] Implement `TEST(LogicOps, EqNotEqual)` verifying EQ(5,6)=0 in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T079 [P] [US6] Implement `TEST(LogicOps, GtComparisons)` verifying GT behavior in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T080 [US6] Implement `TEST(LogicOps, PrognSideEffect)` verifying PROGN(SETROLL,GETROLL) in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T081 [P] [US6] Implement `TEST(LogicOps, Constants)` verifying ZERO, ONE, TWO, OP_PI values in `autoc/tests/gp_evaluator_tests.cc`

### Bytecode Interpreter Tests

- [ ] T082 [US6] Implement `TEST(BytecodeLogic, IfEqGt)` mirroring tree tests in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T083 [US6] Implement `TEST(BytecodeLogic, PrognSideEffect)` mirroring tree test in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T084 [US6] Implement `TEST(BytecodeLogic, Constants)` mirroring tree test in `autoc/tests/gp_evaluator_tests.cc`

**Checkpoint**: Logic/constant operations fully tested - 8 nodes covered

---

## Phase 9: Coverage Validation & Polish

**Purpose**: Verify 100% node coverage and clean up

- [ ] T085 Implement `TEST(Coverage, AllNodesTested)` iterating allNodes[] and verifying each is tested in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T086 Implement `TEST(Coverage, NodeCount)` verifying 42 nodes total in `autoc/tests/gp_evaluator_tests.cc`
- [ ] T087 Run `cd ~/GP && make` and verify all tests compile without errors
- [ ] T088 Run `cd ~/GP/build && ctest --output-on-failure` and verify all tests pass
- [ ] T089 Verify test suite completes in <5 seconds per SC-006
- [ ] T090 Remove or refactor any existing tests that use EXPECT_NEAR for LUT functions (should be EXPECT_FLOAT_EQ)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies - can start immediately
- **Foundational (Phase 2)**: Depends on Setup - BLOCKS all user story tests
- **User Stories (Phase 3-8)**: All depend on Foundational phase completion
  - US1-US4 are all P1 priority - can proceed in parallel
  - US5 (P2) can proceed after foundational, independent of US1-US4
  - US6 (P3) can proceed after foundational, independent of others
- **Coverage (Phase 9)**: Depends on ALL user stories complete

### User Story Dependencies

| Story | Priority | Depends On | Independent |
|-------|----------|------------|-------------|
| US1 (Math) | P1 | Foundational | ✅ Yes |
| US2 (Control) | P1 | Foundational | ✅ Yes |
| US3 (Navigation) | P1 | Foundational | ✅ Yes |
| US4 (Sensors) | P1 | Foundational | ✅ Yes |
| US5 (Temporal) | P2 | Foundational | ✅ Yes |
| US6 (Logic) | P3 | Foundational | ✅ Yes |

All user stories are **independently testable** once Foundational phase completes.

### Within Each User Story

1. Tree evaluator tests first (primary)
2. Bytecode tests mirror tree tests (dual-mode parity)
3. Each test calls `TEST_NODE(op)` for coverage tracking

### Parallel Opportunities

**Phase 2 (Foundational)**:
- T009, T010, T011 can run in parallel (different ExpectedValue functions)

**Phase 3 (US1 - Math)**:
- T015-T021 can all run in parallel (different test cases)

**Phase 4 (US2 - Control)**:
- T026-T028 can run in parallel (different control axes)

**Phase 5 (US3 - Navigation)**:
- T033-T037 can run in parallel (different target positions)

**Phase 6 (US4 - Sensors)**:
- T049-T050, T051, T056-T057 can run in parallel (different sensor types)

**Phase 7 (US5 - Temporal)**:
- T065-T067 can run in parallel (different PREV indices)

**Phase 8 (US6 - Logic)**:
- T075-T079, T081 can run in parallel (different logic operations)

---

## Parallel Example: User Story 1 (Math)

```bash
# Launch tree evaluator tests in parallel:
Task: "T015 [P] [US1] Implement TEST(MathOps, AddSubMulDiv)"
Task: "T016 [P] [US1] Implement TEST(MathOps, DivisionByZero)"
Task: "T017 [P] [US1] Implement TEST(MathOps, SinCosBitAccurate)"
Task: "T018 [P] [US1] Implement TEST(MathOps, SinCosAllQuadrants)"
Task: "T019 [P] [US1] Implement TEST(MathOps, Atan2AllQuadrants)"
Task: "T020 [P] [US1] Implement TEST(MathOps, ClampMinMax)"
Task: "T021 [P] [US1] Implement TEST(MathOps, AbsSqrt)"

# Then bytecode tests (depend on tree tests for pattern):
Task: "T022 [US1] Implement TEST(BytecodeMath, AddSubMulDiv)"
Task: "T023 [US1] Implement TEST(BytecodeMath, SinCosBitAccurate)"
...
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T003)
2. Complete Phase 2: Foundational (T004-T014)
3. Complete Phase 3: User Story 1 - Math (T015-T025)
4. **STOP and VALIDATE**: `./autoc_tests --gtest_filter="MathOps.*:BytecodeMath.*"`
5. Verify 12 math nodes covered

### Incremental Delivery

1. Setup + Foundational → Test infrastructure ready
2. Add US1 (Math) → 12 nodes covered → Validate
3. Add US2 (Control) → 18 nodes covered → Validate
4. Add US3 (Navigation) → 21 nodes covered → Validate
5. Add US4 (Sensors) → 30 nodes covered → Validate
6. Add US5 (Temporal) → 34 nodes covered → Validate
7. Add US6 (Logic) → **42 nodes covered** → Full validation

### Parallel Team Strategy

With multiple developers:

1. Complete Setup + Foundational together
2. Once Foundational done:
   - Dev A: US1 (Math) + US6 (Logic) = 20 nodes
   - Dev B: US2 (Control) + US5 (Temporal) = 10 nodes
   - Dev C: US3 (Navigation) + US4 (Sensors) = 12 nodes
3. All stories complete → Coverage validation

---

## Summary

| Metric | Count |
|--------|-------|
| Total Tasks | 90 |
| Setup Tasks | 3 |
| Foundational Tasks | 11 |
| US1 (Math) Tasks | 11 |
| US2 (Control) Tasks | 7 |
| US3 (Navigation) Tasks | 16 |
| US4 (Sensors) Tasks | 16 |
| US5 (Temporal) Tasks | 10 |
| US6 (Logic) Tasks | 10 |
| Coverage/Polish Tasks | 6 |
| Parallel Opportunities | 48 tasks marked [P] |
| Nodes Covered | 42 (100%) |

---

## Notes

- All tests in single file `autoc/tests/gp_evaluator_tests.cc` per plan.md structure decision
- Each test calls `TEST_NODE(op)` to register coverage
- Tree evaluator tests use `evaluateGPOperator()`
- Bytecode tests use `evaluateBytecodePortable()`
- Use `EXPECT_FLOAT_EQ` for bit-accurate (LUT), `EXPECT_NEAR` for quaternion-based
- Commit after each phase or logical group
- Stop at any checkpoint to validate independently

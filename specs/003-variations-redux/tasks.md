# Tasks: Variations Redux

**Input**: Design documents from `/specs/003-variations-redux/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, quickstart.md

**Tests**: No new unit tests required. Validation is via training runs and log analysis. Existing 82 evaluator tests must continue passing after fitness constant changes.

**Organization**: Tasks grouped by user story (P1-P5), each independently testable via training runs.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Verify both repos build on feature branch, confirm existing test suite passes

- [x] T001 Verify GP repo builds on 003-variations-redux branch: `cd ~/GP/autoc && bash rebuild-perf.sh`
- [x] T002 Verify crrcsim builds on 003-variations-redux branch: `cd ~/crsim/crrcsim-0.9.13/build && cmake -DPERFORMANCE_BUILD=ON .. && make`
- [x] T003 Run existing evaluator tests pass: `cd ~/GP/build && ctest`

**Checkpoint**: Both repos build, all existing tests pass

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Update fitness constants that all subsequent training phases depend on

- [x] T004 Change DISTANCE_NORM from 5.0 to 2.0 in `autoc/autoc.h`
- [x] T005 Change DISTANCE_POWER from 1.5 to 2.0 in `autoc/autoc.h`
- [x] T006 Rebuild GP repo with updated constants: `cd ~/GP && make`
- [x] T007 Re-run evaluator tests to confirm no regressions: `cd ~/GP/build && ctest`

**Checkpoint**: Fitness constants updated, tests pass, ready for training validation

---

## Phase 3: User Story 1 - Tighter Distance Tracking (Priority: P1) MVP

**Goal**: Validate that DISTANCE_NORM=2.0, DISTANCE_POWER=2.0 produces tighter tracking than the 002 baseline (~15-20m)

**Independent Test**: Run 200-gen training with no variations. Compare tracking distance against 002 baseline. Target: <10m on straight segments.

### Implementation for User Story 1

- [x] T008 [US1] Ensure autoc.ini has variations disabled: EnableEntryVariations=0, EnableWindVariations=0, WindScenarios=1, RabbitSpeedSigma=0.0 in `autoc/autoc.ini`
- [x] T009 [US1] Run 200-gen training: `./build/autoc 2>&1 | tee log/autoc-003-p1.log`
- [x] T010 [US1] Analyze fitness convergence from `log/autoc-003-p1.log` — verify downward trend
- [x] T011 [US1] Check `data.stc` for evolved GP tree structure and `data.dat` for final fitness
- [x] T012 [US1] Render best-of-generation to visually verify tracking distance: `./build/renderer -k <keyname>`
- [x] T013 [US1] Document results: fitness value, approximate tracking distance (~10.5m median), comparison to 002 baseline

**Checkpoint**: Tighter tracking confirmed. If tracking is still >10m or crash rate increased, adjust DISTANCE_NORM/POWER before proceeding.

---

## Phase 4: User Story 2 - Wind Direction Variations (Priority: P2)

**Goal**: Enable wind direction variations in crrcsim and validate GP trains against varied wind conditions

**Independent Test**: Enable EnableWindVariations=1 with 9 scenarios. Run training, verify crrcsim applies different wind directions per scenario.

### Implementation for User Story 2

- [x] T014 [US2] Update autoc.ini: set EnableWindVariations=1, WindScenarios=36 in `autoc/autoc.ini`
- [x] T015 [US2] Run validation to verify crrcsim receives and logs wind offsets
- [x] T016 [US2] Verify autoc startup log shows variation table with non-zero wind direction offsets
- [x] T017 [US2] Run 200-gen training (gen 9800, fitness 1.81M): `log/autoc-002-wind1.log`
- [x] T018 [US2] Analyze fitness convergence — 0 divergences, gen23 breakthrough, smooth across all 36 paths
- [x] T019 [US2] Bytecode eval pipeline validated: extract→eval produces identical fitness (1812229.202944)
- [x] T019b [US2] Cross-path eval on aeroStandard/longSequential: median 17.26m (wider but nominal for OOD paths)

**Checkpoint**: Wind variations working end-to-end. Fitness may be higher than P1 but must show improvement over generations.

---

## Phase 5: User Story 3 - Entry Condition Variations (Priority: P3)

**Goal**: Enable entry heading/roll/pitch/speed variations and validate GP evolves recovery behaviors

**Independent Test**: Enable EnableEntryVariations=1 with 36 scenarios. Run training, verify aircraft launches with varied attitudes.

### Implementation for User Story 3

- [ ] T020 [US3] Update autoc.ini: set EnableEntryVariations=1, WindScenarios=36 in `autoc/autoc.ini`
- [ ] T021 [US3] Run short validation (10 gens) to verify entry offsets applied: `./build/autoc 2>&1 | tee log/autoc-003-p3-smoke.log`
- [ ] T022 [US3] Verify autoc startup log shows variation table with non-zero entry offsets (heading, roll, pitch, speed)
- [ ] T023 [US3] If generation time is unacceptable, reduce WindScenarios to 25 or 16
- [ ] T024 [US3] Run 200-gen training: `./build/autoc 2>&1 | tee log/autoc-003-p3.log`
- [ ] T025 [US3] Analyze evolved GP tree from `data.stc` for recovery node usage (GETROLL_RAD, GETPITCH_RAD, GETVEL)
- [ ] T026 [US3] Render results to verify aircraft recovers from off-nominal entries
- [ ] T027 [US3] Document results: fitness value, crash rate, recovery behaviors observed

**Checkpoint**: Entry+wind variations working. Evolved controllers handle off-nominal entries.

---

## Phase 6: User Story 4 - Variable Rabbit Speed (Priority: P4)

**Goal**: Enable variable rabbit speed and validate GP evolves speed-adaptive throttle control

**Independent Test**: Set RabbitSpeedSigma=2.0, run training, verify speed variation in path timestamps.

### Implementation for User Story 4

- [ ] T028 [US4] Update autoc.ini: set RabbitSpeedSigma=2.0, VariationRampStep=0 (disable ramp for non-ramped baseline) in `autoc/autoc.ini`
- [ ] T029 [US4] Run 200-gen training with all variations: `./build/autoc 2>&1 | tee log/autoc-003-p4.log`
- [ ] T030 [US4] Analyze fitness — compare to P3 (entry+wind only)
- [ ] T031 [US4] Document results

**Checkpoint**: Full variation suite enabled (entry + wind + variable rabbit speed)

---

## Phase 7: User Story 5 - Progressive Variation Ramp (Priority: P5) — DEFERRED

**Goal**: Validate that VariationRampStep eases GP into harder conditions for better convergence

**Status**: DEFERRED — Ramp causes false divergence detection (elite re-evaluated under different conditions between generations). VariationRampStep=0 (disabled) produces clean training with zero divergences. Ramp requires divergence detector fix before re-enabling.

### Implementation for User Story 5

- [x] T032 [US5] Tested VariationRampStep=5: caused 2 false divergences at ramp step boundaries
- [x] T033 [US5] Root cause: `bitwiseEqual()` divergence check fails when ramp changes conditions between gens
- [ ] T034 [US5] Fix divergence detector to account for ramp-changed conditions (if ramp re-enabled)
- [ ] T035 [US5] Compare convergence curve to non-ramped full variations
- [ ] T036 [US5] Document whether ramp improves or hinders convergence

**Checkpoint**: Ramp deferred. VariationRampStep=0 is current production config.

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Final documentation, commit, and PR

- [ ] T037 [P] Update spec.md status from Draft to Complete in `specs/003-variations-redux/spec.md`
- [ ] T038 [P] Update BACKLOG.md: mark Variations Redux as DONE in `specs/BACKLOG.md`
- [ ] T039 Commit final autoc.ini with best-performing variation config in `autoc/autoc.ini`
- [ ] T040 Commit fitness constant changes in `autoc/autoc.h`
- [ ] T041 Create PR for GP repo via `gh pr create`
- [ ] T042 Commit any crrcsim changes (if logging adjustments made) and push branch

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Setup — BLOCKS all training phases
- **US1 (Phase 3)**: Depends on Foundational — baseline validation
- **US2 (Phase 4)**: Depends on US1 completion (need baseline to compare)
- **US3 (Phase 5)**: Depends on US2 completion (build on wind-only validation)
- **US4 (Phase 6)**: Depends on US3 completion (add speed on top of entry+wind)
- **US5 (Phase 7)**: Depends on US4 completion (compare ramp vs no-ramp)
- **Polish (Phase 8)**: Depends on all desired user stories

### User Story Dependencies

This feature has **sequential dependencies** between stories because each phase adds complexity on top of the previous:
- US1 (fitness tuning) → US2 (+ wind) → US3 (+ entry) → US4 (+ speed) → US5 (ramp tuning)

This is intentional — progressive validation prevents debugging a 36-scenario failure when the root cause might be fitness constants.

### Parallel Opportunities

- T001, T002, T003 (setup builds/tests) can run in parallel
- T004, T005 (fitness constant changes) can run in parallel (different lines, same file)
- T037, T038 (documentation updates) can run in parallel
- Within each training phase, smoke test and full run are sequential

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (build verification)
2. Complete Phase 2: Foundational (fitness constants)
3. Complete Phase 3: User Story 1 (tighter distance training)
4. **STOP and VALIDATE**: Is tracking distance < 10m? If not, adjust constants.
5. This alone delivers measurable improvement over 002 baseline.

### Incremental Delivery

1. US1 → Tighter tracking validated → Commit constants
2. US2 → Wind variations validated → Update autoc.ini
3. US3 → Entry variations validated → Full variation training
4. US4 → Variable speed added → Complete robustness suite
5. US5 → Ramp tuning → Final production config
6. Each phase builds confidence before adding complexity.

---

## Notes

- Each training run takes ~30-60 minutes (200 gens × 5K population × 36 scenarios)
- Log files in `log/` directory for comparison across phases
- If any phase shows regression, stop and diagnose before proceeding
- crrcsim has NO code changes — all variation application code already exists

### Code changes made during this feature (beyond config):
- `autoc/autoc.h`: DISTANCE_NORM 5.0→2.0, DISTANCE_POWER 1.5→2.0
- `autoc/gp_bytecode.h`: fitness_int (uint32_t) → fitness (double), version 1→2
- `autoc/gp_bytecode.cc`: fixed-point fitness → double, added std::fixed formatting
- `autoc/gpextractor.cc`: fitness storage as double, std::fixed formatting
- `autoc/bytecode2cpp.cc`: updated fitness display for double
- `src/gp.cc`: std::setprecision(17) for full double round-trip fidelity in GP serialization

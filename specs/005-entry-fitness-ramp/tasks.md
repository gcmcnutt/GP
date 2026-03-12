# Tasks: Entry Position Variations with Intercept-Budget Fitness Scaling

**Input**: Design documents from `/specs/005-entry-fitness-ramp/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Constants, data structures, and config plumbing shared across all user stories

- [x] T001 [P] Add intercept budget/scaling constants (INTERCEPT_SCALE_FLOOR, INTERCEPT_SCALE_CEILING, INTERCEPT_BUDGET_MAX, INTERCEPT_TURN_RATE) and entry safe-bounds constants (ENTRY_SAFE_RADIUS, ENTRY_SAFE_ALT_MIN, ENTRY_SAFE_ALT_MAX) to autoc/autoc.h alongside existing DISTANCE_NORM/DISTANCE_POWER constants
- [x] T002 [P] Add positionRadiusSigma and positionAltSigma fields to VariationSigmas struct in autoc/variation_generator.h
- [x] T003 [P] Add entryNorthOffset, entryEastOffset, entryAltOffset fields to VariationOffsets struct in autoc/variation_generator.h
- [x] T004 Add entryNorthOffset, entryEastOffset, entryAltOffset fields to ScenarioMetadata in autoc/minisim.h with version bump from 5 to 6 (backward-compatible: v5 loads set position offsets to 0.0)
- [x] T005 [P] Add entryPositionRadiusSigma (default 0.0) and entryPositionAltSigma (default 0.0) to ExtraConfig in autoc/autoc.h
- [x] T006 Add EntryPositionRadiusSigma and EntryPositionAltSigma config parameters to autoc/config_manager.cc
- [x] T007 [P] Add EntryPositionRadiusSigma and EntryPositionAltSigma (default 0) to autoc/autoc.ini and autoc/autoc-eval.ini in the SCENARIO VARIATIONS section, with comments following existing sigma param style
- [x] T008 [P] Remove stale autoc/autoc-profile.ini (legacy file with old defaults, no variation params)

**Checkpoint**: All data structures and constants in place. Build must pass with no functional change (all new fields default to 0).

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Minisim initial position cleanup and crrcsim position offset plumbing — required before position offsets or budget estimation can work correctly

**⚠️ CRITICAL**: US1 (scaling) can proceed without this phase, but US2 and US3 depend on correct position handling

- [x] T009 Ensure minisim (autoc/minisim.cc) compiles and runs cleanly with ScenarioMetadata v6 (new position offset fields). Minisim does NOT need to apply variation offsets — it just must not break on deserialization. Only crrcsim applies entry variations.
- [x] T010 [P] Add entryNorthOffset, entryEastOffset, entryAltOffset globals to ~/crsim/crrcsim-0.9.13/src/global.h and ~/crsim/crrcsim-0.9.13/src/global.cpp (default 0.0)
- [x] T011 [P] Read entryNorthOffset, entryEastOffset, entryAltOffset from ScenarioMetadata into globals in ~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp
- [x] T012 Apply position offsets (posX += entryNorthOffset, posY += entryEastOffset, Altitude += entryAltOffset) after VARIATIONS1 attitude block in ~/crsim/crrcsim-0.9.13/src/crrc_main.cpp

**Checkpoint**: crrcsim applies position offsets from ScenarioMetadata. Minisim compiles and runs cleanly with v6 metadata (does not apply variations). Build both repos (GP make, crrcsim make). Verify no behavioral change with offsets=0.

---

## Phase 3: User Story 1 — Intercept-Budget Distance Scaling (Priority: P1) 🎯 MVP

**Goal**: Per-step fitness scaling that ramps distance and attitude penalties from floor (~0.1) to ceiling (1.0) over the estimated intercept budget, giving GP gradient signal during approach phase.

**Independent Test**: Unit tests for computeInterceptScale() and computeInterceptBudget(). With budget=0 (no offset), scaling should be 1.0 from step 1 (backward compatible). With budget=7s, at t=1s scaling should be ~0.12.

### Tests for User Story 1

- [x] T013 [P] [US1] Add unit test for computeInterceptScale(): verify floor at t=0, quadratic ramp, ceiling at t=budget, ceiling beyond budget, degenerate budget=0 case in autoc/tests/gp_evaluator_tests.cc

### Implementation for User Story 1

- [x] T014 [US1] Implement computeInterceptScale(stepTime, budget) function in autoc/autoc.cc: `scale = FLOOR + (CEILING - FLOOR) * min(1, (t / budget))²`, handle budget≤0 → return CEILING
- [x] T015 [US1] Apply intercept scaling in GP-tree fitness loop (autoc/autoc.cc ~line 1192): compute interceptBudget from initial displacement+heading, then per-step `interceptScale * distance / DISTANCE_NORM` and `interceptScale * attitude_delta / ATTITUDE_NORM` before pow()
- [x] T016 [US1] Apply identical intercept scaling in bytecode fitness loop (autoc/autoc.cc ~line 1713) to maintain dual-mode parity
- [x] T017 [US1] Add unit test verifying backward compatibility: with zero position offset and zero heading offset, fitness values are identical to current (budget≈0 → immediate full penalty) in autoc/tests/gp_evaluator_tests.cc

**Checkpoint**: Intercept scaling works with existing attitude-only variations. Budget computed from displacement+heading. Full backward compatibility when no offsets applied.

---

## Phase 4: User Story 3 — Intercept Budget Estimation (Priority: P2)

**Goal**: Geometric time-to-intercept estimate from initial displacement, heading offset, aircraft speed, and rabbit speed. Crude "hacktor" — ±50% accuracy is fine.

**Independent Test**: Unit test with known inputs (60m offset, 90° heading, 20 m/s speed, 15 m/s rabbit) produces estimate in 5–10s range. Zero displacement → near-zero budget.

### Tests for User Story 3

- [x] T018 [P] [US3] Add unit test for computeInterceptBudget(): known-offset scenarios (60m/90°, 0m/0°, large offset capped at INTERCEPT_BUDGET_MAX) in autoc/tests/gp_evaluator_tests.cc

### Implementation for User Story 3

- [x] T019 [US3] Implement computeInterceptBudget(displacement, headingOffset, aircraftSpeed, rabbitSpeed) in autoc/autoc.cc: turn_time + closure_time + rabbit_compensation, clamped to INTERCEPT_BUDGET_MAX

**Checkpoint**: Budget estimation function tested independently. Feeds into US1 scaling (already wired in T015/T016).

---

## Phase 5: User Story 2 — Entry Position Variations (Priority: P2)

**Goal**: Cylindrical position generation (Gaussian radius + uniform angle → N/E, Gaussian altitude), clamped to safe arena bounds, flowing through the variation pipeline to simulators.

**Independent Test**: Run with EntryPositionRadiusSigma=30, EntryPositionAltSigma=10. Verify ScenarioMetadata carries non-zero position offsets, aircraft starts displaced, positions stay within safe bounds.

**Depends on**: Phase 2 (position offset plumbing in minisim/crrcsim)

### Tests for User Story 2

- [x] T020 [P] [US2] Add unit test for cylindrical position generation: verify Gaussian half-normal radius, uniform angle, Gaussian altitude offset, clamping to ENTRY_SAFE_RADIUS and ENTRY_SAFE_ALT bounds from EntryPositionRadiusSigma and EntryPositionAltSigma in autoc/tests/gp_evaluator_tests.cc

### Implementation for User Story 2

- [x] T021 [US2] Add cylindrical position generation to generateVariationsFromGPrand() in autoc/variation_generator.h: half-normal radius from positionRadiusSigma, uniform angle → N/E, Gaussian altitude from positionAltSigma, clamp to safe bounds
- [x] T022 [US2] Add position offset propagation in populateVariationOffsets() in autoc/autoc.cc: copy entryNorthOffset/entryEastOffset/entryAltOffset to ScenarioMetadata with RAMP_LANDSCAPE scaling. Note: RAMP_LANDSCAPE is the existing per-generation variation sigma ramp (computeVariationScale()) — distinct from the per-step intercept budget scaling in US1.
- [x] T023 [US2] Update VariationSigmas::fromDegrees() or add factory that accepts positionRadiusSigma and positionAltSigma in autoc/variation_generator.h
- [x] T024 [US2] Wire entryPositionRadiusSigma and entryPositionAltSigma from ExtraConfig into VariationSigmas at variation generation call site in autoc/autoc.cc

**Checkpoint**: Full position variation pipeline working. Aircraft starts at random positions within safe arena bounds. Fitness scaling from US1 handles the intercept phase.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Validation, build verification, integration testing

- [ ] T029 [US2] Add entry condition sanity clamping to variation_generator.h: constrain combined entry position and orientation to survivable envelopes. Eval run showed 5/49 crashes from extreme combos (e.g., -66° heading + -32° roll + 39m displaced, or 10m high + pitch-up + 36m displaced). Specific constraints:
  - Clamp pitch to ±15° (currently σ=7.5° but tails reach ±21°)
  - Clamp roll to ±45° (currently σ=22.5° but tails reach ±67°)
  - Reduce position altitude sigma or clamp altitude offset to ±8m (currently σ=10 allows ±17m tails, and starting low + pitch-down = ground impact)
  - Consider coupling: when heading offset > 45°, reduce position displacement proportionally (large heading + large displacement = unrecoverable)
  - Add compile-time constants for clamp bounds in autoc.h (ENTRY_MAX_PITCH, ENTRY_MAX_ROLL, ENTRY_MAX_ALT_OFFSET)
  - Apply clamps after Gaussian generation, before safe-bounds check

- [ ] T025 Verify build stability: GP `make` (from ~/GP/build), crrcsim `make` (from ~/crsim/crrcsim-0.9.13), xiao-gp `pio run` (from ~/xiao-gp) — MANUAL: user will verify crrcsim and xiao-gp builds
- [x] T026 Run all GoogleTest tests: `./build/autoc_tests` passes (101 tests)
- [ ] T027 Integration validation: run short evolution (5 generations) with EntryPositionRadiusSigma=30, EntryPositionAltSigma=10, and verify fitness differentiates between good and bad intercept behavior (SC-004: variance > 10% in gen 1) — MANUAL
- [ ] T028 Backward compatibility validation: run with EntryPositionRadiusSigma=0, EntryPositionAltSigma=0, and verify identical fitness to baseline (SC-005) — MANUAL

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (T004 for ScenarioMetadata fields)
- **US1 Scaling (Phase 3)**: Depends on Phase 1 (constants). Can proceed WITHOUT Phase 2 — uses existing displacement from initial state
- **US3 Budget (Phase 4)**: Depends on Phase 1. Independent of Phase 2
- **US2 Position (Phase 5)**: Depends on Phase 1 AND Phase 2 (position offset plumbing)
- **Polish (Phase 6)**: Depends on all previous phases

### User Story Dependencies

- **US1 (P1)**: Can start after Phase 1 — no dependency on other stories. **This is the MVP.**
- **US3 (P2)**: Can start after Phase 1 — the budget function is used by US1 but can be developed/tested independently
- **US2 (P2)**: Depends on Phase 2 (simulators must apply position offsets). Uses US1 scaling and US3 budget estimation

### Within Each User Story

- Tests written and verified to fail before implementation
- Implementation follows data flow: generation → propagation → application → fitness

### Parallel Opportunities

- T001, T002, T003, T005, T007, T008 can all run in parallel (different files/sections)
- T010, T011 can run in parallel (different crrcsim files)
- T013, T018, T020 (test stubs) can run in parallel
- US1 and US3 implementation can proceed in parallel after Phase 1

---

## Parallel Example: Phase 1 Setup

```bash
# Setup tasks touch different files/sections — run in parallel:
Task: "Add intercept + safe-bounds constants to autoc/autoc.h"        # T001
Task: "Add position sigma fields to VariationSigmas"                   # T002
Task: "Add position offset fields to VariationOffsets"                 # T003
Task: "Add position sigmas to ExtraConfig in autoc/autoc.h"           # T005
Task: "Add position sigma params to autoc.ini + autoc-eval.ini"       # T007
Task: "Remove stale autoc-profile.ini"                                 # T008
```

## Parallel Example: US1 + US3

```bash
# After Phase 1, US1 scaling and US3 budget can proceed simultaneously:
Task: "Implement computeInterceptScale() in autoc/autoc.cc"
Task: "Implement computeInterceptBudget() in autoc/autoc.cc"
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (constants, data structures, config)
2. Complete Phase 3: US1 — Intercept scaling with simple displacement-based budget
3. **STOP and VALIDATE**: Run tests, verify backward compatibility, check fitness differentiation
4. This alone provides the core value: GP gets gradient signal during intercept

### Incremental Delivery

1. Phase 1 → Setup complete
2. US1 (scaling) + US3 (budget) → Core intercept fitness working with existing attitude offsets
3. Phase 2 (position plumbing) + US2 (position generation) → Full entry position variations
4. Phase 6 → Integration validation, build stability confirmed
5. Each increment adds value without breaking previous functionality

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- US1 is independently valuable — scaling works with existing attitude offsets even without position offsets
- US3 (budget estimation) is a helper for US1 but can be developed/tested separately
- US2 (position generation) is the final piece that enables full entry condition training
- Constitution: testing-first, build stability (3 repos), dual-mode parity (GP tree + bytecode fitness loops)
- Commit after each task or logical group

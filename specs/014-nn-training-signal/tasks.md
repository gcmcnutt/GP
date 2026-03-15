# Tasks: NN Training Signal Improvement

**Input**: Design documents from `/specs/014-nn-training-signal/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md, contracts/

**Tests**: Contract tests are included per spec clarification ("contract tests first, before ripping things out").

**Organization**: Tasks follow the refactor-first implementation order from the spec clarifications. Foundational phases (GP removal, Boost removal, source reorg) are prerequisites. User story phases (sigma floor, curriculum, CMA-ES, etc.) follow.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Dependency analysis and contract tests that define surviving behavior

- [ ] T001 Run `nm ~/GP/lib/libgp.a | grep ' T '` and document all 103 exported symbols used by autoc in specs/014-nn-training-signal/analysis/libgp-symbols.md
- [ ] T002 Scan all `#include "gp.h"` and `#include "gpconfig.h"` across 11 autoc files, document transitive usage in specs/014-nn-training-signal/analysis/gp-includes.md
- [ ] T003 Scan all `#include <boost/` across autoc files, document per-file usage in specs/014-nn-training-signal/analysis/boost-includes.md
- [ ] T004 [P] Write contract test for NN evaluator (sensor-in/control-out: 22 inputs → 3 outputs) in autoc/tests/contract_evaluator_tests.cc
- [ ] T005 [P] Write contract test for config parsing (load autoc.ini, verify key types and defaults) in autoc/tests/contract_config_tests.cc
- [ ] T006 [P] Write contract test for NN evolution loop (1 generation: init → evaluate → select → reproduce → verify fitness improves) in autoc/tests/contract_evolution_tests.cc
- [ ] T006a [P] Write contract test for RPC transport (serialize EvalRequest/EvalResponse, round-trip, verify fields match) in autoc/tests/rpc_transport_tests.cc
- [ ] T007 Update autoc/CMakeLists.txt to add new contract test targets and ensure tests build into build/tests/
- [ ] T008 Verify all existing + new tests pass: `cd ~/GP/build && ctest --output-on-failure`

**Checkpoint**: Dependency surface fully documented, contract tests pass, safe to begin surgery

---

## Phase 2: Foundational — Strip GP Dependencies (FR-012, FR-014)

**Purpose**: Remove all GP code and libgp.a dependency. BLOCKS all subsequent work.

**⚠️ CRITICAL**: No user story work can begin until this phase is complete

### 2a: Vendor inih config parser

- [ ] T009 Vendor inih library (ini.h, ini.c, INIReader.h, INIReader.cpp) into autoc/third_party/inih/
- [ ] T010 Add inih to autoc/CMakeLists.txt as a source dependency
- [ ] T011 Rewrite autoc/config_manager.h to use inih INIReader instead of GPConfiguration/GPVariables/GPConfigVarInformation — remove `#include "gp.h"` and `#include "gpconfig.h"`
- [ ] T012 Rewrite autoc/config_manager.cc to use INIReader::GetInteger/GetReal/GetString with empty section — remove all DATAINT/DATADOUBLE/DATASTRING usage
- [ ] T013 Verify contract_config_tests pass with new parser

### 2b: Replace PRNG

- [ ] T014 [P] Create autoc/rng.h with std::mt19937 wrapper replacing GPrand()/GPsrand()/GPRandomPercent()
- [ ] T015 Replace GPrand() calls in autoc/pathgen.cc with rng.h wrapper — remove `#include "gp.h"`
- [ ] T016 Replace GPrand() calls in autoc/gp_math_utils.h with rng.h wrapper — remove `#include "gp.h"`

### 2c: Break GP inheritance and remove GP code

- [ ] T017 Rewrite autoc/autoc.h: remove `#include "gp.h"`, remove MyGP : GP / MyGene : GPGene / MyPopulation : GPPopulation inheritance, replace with standalone NN-only classes wrapping NNPopulation
- [ ] T018 Rewrite autoc/autoc.cc: remove GPInit(), GPRegisterClass(), GP tree creation, GPVariables usage — replace with direct NNPopulation init using config from inih. Remove `#include "gp.h"` and `#include "gpconfig.h"`
- [ ] T019 Delete autoc/autoc-eval.cc (GP tree evaluation — entirely replaced by nn_evaluator_portable)
- [ ] T020 [P] Delete autoc/gp_bytecode.h and autoc/gp_bytecode.cc (bytecode interpreter)
- [ ] T021 [P] Delete autoc/gp_evaluator_portable.h and autoc/gp_evaluator_portable.cc (GP portable evaluator)
- [ ] T022 Remove GP tree deserialization from autoc/minisim.cc — remove `#include "gp.h"`, keep only NN weight handling
- [ ] T023 Remove GP tree deserialization from autoc/renderer.cc — remove `#include "gp.h"`, keep only NN archive replay
- [ ] T024 Remove GPConfiguration usage from autoc/nnextractor.cc — use inih instead, remove `#include "gp.h"`
- [ ] T025 Remove autoc/gpextractor.cc from build (GP-only tool) — delete or move to archive
- [ ] T026 Strip all if(ControllerType=="GP") / if(ControllerType=="NN") modal branches across autoc.cc, minisim.cc — happy path NN only

### 2d: Sever CMake dependency on parent repo

- [ ] T027 Update autoc/CMakeLists.txt: remove `include_directories(../include)`, `link_directories(../lib)`, remove `gp` from all target_link_libraries
- [ ] T028 Remove references to deleted source files (autoc-eval.cc, gp_bytecode.cc, gp_evaluator_portable.cc, gpextractor.cc) from autoc/CMakeLists.txt
- [ ] T029 Full rebuild: `cd ~/GP/autoc && bash rebuild.sh` — verify zero GP references remain
- [ ] T030 Run all tests: `cd ~/GP/build && ctest --output-on-failure` — all contract + existing tests pass
- [ ] T030a Integration smoke test: run `./build/autoc` for 3 generations with NN config, verify data.dat and data.stc are produced and fitness values are present in console output

**Checkpoint**: autoc builds with zero dependency on libgp.a or parent GP repo. All `#include "gp.h"` removed.

---

## Phase 3: Foundational — Replace Boost (FR-013, phases C/D/E)

**Purpose**: Remove Boost dependency entirely

### 3a: Trivial Boost replacements

- [ ] T031 [P] Replace boost::thread with std::thread in autoc/threadpool.h
- [ ] T032 [P] Replace boost::mutex/unique_lock/condition_variable with std:: equivalents in autoc/threadpool.h and autoc/autoc.cc
- [ ] T033 [P] Replace boost::format with sprintf or std::stringstream in autoc/minisim.h and autoc/aircraft_state.h
- [ ] T034 [P] Replace boost::date_time with std::chrono in autoc/autoc.cc logging setup
- [ ] T035 Replace boost::process with fork()/execv() in autoc/threadpool.h

### 3b: Medium Boost replacements

- [ ] T036 Create autoc/socket_wrapper.h with POSIX TCP socket wrapper (~150 LOC) replacing boost::asio usage
- [ ] T037 Replace boost::asio in autoc/threadpool.h (tcp::acceptor, io_service) with socket_wrapper.h
- [ ] T038 Replace boost::asio in autoc/minisim.cc (tcp::resolver, connect) with socket_wrapper.h
- [ ] T039 Replace boost::asio in autoc/minisim.h (read/write) with socket_wrapper.h
- [ ] T040 Rewrite autoc/logger.h and autoc/logger.cc: replace boost::log with simple custom logger (~30 LOC, severity levels, file + stderr output)
- [ ] T041 Replace boost::iostreams with std::stringstream in autoc/minisim.h and autoc/autoc.cc

### 3c: Boost serialization replacement + unified transport

- [ ] T042 Define EvalRequest and EvalResponse structs with manual serialize/deserialize methods per contracts/rpc-transport.md in autoc/rpc_protocol.h
- [ ] T043 Rewrite sendRPC/receiveRPC in autoc/minisim.h to use manual binary serialization instead of boost::archive
- [ ] T044 Remove BOOST_SERIALIZATION_ACCESS, BOOST_CLASS_VERSION, and all boost::serialization friend classes from autoc/minisim.h, autoc/aircraft_state.h, autoc/gp_bytecode.h (if still present)
- [ ] T045 Remove all `#include <boost/` directives from all autoc source files
- [ ] T046 Remove Boost from find_package and target_link_libraries in autoc/CMakeLists.txt
- [ ] T047 Full rebuild: `cd ~/GP/autoc && bash rebuild.sh` — verify zero Boost references remain
- [ ] T048 Run all tests: `cd ~/GP/build && ctest --output-on-failure`
- [ ] T048a Integration smoke test: run `./build/autoc` for 3 generations, verify RPC communication works (autoc ↔ minisim), data.dat produced, fitness values in output

**Checkpoint**: autoc builds with zero Boost dependency. RPC uses plain binary protocol.

---

## Phase 4: Foundational — Source Reorganization (Phase F)

**Purpose**: Clean C++ project layout enabling standalone repo extraction

- [ ] T049 Create directory structure: autoc/include/autoc/{nn,eval,util}/, autoc/src/{nn,eval,util}/, autoc/tools/
- [ ] T050 [P] `git mv` NN headers to include/autoc/nn/: nn_topology.h, nn_genome.h (from nn_population.h), nn_population.h, nn_forward.h, nn_serialization.h
- [ ] T051 [P] `git mv` eval headers to include/autoc/eval/: nn_evaluator_portable.h, aircraft_state.h
- [ ] T052 [P] `git mv` util headers to include/autoc/util/: config_manager.h (→ config.h), logger.h, s3_archive.h, socket_wrapper.h, rng.h
- [ ] T053 [P] `git mv` NN source files to src/nn/
- [ ] T054 [P] `git mv` eval source files to src/eval/
- [ ] T055 [P] `git mv` util source files to src/util/
- [ ] T056 [P] `git mv` tool executables to tools/: nnextractor.cc, nn2cpp.cc, renderer.cc
- [ ] T057 [P] `git mv` test files to tests/
- [ ] T058 Update all `#include` directives across all moved files to use `"autoc/nn/..."`, `"autoc/eval/..."`, `"autoc/util/..."` style
- [ ] T059 Rewrite autoc/CMakeLists.txt for new directory structure: include paths, source file lists, add_subdirectory(tests), test output to build/tests/
- [ ] T060 Full rebuild and test: `cd ~/GP/autoc && bash rebuild.sh && cd ~/GP/build && ctest --output-on-failure`
- [ ] T060a Integration smoke test: run `./build/autoc` for 3 generations, verify everything still works after file moves (same output as T048a)

**Checkpoint**: Clean source layout, all tests pass, ready for cross-repo updates

---

## Phase 5: Foundational — Cross-Repo & Output Cleanup (Phases G/H)

**Purpose**: Update CRRCSim and xiao-gp, clean up run output

- [ ] T061 Update CRRCSim minisim RPC to match new binary transport protocol in ~/crsim/crrcsim-0.9.13 (matching autoc/rpc_protocol.h)
- [ ] T062 Remove `#ifdef GP_BUILD` guards from shared NN evaluator code in ~/crsim/crrcsim-0.9.13
- [ ] T063 Build CRRCSim: `cd ~/crsim/crrcsim-0.9.13/build && cmake .. && make`
- [ ] T064 Verify xiao-gp export pipeline: `./build/nnextractor` → `./build/nn2cpp` → `cd ~/xiao-gp && pio run`
- [ ] T065 Remove `#ifdef GP_BUILD` guards from shared NN evaluator code in ~/xiao-gp
- [ ] T066 Add OutputDir config key to autoc.ini (default: current directory for backward compat)
- [ ] T067 Implement auto-created run subdirectory in autoc/src/autoc.cc: create OutputDir/{timestamp}/ at startup, route data.dat, data.stc, logs to it
- [ ] T068 Remove hardcoded `eval-` prefix from data file naming in autoc/src/autoc.cc
- [ ] T069 Run all 3 repo builds and verify: GP autoc, CRRCSim, xiao-gp
- [ ] T069a Integration smoke test: run `./build/autoc` for 3 generations with OutputDir set, verify run artifacts land in auto-created subdirectory
- [ ] T069b Cross-repo integration test: run autoc with minisim worker from CRRCSim build, verify RPC works end-to-end

**Checkpoint**: All 3 repos build, cross-repo contracts aligned, run output goes to subdirectories

---

## Phase 6: User Story 1 — Sigma Floor (Priority: P1) 🎯 MVP

**Goal**: Prevent search freeze by clamping mutation sigma to a configurable minimum

**Independent Test**: Run short evolution, verify sigma never drops below floor, fitness continues improving

### Contract Tests

- [ ] T070 [US1] Write test: sigma clamped to floor when it would decay below in autoc/tests/sigma_floor_tests.cc
- [ ] T071 [US1] Write test: sigma floor = 0 produces identical behavior to unclamped in autoc/tests/sigma_floor_tests.cc

### Implementation

- [ ] T072 [US1] Add NNSigmaFloor config key (default: 0) to config parser in autoc/src/util/config.cc
- [ ] T073 [US1] Implement sigma floor clamping in NNPopulation::mutate() in autoc/src/nn/nn_population.cc — clamp sigma after self-adaptive update, log clamping event
- [ ] T074 [US1] Add sigma floor edge case: warn if floor > initial sigma in autoc/src/nn/nn_population.cc
- [ ] T075 [US1] Run tests: `cd ~/GP/build && ctest -R sigma_floor --output-on-failure`
- [ ] T075a [US1] Integration test: run `./build/autoc` for 20 generations with NNSigmaFloor=0.05, verify sigma never drops below 0.05 in console output

**Checkpoint**: Sigma floor working. SC-001 can be validated with a long evolution run.

---

## Phase 7: User Story 2 — Curriculum Scenario Ramp (Priority: P1)

**Goal**: Progressive scenario difficulty ramping for NN evolution

**Independent Test**: Run NN evolution, observe early gens use fewer scenarios, later gens use full suite

### Contract Tests

- [ ] T076 [US2] Write test: CurriculumSchedule parses "1:50,7:150,49:0" into 3 stages in autoc/tests/curriculum_tests.cc
- [ ] T077 [US2] Write test: stage transitions occur at correct generation boundaries in autoc/tests/curriculum_tests.cc
- [ ] T078 [US2] Write test: disabled curriculum uses all scenarios from gen 0 in autoc/tests/curriculum_tests.cc

### Implementation

- [ ] T079 [US2] Create CurriculumSchedule class in autoc/include/autoc/eval/curriculum.h and autoc/src/eval/curriculum.cc
- [ ] T080 [US2] Add CurriculumEnabled and CurriculumSchedule config keys to config parser in autoc/src/util/config.cc
- [ ] T081 [US2] Wire CurriculumSchedule into NN evolution loop in autoc/src/autoc.cc — query active scenario count per generation, log stage transitions
- [ ] T082 [US2] Run tests: `cd ~/GP/build && ctest -R curriculum --output-on-failure`
- [ ] T082a [US2] Integration test: run `./build/autoc` for 20 generations with CurriculumSchedule=1:5,7:15,49:0, verify stage transitions in console output

**Checkpoint**: Curriculum ramping working. SC-002 can be validated with comparative evolution runs.

---

## Phase 8: User Story 3 — Per-Scenario Fitness Decomposition (Priority: P2)

**Goal**: Minimax or percentile fitness aggregation across scenarios

**Independent Test**: Compare population rankings under sum vs minimax, verify minimax penalizes high-variance individuals

### Contract Tests

- [ ] T083 [US3] Write test: minimax returns worst-case scenario score in autoc/tests/fitness_aggregator_tests.cc
- [ ] T084 [US3] Write test: percentile returns correct value at 95th percentile in autoc/tests/fitness_aggregator_tests.cc
- [ ] T085 [US3] Write test: sum mode matches current behavior in autoc/tests/fitness_aggregator_tests.cc

### Implementation

- [ ] T086 [US3] Create FitnessAggregator class in autoc/include/autoc/eval/fitness_aggregator.h and autoc/src/eval/fitness_aggregator.cc
- [ ] T087 [US3] Add FitnessAggregation and FitnessPercentile config keys to config parser in autoc/src/util/config.cc
- [ ] T088 [US3] Wire FitnessAggregator into fitness computation in autoc/src/autoc.cc — replace sum-over-scenarios with aggregator call
- [ ] T089 [US3] Run tests: `cd ~/GP/build && ctest -R fitness_aggregator --output-on-failure`
- [ ] T089a [US3] Integration test: run `./build/autoc` for 20 generations with FitnessAggregation=minimax, verify worst-case scenario fitness drives selection in console output

**Checkpoint**: Fitness aggregation working. SC-003 can be validated with comparative runs.

---

## Phase 9: User Story 4 — sep-CMA-ES Optimizer (Priority: P2)

**Goal**: Replace GA with sep-CMA-ES for efficient 531-dimensional optimization

**Independent Test**: Run sep-CMA-ES on Rosenbrock benchmark, then on NN flight controller, compare against GA baseline

### Contract Tests

- [ ] T090 [US4] Write test: CMA-ES converges on Rosenbrock function (N=10) in autoc/tests/cmaes_tests.cc
- [ ] T091 [US4] Write test: ask() generates lambda candidates from distribution in autoc/tests/cmaes_tests.cc
- [ ] T092 [US4] Write test: tell() updates mean/sigma/covariance correctly in autoc/tests/cmaes_tests.cc
- [ ] T093 [US4] Write test: CMA-ES state serialization round-trips in autoc/tests/cmaes_tests.cc

### Implementation

- [ ] T094 [US4] Implement SepCMAES class in autoc/include/autoc/nn/sep_cmaes.h and autoc/src/nn/sep_cmaes.cc — ask/tell interface, Eigen vectors, hyperparams from research.md (lambda=50, mu=25, sep-scaled learning rates)
- [ ] T095 [US4] Add OptimizerType config key ("ga" or "sep-cma-es") to config parser in autoc/src/util/config.cc
- [ ] T096 [US4] Wire SepCMAES into evolution loop in autoc/src/autoc.cc — replace mutate/crossover/select with ask/tell when OptimizerType=sep-cma-es
- [ ] T097 [US4] Add CMA-ES logging: sigma, best fitness, mean fitness per generation in autoc/src/autoc.cc
- [ ] T098 [US4] Add CMA-ES state to S3 archive format in autoc/src/nn/nn_serialization.cc
- [ ] T099 [US4] Run tests: `cd ~/GP/build && ctest -R cmaes --output-on-failure`
- [ ] T099a [US4] Integration test: run `./build/autoc` for 20 generations with OptimizerType=sep-cma-es and PopulationSize=50, verify sigma/fitness logging and fitness improvement

**Checkpoint**: sep-CMA-ES working. SC-004 can be validated with comparative evolution runs.

---

## Phase 10: User Story 7 — Per-Timestep Fitness Streaming (Priority: P3)

**Goal**: Return per-timestep data from simulator instead of just aggregate scalar

**Independent Test**: Run single evaluation, verify returned data includes per-timestep distance/attitude/commands

**Note**: Moved before US5/US6 because both depend on per-timestep data

### Implementation

- [ ] T100 [US7] Define TimestepRecord struct in autoc/include/autoc/eval/timestep_record.h
- [ ] T101 [US7] Add per-timestep data collection in minisim evaluation loop in autoc/src/minisim.cc — populate TimestepRecord vector during simulation
- [ ] T102 [US7] Extend EvalResponse in autoc/rpc_protocol.h to include per-timestep data (timestep_count + records per scenario)
- [ ] T103 [US7] Update RPC serialization in autoc/src/minisim.cc to send per-timestep data when enabled
- [ ] T104 [US7] Update CRRCSim minisim to match new response format in ~/crsim/crrcsim-0.9.13
- [ ] T105 [US7] Verify aggregate fitness computed from streamed data matches legacy scalar in autoc/tests/timestep_streaming_tests.cc
- [ ] T106 [US7] Run tests: `cd ~/GP/build && ctest -R timestep --output-on-failure`

**Checkpoint**: Per-timestep data flowing. Enables US5 (segment scoring) and US6 (behavioral cloning data).

---

## Phase 11: User Story 5 — Per-Segment Credit Assignment (Priority: P3)

**Goal**: Score trajectory segments by error reduction to amplify fitness signal

**Independent Test**: Compute segment scores on recorded trajectory, verify controllers with good local corrections score higher

**Depends on**: US7 (per-timestep streaming)

### Implementation

- [ ] T107 [US5] Create TrajectorySegment struct and segment scoring function in autoc/include/autoc/eval/segment_scorer.h and autoc/src/eval/segment_scorer.cc
- [ ] T108 [US5] Implement difficulty weighting (turn rate, crosswind) in segment_scorer.cc
- [ ] T109 [US5] Wire segment scoring into fitness aggregation pipeline in autoc/src/autoc.cc — use TimestepRecord data from US7
- [ ] T110 [US5] Write tests: segment scoring produces expected scores for known trajectories in autoc/tests/segment_scorer_tests.cc
- [ ] T111 [US5] Run tests: `cd ~/GP/build && ctest -R segment --output-on-failure`

**Checkpoint**: Segment scoring working. Richer fitness signal for evolution.

---

## Phase 12: User Story 6 — Behavioral Cloning Bootstrap (Priority: P3) — DEFERRED

**Status**: Deferred to backlog. Only pursue if direct NN training is exhausted and GP→NN weight transfer is attempted.

~~- [ ] T112-T117: Behavioral cloning tasks~~

**Checkpoint**: Skipped — proceed to Phase 13.

---

## Phase 13: User Story 8 — Checkpoint/Resume (Priority: P3)

**Goal**: Save full evolution state per generation for crash recovery and resume

**Independent Test**: Run 10 gens, kill, resume, verify gen 11 continues correctly

### Implementation

- [ ] T118 [US8] Define EvolutionCheckpoint struct with serialization in autoc/include/autoc/nn/checkpoint.h and autoc/src/nn/checkpoint.cc
- [ ] T119 [US8] Include CMAESState in checkpoint when optimizer_type=sep-cma-es (requires US4/Phase 9 complete)
- [ ] T120 [US8] Implement checkpoint writing at end of each generation in autoc/src/autoc.cc
- [ ] T121 [US8] Implement checkpoint loading and resume at startup in autoc/src/autoc.cc — detect existing checkpoint, resume from saved generation
- [ ] T122 [US8] Add checkpoint corruption detection (magic bytes + checksum) in autoc/src/nn/checkpoint.cc
- [ ] T123 [US8] Write test: checkpoint round-trip produces identical state in autoc/tests/checkpoint_tests.cc
- [ ] T124 [US8] Run tests: `cd ~/GP/build && ctest -R checkpoint --output-on-failure`
- [ ] T124a [US8] Integration test: run `./build/autoc` for 5 generations, kill, resume, verify generation 6 continues with correct population state

**Checkpoint**: Checkpoint/resume working. SC-006 can be validated by comparing interrupted vs uninterrupted runs.

---

## Phase 14: Polish & Cross-Cutting Concerns

**Purpose**: Repo extraction preparation and final cleanup

- [ ] T125 Verify all 3 repo builds pass: GP autoc (`rebuild.sh`), CRRCSim (`make`), xiao-gp (`pio run`)
- [ ] T126 Run full test suite: `cd ~/GP/build && ctest --output-on-failure`
- [ ] T127 Verify export pipeline end-to-end: nnextractor → nn2cpp → xiao-gp build
- [ ] T128 Update autoc.ini with all new config keys and sensible defaults
- [ ] T129 Update constitution build commands to reflect standalone autoc
- [ ] T130 Prepare autoc/ for standalone repo extraction: verify CMakeLists.txt has no parent references, all includes use autoc/ prefix
- [ ] T131 Run quickstart.md validation: follow all steps and verify they work

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Strip GP)**: Depends on Phase 1 contract tests — BLOCKS all subsequent work
- **Phase 3 (Replace Boost)**: Depends on Phase 2 — can overlap with late Phase 2 work
- **Phase 4 (Source Reorg)**: Depends on Phases 2+3 — all code changes done before moving files
- **Phase 5 (Cross-Repo)**: Depends on Phase 4 — transport contract finalized
- **Phases 6-13 (User Stories)**: Depend on Phase 5 — can proceed in priority order or parallel

### User Story Dependencies

- **US1 (Sigma Floor, P1)**: Independent — can start after Phase 5
- **US2 (Curriculum, P1)**: Independent — can start after Phase 5
- **US3 (Fitness Decomposition, P2)**: Independent — can start after Phase 5
- **US4 (sep-CMA-ES, P2)**: Independent — can start after Phase 5
- **US7 (Per-Timestep Streaming, P3)**: Independent — moved early because US5 and US6 depend on it
- **US5 (Segment Scoring, P3)**: Depends on US7 (needs per-timestep data)
- **US6 (Behavioral Cloning, P3)**: Depends on US7 (for data collection format)
- **US8 (Checkpoint/Resume, P3)**: T119 depends on US4 (CMAESState must exist). Rest of US8 is independent.

### Within Each User Story

- Contract tests first → implementation → integration → verify

### Parallel Opportunities

- T004/T005/T006: All contract tests in Phase 1
- T031/T032/T033/T034: All trivial Boost replacements
- T050/T051/T052/T053/T054/T055/T056/T057: All git mv operations
- US1 and US2 can run in parallel (both P1, independent)
- US3 and US4 can run in parallel (both P2, independent)

---

## Parallel Example: Phase 2c (GP removal)

```bash
# These can run in parallel (different files, no dependencies):
Task T019: "Delete autoc/autoc-eval.cc"
Task T020: "Delete autoc/gp_bytecode.h and autoc/gp_bytecode.cc"
Task T021: "Delete autoc/gp_evaluator_portable.h and autoc/gp_evaluator_portable.cc"
```

## Parallel Example: User Stories 1 + 2 (both P1)

```bash
# After Phase 5, these can run in parallel:
Task T072-T075: "Sigma Floor implementation"
Task T079-T082: "Curriculum Ramp implementation"
```

---

## Implementation Strategy

### MVP First (Phases 1-5 + User Story 1)

1. Complete Phase 1: Setup + contract tests
2. Complete Phase 2: Strip GP (CRITICAL — largest refactor)
3. Complete Phase 3: Replace Boost
4. Complete Phase 4: Source reorganization
5. Complete Phase 5: Cross-repo updates
6. Complete Phase 6: Sigma Floor (US1)
7. **STOP and VALIDATE**: Run evolution with sigma floor, verify SC-001

### Incremental Delivery

1. Phases 1-5 → Standalone NN-only autoc (major milestone)
2. Add US1 (Sigma Floor) → Immediate fix for nn13 stall
3. Add US2 (Curriculum) → Progressive difficulty
4. Add US3 (Fitness Decomposition) → Robust controllers
5. Add US4 (sep-CMA-ES) → Efficient optimization
6. Add US7 → US5 → US6 → Advanced scoring + warm start
7. Add US8 (Checkpoint) → Crash recovery

### Repo Extraction (after Phase 14)

Separate follow-up: extract autoc/ to standalone repo, migrate relevant specs/docs.

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- All file moves use `git mv` to preserve history
- Happy path only — strip all defensive fallback branches
- Contract tests define surviving behavior — don't break them during refactoring
- Commit after each logical group of tasks
- Stop at any checkpoint to validate independently

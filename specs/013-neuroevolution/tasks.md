# Tasks: 013-neuroevolution

**Input**: Design documents from `/specs/013-neuroevolution/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md

**Tests**: Required by project constitution (Testing-First principle). Tests written before implementation.

**Organization**: Tasks grouped by spec implementation phase. Each phase is an independently testable increment.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which implementation phase (US1=Unify Eval, US2=NN Core, US3=Serialization, US4=Evolution, US5=Embedded Codegen, US6=Renderer/S3)

---

## Phase 1: Setup (Build System & File Scaffolding)

**Purpose**: Create new source files and update CMakeLists.txt for NN support

- [x] T001 Create stub header autoc/eval_backend.h with ControllerBackend interface declaration (virtual evaluate, virtual getName)
- [x] T002 [P] Create stub files autoc/fitness_computer.h and autoc/fitness_computer.cc
- [x] T003 [P] Create stub files autoc/eval_logger.h and autoc/eval_logger.cc
- [x] T004 [P] Create stub files autoc/nn_evaluator_portable.h and autoc/nn_evaluator_portable.cc with normalization constants (NORM_ANGLE, NORM_DIST, NORM_VEL, NORM_RATE)
- [x] T005 [P] Create stub files autoc/nn_serialization.h and autoc/nn_serialization.cc
- [x] T006 [P] Create stub files autoc/nn_population.h and autoc/nn_population.cc
- [x] T007 Add all new source files to autoc/CMakeLists.txt: nn_evaluator_portable, fitness_computer, eval_logger, nn_serialization, nn_population
- [x] T008 Add test targets to autoc/CMakeLists.txt: nn_evaluator_tests.cc, nn_serialization_tests.cc, nn_population_tests.cc, fitness_computer_tests.cc in autoc/tests/
- [x] T009 Verify build compiles with stubs (full rebuild required — new files added to CMakeLists.txt): `cd ~/GP/autoc && bash rebuild.sh`

---

## Phase 2: Unify Evaluation Pipeline — US1 (Absorbs 007-unify-eval)

**Goal**: Extract ~350 lines of duplicated eval logic into shared components. Define ControllerBackend interface so GP tree, bytecode, and NN all use the same fitness/logging/scenario paths.

**Independent Test**: Both GP tree mode (EvaluateMode=0) and bytecode mode (EvaluateMode=1) produce identical fitness values after refactor.

**Functional Requirements**: FR-008b, FR-012

### Tests for US1

> **Write tests FIRST, ensure they FAIL before implementation**

- [x] T010 [US1] Write parity test in autoc/tests/fitness_computer_tests.cc: verify FitnessComputer::computeStepPenalty matches existing inline fitness calc for known distance/attitude values
- [x] T011 [P] [US1] Write parity test: verify FitnessComputer::computeCrashPenalty matches existing crash penalty for known fraction_completed values
- [x] T012 [P] [US1] Write parity test: verify FitnessComputer::computeAttitudeScale matches existing attitude scale computation

### Implementation for US1

- [x] T013 [US1] Define ControllerBackend interface in autoc/eval_backend.h: virtual evaluate(AircraftState&, PathProvider&) and virtual getName() methods
- [x] T014 [US1] Add ControllerType enum to EvalData in autoc/minisim.h: GP_TREE=0, BYTECODE=1, NEURAL_NET=2 — explicit type tag replaces magic byte heuristics in minisim
- [x] T015 [US1] Implement FitnessComputer in autoc/fitness_computer.cc: extract computeStepPenalty, computeCrashPenalty, computeAttitudeScale from autoc.cc lines ~1278-1350
- [x] T016 [US1] Implement EvalLogger in autoc/eval_logger.cc: stubs in place, actual extraction deferred (NN will have own logging)
- [ ] ~~T017 [US1] Extract scenario metadata builder~~ — deferred (too many globals, NN will have own eval path)
- [x] T018 [US1] Extract GPrandGaussian() utility from variation_generator.h: Box-Muller transform over GPrand(), replaces 4 duplicate inline lambdas. Shared by variation_generator.h and nn_population.cc. Place in a shared header (e.g. autoc/gp_math_utils.h).
- [ ] ~~T019 [US1] Implement GPTreeBackend wrapper~~ — deferred (NN will have own eval path, not needed for critical path)
- [ ] ~~T020 [P] [US1] Implement BytecodeBackend wrapper~~ — deferred (NN will have own eval path, not needed for critical path)
- [ ] ~~T021 [US1] Merge evalTask() into single implementation~~ — deferred (too risky, NN will have own evalTask)
- [x] T022 [US1] Update minisim.cc format detection: use ControllerType tag from EvalData instead of magic byte heuristics. ControllerType switch dispatch added.
- [x] T023 [US1] Verify parity: run GP tree mode and bytecode mode with same inputs, confirm identical fitness values — verified
- [x] T024 [US1] Verify build and tests: `cd ~/GP && make && cd build && ctest --output-on-failure`

**Checkpoint**: Evaluation pipeline unified. Single evalTask() with shared scenario builder, FitnessComputer, EvalLogger, and ControllerBackend interface. BytecodeEvaluationGP inline class eliminated. Minisim uses explicit type tags. All existing tests pass. Ready for NN backend.

---

## Phase 3: NN Evaluator Core — US2

**Goal**: Implement portable NN forward pass with tanh LUT. Create NNControllerBackend conforming to Phase 2 interface.

**Independent Test**: Forward pass produces correct outputs for known weight matrices. Output values always in [-1, 1]. Deterministic for same inputs.

**Functional Requirements**: FR-001, FR-004, FR-005, FR-011, FR-014

### Tests for US2

> **Write tests FIRST, ensure they FAIL before implementation**

- [x] T025 [P] [US2] Write forward pass test in autoc/tests/nn_evaluator_tests.cc: single-layer identity weights produce expected output
- [x] T026 [P] [US2] Write forward pass test: multi-layer 14-16-8-3 topology with known weights produces expected outputs (hand-computed reference)
- [x] T027 [P] [US2] Write output range test: random weights always produce outputs in [-1, 1] via tanh
- [x] T028 [P] [US2] Write tanh LUT test: fast_tanh() matches std::tanh() within 1e-2 tolerance across [-5, 5]
- [x] T029 [P] [US2] Write Xavier init test: initialized weights have zero mean and expected variance (1/fan_in) within statistical bounds
- [x] T030 [P] [US2] Write determinism test: same weights + same inputs produce bit-exact same outputs across repeated calls
- [x] T031 [P] [US2] Write weight count test: topology {14,16,8,3} produces exactly 403 weights
- [x] T032 [P] [US2] Write input normalization test: verify sensor values are divided by correct NORM_* constants before forward pass

### Implementation for US2

- [x] T033 [US2] Define NNGenome struct in autoc/nn_evaluator_portable.h: weights vector, topology vector, fitness, generation, mutation_sigma per data-model.md
- [x] T034 [US2] Implement nn_forward() in autoc/nn_evaluator_portable.cc: feedforward pass with row-major layer-sequential weight layout, tanh activation
- [x] T035 [US2] Implement fast_tanh() LUT in autoc/nn_evaluator_portable.cc: 512-entry table, domain [-5, 5], linear interpolation (same pattern as existing sin LUT in gp_evaluator_portable.cc)
- [x] T036 [US2] Implement nn_weight_count() utility: compute total weights+biases from topology vector
- [x] T037 [US2] Implement nn_xavier_init() in autoc/nn_evaluator_portable.cc: Xavier/Glorot initialization per topology, using local Gaussian in test/embedded builds, GPrandGaussian() in full GP build
- [x] T038 [US2] Implement nn_gather_inputs() in autoc/nn_evaluator_portable.cc: call 14 sensor functions (executeGetDPhi, executeGetDist, etc.), apply NORM_ANGLE/NORM_DIST/NORM_VEL/NORM_RATE normalization, build input vector
- [x] T039 [US2] Implement NNControllerBackend in autoc/nn_evaluator_portable.cc: wraps nn_gather_inputs + nn_forward + setPitch/setRoll/setThrottle, conforms to ControllerBackend interface from eval_backend.h
- [x] T040 [US2] Verify build and tests: `cd ~/GP && make && cd build && ctest --output-on-failure`

**Checkpoint**: NN forward pass works, tanh LUT accurate, NNControllerBackend plugs into unified eval pipeline. Can manually test NN eval with hardcoded weights.

---

## Phase 4: NN Serialization & Archive — US3

**Goal**: Binary serialization format for NN genomes (RPC transport, S3 storage, minisim detection). Incompatible with GP/bytecode formats by design. NN runs use `nn-{timestamp}/` S3 prefix.

**Independent Test**: Round-trip serialize/deserialize preserves all NNGenome fields exactly. Minisim detects NN format by magic number.

**Functional Requirements**: FR-006, FR-007, FR-008a, FR-017

### Tests for US3

> **Write tests FIRST, ensure they FAIL before implementation**

- [x] T041 [P] [US3] Write round-trip test in autoc/tests/nn_serialization_tests.cc: serialize NNGenome → deserialize → all fields match exactly (weights, topology, fitness, generation)
- [x] T042 [P] [US3] Write magic number test: serialized data starts with "NN01" magic bytes
- [x] T043 [P] [US3] Write format detection test: NN data detected as NN, GP tree data not detected as NN, bytecode data not detected as NN
- [x] T044 [P] [US3] Write corrupt data test: truncated/invalid data returns error, does not crash

### Implementation for US3

- [x] T045 [US3] Implement nn_serialize() in autoc/nn_serialization.cc: write NNGenome to binary format per data-model NNSerializationFormat (magic "NN01", topology, weights, metadata)
- [x] T046 [US3] Implement nn_deserialize() in autoc/nn_serialization.cc: read binary format back to NNGenome with validation
- [x] T047 [US3] Implement nn_detect_format() in autoc/nn_serialization.cc: check first 4 bytes for "NN01" magic
- [x] T048 [US3] Boost adapter not needed — EvalData.gp (vector<char>) carries serialized NN blob natively. autoc fills it via nn_serialize, minisim reads via nn_deserialize.
- [x] T049 [US3] Integrate NN format detection in autoc/minisim.cc: ControllerType::NEURAL_NET case deserializes NN genome and creates NNControllerBackend
- [x] T050 [US3] Implement NN evaluation in autoc/minisim.cc: NNControllerBackend wraps nn_gather_inputs + nn_forward + setPitch/setRoll/setThrottle
- [x] T051 [US3] Update S3 archive key generation in autoc/autoc.cc: generate_iso8601_timestamp() uses "nn-" prefix when ControllerType=NN
- [x] T052 [US3] Implement nnextractor tool — extracts best NN genome from S3 EvalResults archives, writes NN01 binary weight file. CLI: -k keyname, -g gen, -o output, -i config
- [x] T053 [US3] Add nnextractor build target in autoc/CMakeLists.txt
- [x] T054 [US3] Verify build and tests: `cd ~/GP && make` — all 162 tests pass

**Checkpoint**: NN genomes can be serialized/deserialized for RPC and S3. Minisim auto-detects NN format. S3 archives use `nn-` prefix. nnextractor extracts best weights.

---

## Phase 5: Evolution Integration — US4

**Goal**: Wire NNPopulation into existing evolution engine. Arithmetic crossover, Gaussian mutation, tournament selection on weight vectors. Full ControllerType config with all 6 NN-specific params.

**Independent Test**: Run NN evolution for 5 generations on a trivial problem (constant target) and verify fitness improves.

**Functional Requirements**: FR-002, FR-003, FR-009, FR-010, FR-014, FR-015

### Tests for US4

> **Write tests FIRST, ensure they FAIL before implementation**

- [x] T055 [P] [US4] Write crossover test in autoc/tests/nn_population_tests.cc: arithmetic crossover of two known genomes produces expected blended weights
- [x] T056 [P] [US4] Write mutation test: Gaussian mutation changes weights, mean change is near zero over many trials, sigma controls spread
- [x] T057 [P] [US4] Write self-adaptive sigma test: mutation_sigma itself mutates and stays positive
- [x] T058 [P] [US4] Write population init test: all individuals have correct weight count, all weights finite, no NaN/Inf
- [x] T059 [P] [US4] Write tournament selection test: higher fitness individuals selected more often over many tournaments

### Implementation for US4

- [x] T060 [US4] Implement NNPopulation struct in autoc/nn_population.h: vector of NNGenome, shared topology, generation counter, best_fitness/best_index per data-model.md
- [x] T061 [US4] Implement nn_arithmetic_crossover() in autoc/nn_population.cc: BLX-alpha blend of two parent weight vectors, alpha from NNCrossoverAlpha config (-1 = uniform [0,1] via GPrand())
- [x] T062 [US4] Implement nn_gaussian_mutation() in autoc/nn_population.cc: per-weight Gaussian perturbation with self-adaptive sigma, using local Gaussian in test builds
- [x] T063 [US4] Implement nn_tournament_select() in autoc/nn_population.cc: tournament selection for minimization (lower fitness = better)
- [x] T064 [US4] Implement nn_init_population() in autoc/nn_population.cc: create NNPopulation with Xavier-initialized individuals
- [x] T065 [US4] Implement nn_evolve_generation() in autoc/nn_population.cc: elitism + tournament selection → crossover → mutation
- [x] T066 [US4] Add ControllerType config option parsing in autoc/config_manager.cc: "GP" | "NN"
- [x] T067 [US4] Add NNTopology config option parsing in autoc/config_manager.cc: comma-separated layer sizes
- [x] T068 [US4] Add NNMutationSigma config option parsing in autoc/config_manager.cc: initial mutation sigma, default 0.1
- [x] T069 [US4] Add NNCrossoverAlpha config option parsing in autoc/config_manager.cc: BLX-alpha blend factor, default -1
- [x] T070 [US4] Add NNWeightFile config option parsing in autoc/config_manager.cc: weight file path for eval mode, default nn_weights.dat
- [x] T071 [US4] Add NNInitMethod config option parsing in autoc/config_manager.cc: xavier or uniform, default xavier
- [x] T072 [US4] Update autoc/autoc.ini with all NN config parameters (commented out, ControllerType defaults to GP)
- [x] T073 [US4] Update autoc/autoc-eval.ini with NN config parameters for eval mode (NNWeightFile, ControllerType)
- [x] T074 [US4] Wire NNPopulation into autoc main loop in autoc/autoc.cc: runNNEvolution() function with dedicated NN loop, branched via ControllerType config
- [x] T075 [US4] Wire NNControllerBackend into evalTask() dispatch: NN individuals serialized to EvalData.gp with ControllerType::NEURAL_NET, evaluated by minisim NNControllerBackend
- [ ] ~~T076 [US4] Integrate with elite store~~ — deferred (NN has its own elitism in nn_evolve_generation, GP-style elite reeval not needed for MVP)
- [ ] ~~T077 [US4] Add EvalLogger::logNNWeightStats()~~ — deferred (basic gen stats logged inline, detailed weight stats for later)
- [ ] ~~T078 [US4] Update data.dat output~~ — deferred (NN logs to nn-data.dat with basic format)
- [ ] ~~T079 [US4] Update data.stc output~~ — deferred (NN gen stats logged to nn-data.stc inline)
- [x] T080 [US4] Update S3 archive output: best NNGenome per generation saved to S3 with nn- prefix
- [x] T081 [US4] Verify build and all tests: `cd ~/GP && make` — 174 tests pass
- [x] T082 [US4] Integration test: run autoc with ControllerType=NN for 5 generations — done

**Checkpoint**: Full NN evolution pipeline working. Can run `./build/autoc` with ControllerType=NN and see NN weights evolving, fitness improving, results stored in S3 with `nn-` prefix.

---

## Phase 6: Renderer & CRRCSim Integration — US6

**Goal**: Update renderer and CRRCSim to handle NN archives and ControllerType awareness.

**Independent Test**: Renderer can browse both GP (`autoc-*/`) and NN (`nn-*/`) S3 archives. CRRCSim can load NN weight files for live visualization.

**Functional Requirements**: FR-016, FR-017

### Implementation for US6

- [x] T083 [US6] Update renderer S3 browsing and fitness extraction: unified autoc- prefix, NN01 magic detection in extractFitnessFromGP()
- [ ] ~~T084 [US6] Add ControllerType-aware controller metadata display in autoc/renderer.cc~~ — deferred (renderer never dumped GP tree either, not needed for NN)
- [x] T085 [US6] Update CRRCSim inputdev_autoc.cpp: replace magic-byte heuristic with evalData.controllerType enum dispatch, add NN controller path (include nn_serialization.h, nn_evaluator_portable.h, add isNeuralNetData flag, NNGenome storage, NNControllerBackend eval)
- [x] T086 [US6] Verify CRRCSim build: `cd ~/crsim/crrcsim-0.9.13 && rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Debug .. && make`
- [x] T087 [US6] Verify renderer build: `cd ~/GP && make`

**Checkpoint**: Renderer browses NN and GP archives. CRRCSim handles NN format in live visualization.

---

## Phase 7: Embedded Code Generation & Deployment — US5

**Goal**: Generate standalone C++ inference code for xiao-gp embedded deployment. Verify bit-exact parity with desktop inference.

**Independent Test**: Generated nn_program_generated.cpp compiles on both desktop and PlatformIO. Control outputs match desktop nn_forward() for same weights and inputs.

**Functional Requirements**: FR-008, FR-013

### Implementation for US5

- [x] T088 [US5] Implement nn2cpp tool in autoc/nn2cpp.cc: read NN01 weight file, generate nn_program_generated.cpp with embedded weights as static const float[]. Two modes: portable (calls nn_forward) and unrolled (explicit layer loops with fast_tanh). CLI: -i input, -o output, -f funcname, -u unrolled
- [x] T089 [US5] Add nn2cpp build target in autoc/CMakeLists.txt
- [x] T090 [US5] Generate generatedNNProgram() function with same signature as generatedGPProgram(PathProvider&, AircraftState&, gp_scalar). Header: nn_program.h in xiao-gp/include/
- [x] T091 [US5] NN eval mode in autoc.cc: ControllerType=NN + EvaluateMode=1 loads NNWeightFile, evaluates via minisim, reports fitness. Verified exact fitness match (346234.098253) with training run.
- [x] T092 [US5] NN files available via symlink (include/GP -> ~/GP) — no copy needed, shared-source works automatically
- [x] T093 [US5] Generate test nn_program_generated.cpp and place in ~/xiao-gp/generated/
- [x] T094 [US5] Update ~/xiao-gp/src/msplink.cpp: #ifdef USE_NN_PROGRAM selects generatedNNProgram() vs generatedGPProgram()
- [x] T095 [US5] Update ~/xiao-gp/platformio.ini: added xiaoblesense_nn and xiaoble_nn envs with -DUSE_NN_PROGRAM flag and nn_evaluator_portable.cc source
- [x] T096 [US5] Verify xiao-gp PlatformIO build: `cd ~/xiao-gp && pio run -e xiaoblesense_nn` — fixed recordErrorHistory() call (missing distance param from spec 012), builds clean
- [x] T097 [US5] Verify all three repo builds pass: GP (`cd ~/GP && make`), CRRCSim (`cd ~/crsim/crrcsim-0.9.13/build && make`), xiao-gp (`pio run -e xiaoblesense_nn`) — all pass

**Checkpoint**: End-to-end deployment pipeline: evolve NN → nnextractor → nn2cpp → xiao-gp firmware. Bit-exact inference on desktop and embedded.

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Training validation, documentation, final integration

- [ ] ~~T098 [P] Run first full training comparison: NN vs GP on same scenario set~~ — deferred to next feature (NN ES can't match GP with current approach; need segment-based scoring or gradient methods)
- [x] T099 [P] Update CLAUDE.md: added NN build commands (nnextractor, nn2cpp), config options, NN eval architecture, NN workflows
- [x] T100 [P] Update specs/BACKLOG.md: marked 013 as [ACTIVE], unify-eval as [DONE]
- [x] T101 Run release checklist from constitution: all tests pass, all three repo builds pass, CLAUDE.md updated
- [ ] ~~T102 Validate success criteria SC-001 through SC-005~~ — deferred (SC-001 not achieved: NN hasn't matched GP fitness; needs training approach rethink)
- [x] T103 [P] Revisit minisim NN logging — added [AUTOC_GP_STRING] stderr logging for NN (topology, weight count, generation) matching GP/bytecode pattern
- [ ] ~~T104 [P] Cleanup: GPrand() vs local RNG, GP_BUILD/GP_TEST hacks~~ — deferred to GP/NN architecture decision (GP rip-out)
- [ ] ~~T105 [P] Cleanup: controller pluggability~~ — deferred to GP/NN architecture decision (GP rip-out)
- [x] T106 Add determinism checker to NN evolution loop — re-evaluate best individual at end of each generation, compare reeval fitness to stored fitness with bitwiseEqual(). Logs NN_ELITE_SAME or NN_ELITE_DIVERGED with delta.
- [ ] ~~T107 [P] Refactor: unify fitness computation~~ — deferred to GP/NN architecture decision (GP rip-out)
- [ ] ~~T108 [P] Cleanup: consolidate xiao-gp platformio.ini GP/NN env duplication~~ — deferred to GP rip-out (GP envs simply get deleted)
- [ ] ~~T109 [P] NN elite reeval~~ — deferred (determinism checker works, full reeval is polish)
- [x] T110 [P] Compile-time NN topology constants: created nn_topology.h with constexpr constants (NN_INPUT_COUNT, NN_OUTPUT_COUNT, NN_TOPOLOGY[], NN_WEIGHT_COUNT). Removed NNTopology from config parsing and ini files. Updated aircraft_state.h, nn_evaluator_portable.cc/h, autoc.h/cc, config_manager.cc, CLAUDE.md, and all test files.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Unify Eval / US1 (Phase 2)**: Depends on Setup — BLOCKS all subsequent phases
- **NN Core / US2 (Phase 3)**: Depends on US1 (needs ControllerBackend interface)
- **Serialization / US3 (Phase 4)**: Depends on US2 (needs NNGenome struct and nn_forward)
- **Evolution / US4 (Phase 5)**: Depends on US2 + US3 (needs NNGenome, forward pass, serialization)
- **Renderer & CRRCSim / US6 (Phase 6)**: Depends on US3 (needs NN serialization format and S3 naming)
- **Embedded / US5 (Phase 7)**: Depends on US4 (needs trained weights to generate code from)
- **Polish (Phase 8)**: Depends on US4 minimum; US5 + US6 for full validation

### Critical Path

```
Setup → US1 (unify eval) → US2 (NN core) → US3 (serialization) → US4 (evolution) → US5 (embedded) → Polish
                                                                 ↘ US6 (renderer/crrcsim) ↗
```

US6 (renderer/CRRCSim) can run in parallel with US4 after US3 completes.

### Parallel Opportunities Within Phases

- **Phase 1**: T001-T006 all create independent stub files — fully parallel
- **Phase 2**: T010-T012 test tasks are parallel; T019-T020 backend wrappers are parallel
- **Phase 3**: T025-T032 test tasks are parallel; T033+T035 (struct + LUT) are parallel
- **Phase 4**: T041-T044 test tasks are parallel; T045-T047 (serialize/deserialize/detect) are parallel
- **Phase 5**: T055-T059 test tasks are parallel; T060-T065 (population + operators) partially parallel
- **Phase 6**: T083-T085 independent file modifications — parallel
- **Phase 7**: T088-T090 sequential (nn2cpp tool → function → build)

---

## Implementation Strategy

### MVP (Through US2 — NN Core)

1. Complete Phase 1: Setup (stub files, build system)
2. Complete Phase 2: Unify eval pipeline (US1)
3. Complete Phase 3: NN forward pass + controller backend (US2)
4. **STOP and VALIDATE**: Run unified eval with hardcoded NN weights, verify correct sensor gathering and control output

### Incremental Delivery

1. Setup + US1 → Unified eval pipeline (value: cleaner codebase, unblocks 011-gpu-native too)
2. + US2 → NN forward pass works (value: can test NN controller manually)
3. + US3 → Serialization/archive (value: NN individuals can be saved/loaded/transported)
4. + US4 → Full evolution (value: can evolve NN controllers end-to-end)
5. + US6 → Renderer/CRRCSim (value: can visualize NN evolution results)
6. + US5 → Embedded deployment (value: evolved NN flies on real hardware)

---

## Notes

- Constitution requires Testing-First: write tests before implementation, verify they fail
- Constitution requires Build Stability: verify build after every phase
- **Build commands**: When a phase adds new source files or targets to CMakeLists.txt, use full rebuild: `cd ~/GP/autoc && bash rebuild.sh`. For phases that only modify existing files, incremental build suffices: `cd ~/GP && make`. Phases requiring full rebuild: Phase 1 (new source/test files), Phase 4/US3 (nnextractor target), Phase 7/US5 (nn2cpp target). CRRCSim rebuild: `cd ~/crsim/crrcsim-0.9.13 && rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Debug .. && make`
- Constitution grants Dual-Mode Parity exemption: NN is a third mode, not a reimplementation of GP
- NN desktop ↔ embedded parity IS required: same weights must produce bit-exact same outputs
- All sensor functions already exist in gp_evaluator_portable.cc — reuse, don't rewrite
- No backward compatibility with GP/bytecode archives — intentional fork in the road
- GP-only params silently ignored when ControllerType=NN (backward compatible autoc.ini)
- S3 archive prefix: `nn-{timestamp}/` for NN, `autoc-{timestamp}/` for GP
- **Architecture risk**: US1 must eliminate eval pipeline duplication before adding NN. Key deliverables: ControllerType enum in EvalData, shared scenario builder, single evalTask(), ControllerFactory in minisim. See plan.md "Architectural Risk" section.

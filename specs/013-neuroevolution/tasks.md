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

- [ ] T001 Create stub header autoc/eval_backend.h with ControllerBackend interface declaration (virtual evaluate, virtual getName)
- [ ] T002 [P] Create stub files autoc/fitness_computer.h and autoc/fitness_computer.cc
- [ ] T003 [P] Create stub files autoc/eval_logger.h and autoc/eval_logger.cc
- [ ] T004 [P] Create stub files autoc/nn_evaluator_portable.h and autoc/nn_evaluator_portable.cc with normalization constants (NORM_ANGLE, NORM_DIST, NORM_VEL, NORM_RATE)
- [ ] T005 [P] Create stub files autoc/nn_serialization.h and autoc/nn_serialization.cc
- [ ] T006 [P] Create stub files autoc/nn_population.h and autoc/nn_population.cc
- [ ] T007 Add all new source files to autoc/CMakeLists.txt: nn_evaluator_portable, fitness_computer, eval_logger, nn_serialization, nn_population
- [ ] T008 Add test targets to autoc/CMakeLists.txt: nn_evaluator_tests.cc, nn_serialization_tests.cc, nn_population_tests.cc, fitness_computer_tests.cc in autoc/tests/
- [ ] T009 Verify build compiles with stubs: `cd ~/GP && make`

---

## Phase 2: Unify Evaluation Pipeline — US1 (Absorbs 007-unify-eval)

**Goal**: Extract ~350 lines of duplicated eval logic into shared components. Define ControllerBackend interface so GP tree, bytecode, and NN all use the same fitness/logging/scenario paths.

**Independent Test**: Both GP tree mode (EvaluateMode=0) and bytecode mode (EvaluateMode=1) produce identical fitness values after refactor.

**Functional Requirements**: FR-008b, FR-012

### Tests for US1

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T010 [US1] Write parity test in autoc/tests/fitness_computer_tests.cc: verify FitnessComputer::computeStepPenalty matches existing inline fitness calc for known distance/attitude values
- [ ] T011 [P] [US1] Write parity test: verify FitnessComputer::computeCrashPenalty matches existing crash penalty for known fraction_completed values
- [ ] T012 [P] [US1] Write parity test: verify FitnessComputer::computeAttitudeScale matches existing attitude scale computation

### Implementation for US1

- [ ] T013 [US1] Define ControllerBackend interface in autoc/eval_backend.h: virtual evaluate(AircraftState&, PathProvider&) and virtual getName() methods
- [ ] T014 [US1] Add ControllerType enum to EvalData in autoc/autoc.h: GP_TREE=0, BYTECODE=1, NEURAL_NET=2 — explicit type tag replaces magic byte heuristics in minisim
- [ ] T015 [US1] Implement FitnessComputer in autoc/fitness_computer.cc: extract computeStepPenalty, computeCrashPenalty, computeAttitudeScale from autoc.cc lines ~1278-1350
- [ ] T016 [US1] Implement EvalLogger in autoc/eval_logger.cc: extract logStepHeader, logStep, logGenerationStats, logBestController from autoc.cc and autoc-eval.cc
- [ ] T017 [US1] Extract scenario metadata builder: shared function buildScenarioMetadata() from the ~100 lines duplicated between MyGP::evalTask() and BytecodeEvaluationGP::evalTask() (wind/entry variations, demetic grouping, path timing)
- [ ] T018 [US1] Extract GPrandGaussian() utility from variation_generator.h: Box-Muller transform over GPrand(), replaces 4 duplicate inline lambdas. Shared by variation_generator.h and nn_population.cc. Place in a shared header (e.g. autoc/gp_math_utils.h).
- [ ] T019 [US1] Implement GPTreeBackend in autoc/autoc-eval.cc: wrap existing MyGene::evaluate() behind ControllerBackend interface
- [ ] T020 [P] [US1] Implement BytecodeBackend in autoc/autoc-eval.cc: wrap existing GPBytecodeInterpreter::evaluate() behind ControllerBackend interface
- [ ] T021 [US1] Merge evalTask() into single implementation: uses buildScenarioMetadata(), ControllerBackend for eval dispatch, FitnessComputer, EvalLogger. Eliminate BytecodeEvaluationGP inline subclass.
- [ ] T022 [US1] Update minisim.cc format detection: use ControllerType tag from EvalData instead of magic byte heuristics. Add ControllerFactory::create(type, payload) dispatch.
- [ ] T023 [US1] Verify parity: run GP tree mode and bytecode mode with same inputs, confirm identical fitness values
- [ ] T024 [US1] Verify build and tests: `cd ~/GP && make && cd build && ctest --output-on-failure`

**Checkpoint**: Evaluation pipeline unified. Single evalTask() with shared scenario builder, FitnessComputer, EvalLogger, and ControllerBackend interface. BytecodeEvaluationGP inline class eliminated. Minisim uses explicit type tags. All existing tests pass. Ready for NN backend.

---

## Phase 3: NN Evaluator Core — US2

**Goal**: Implement portable NN forward pass with tanh LUT. Create NNControllerBackend conforming to Phase 2 interface.

**Independent Test**: Forward pass produces correct outputs for known weight matrices. Output values always in [-1, 1]. Deterministic for same inputs.

**Functional Requirements**: FR-001, FR-004, FR-005, FR-011, FR-014

### Tests for US2

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T025 [P] [US2] Write forward pass test in autoc/tests/nn_evaluator_tests.cc: single-layer identity weights produce expected output
- [ ] T026 [P] [US2] Write forward pass test: multi-layer 14-16-8-3 topology with known weights produces expected outputs (hand-computed reference)
- [ ] T027 [P] [US2] Write output range test: random weights always produce outputs in [-1, 1] via tanh
- [ ] T028 [P] [US2] Write tanh LUT test: fast_tanh() matches std::tanh() within 1e-3 tolerance across [-5, 5]
- [ ] T029 [P] [US2] Write Xavier init test: initialized weights have zero mean and expected variance (1/fan_in) within statistical bounds
- [ ] T030 [P] [US2] Write determinism test: same weights + same inputs produce bit-exact same outputs across repeated calls
- [ ] T031 [P] [US2] Write weight count test: topology {14,16,8,3} produces exactly 403 weights
- [ ] T032 [P] [US2] Write input normalization test: verify sensor values are divided by correct NORM_* constants before forward pass

### Implementation for US2

- [ ] T033 [US2] Define NNGenome struct in autoc/nn_evaluator_portable.h: weights vector, topology vector, fitness, generation, mutation_sigma per data-model.md
- [ ] T034 [US2] Implement nn_forward() in autoc/nn_evaluator_portable.cc: feedforward pass with row-major layer-sequential weight layout, tanh activation
- [ ] T035 [US2] Implement fast_tanh() LUT in autoc/nn_evaluator_portable.cc: 512-entry table, domain [-5, 5], linear interpolation (same pattern as existing sin LUT in gp_evaluator_portable.cc)
- [ ] T036 [US2] Implement nn_weight_count() utility: compute total weights+biases from topology vector
- [ ] T037 [US2] Implement nn_xavier_init() in autoc/nn_evaluator_portable.cc: Xavier/Glorot initialization per topology, using GPrand() uniform samples (not std::random)
- [ ] T038 [US2] Implement nn_gather_inputs() in autoc/nn_evaluator_portable.cc: call 14 sensor functions (executeGetDPhi, executeGetDist, etc.), apply NORM_ANGLE/NORM_DIST/NORM_VEL/NORM_RATE normalization, build input vector
- [ ] T039 [US2] Implement NNControllerBackend in autoc/nn_evaluator_portable.cc: wraps nn_gather_inputs + nn_forward + setPitch/setRoll/setThrottle, conforms to ControllerBackend interface from eval_backend.h
- [ ] T040 [US2] Verify build and tests: `cd ~/GP && make && cd build && ctest --output-on-failure`

**Checkpoint**: NN forward pass works, tanh LUT accurate, NNControllerBackend plugs into unified eval pipeline. Can manually test NN eval with hardcoded weights.

---

## Phase 4: NN Serialization & Archive — US3

**Goal**: Binary serialization format for NN genomes (RPC transport, S3 storage, minisim detection). Incompatible with GP/bytecode formats by design. NN runs use `nn-{timestamp}/` S3 prefix.

**Independent Test**: Round-trip serialize/deserialize preserves all NNGenome fields exactly. Minisim detects NN format by magic number.

**Functional Requirements**: FR-006, FR-007, FR-008a, FR-017

### Tests for US3

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T041 [P] [US3] Write round-trip test in autoc/tests/nn_serialization_tests.cc: serialize NNGenome → deserialize → all fields match exactly (weights, topology, fitness, generation)
- [ ] T042 [P] [US3] Write magic number test: serialized data starts with "NN01" magic bytes
- [ ] T043 [P] [US3] Write format detection test: NN data detected as NN, GP tree data not detected as NN, bytecode data not detected as NN
- [ ] T044 [P] [US3] Write corrupt data test: truncated/invalid data returns error, does not crash

### Implementation for US3

- [ ] T045 [US3] Implement nn_serialize() in autoc/nn_serialization.cc: write NNGenome to binary format per data-model NNSerializationFormat (magic "NN01", topology, weights, metadata)
- [ ] T046 [US3] Implement nn_deserialize() in autoc/nn_serialization.cc: read binary format back to NNGenome with validation
- [ ] T047 [US3] Implement nn_detect_format() in autoc/nn_serialization.cc: check first 4 bytes for "NN01" magic
- [ ] T048 [US3] Add Boost binary serialization adapter for NNGenome in autoc/nn_serialization.cc (for RPC transport via EvalData.gp blob)
- [ ] T049 [US3] Integrate NN format detection in autoc/minisim.cc: extend ControllerFactory from US1 to handle NEURAL_NET type
- [ ] T050 [US3] Implement NNInterpreter in autoc/minisim.cc: load NN weights from payload, evaluate via nn_forward() + sensor gathering, registered with ControllerFactory
- [ ] T051 [US3] Update S3 archive key generation in autoc/autoc.cc: use `nn-{timestamp}/` prefix when ControllerType=NN (vs existing `autoc-{timestamp}/` for GP)
- [ ] T052 [US3] Implement nnextractor tool in autoc/nnextractor.cc: extract best NNGenome from S3 archive → standalone weight file
- [ ] T053 [US3] Add nnextractor build target in autoc/CMakeLists.txt
- [ ] T054 [US3] Verify build and tests: `cd ~/GP && make && cd build && ctest --output-on-failure`

**Checkpoint**: NN genomes can be serialized/deserialized for RPC and S3. Minisim auto-detects NN format. S3 archives use `nn-` prefix. nnextractor extracts best weights.

---

## Phase 5: Evolution Integration — US4

**Goal**: Wire NNPopulation into existing evolution engine. Arithmetic crossover, Gaussian mutation, tournament selection on weight vectors. Full ControllerType config with all 6 NN-specific params.

**Independent Test**: Run NN evolution for 5 generations on a trivial problem (constant target) and verify fitness improves.

**Functional Requirements**: FR-002, FR-003, FR-009, FR-010, FR-014, FR-015

### Tests for US4

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T055 [P] [US4] Write crossover test in autoc/tests/nn_population_tests.cc: arithmetic crossover of two known genomes produces expected blended weights
- [ ] T056 [P] [US4] Write mutation test: Gaussian mutation changes weights, mean change is near zero over many trials, sigma controls spread
- [ ] T057 [P] [US4] Write self-adaptive sigma test: mutation_sigma itself mutates and stays positive
- [ ] T058 [P] [US4] Write population init test: all individuals have correct weight count, all weights finite, no NaN/Inf
- [ ] T059 [P] [US4] Write tournament selection test: higher fitness individuals selected more often over many tournaments

### Implementation for US4

- [ ] T060 [US4] Implement NNPopulation struct in autoc/nn_population.h: vector of NNGenome, shared topology, generation counter, best_fitness/best_index per data-model.md
- [ ] T061 [US4] Implement nn_arithmetic_crossover() in autoc/nn_population.cc: BLX-alpha blend of two parent weight vectors, alpha from NNCrossoverAlpha config (-1 = uniform [0,1] via GPrand())
- [ ] T062 [US4] Implement nn_gaussian_mutation() in autoc/nn_population.cc: per-weight Gaussian perturbation with self-adaptive sigma, using GPrandGaussian() from T018
- [ ] T063 [US4] Implement nn_tournament_select() in autoc/nn_population.cc: reuse existing tournament selection logic adapted for NNGenome fitness
- [ ] T064 [US4] Implement nn_init_population() in autoc/nn_population.cc: create NNPopulation with Xavier-initialized individuals (or uniform per NNInitMethod config), all via GPrand()
- [ ] T065 [US4] Implement nn_evolve_generation() in autoc/nn_population.cc: selection → crossover → mutation → evaluation loop for one generation
- [ ] T066 [US4] Add ControllerType config option parsing in autoc/config_manager.cc: GP | NN (FR-009)
- [ ] T067 [US4] Add NNTopology config option parsing in autoc/config_manager.cc: comma-separated layer sizes (FR-010)
- [ ] T068 [US4] Add NNMutationSigma config option parsing in autoc/config_manager.cc: initial mutation sigma, default 0.1
- [ ] T069 [US4] Add NNCrossoverAlpha config option parsing in autoc/config_manager.cc: BLX-alpha blend factor, default -1
- [ ] T070 [US4] Add NNWeightFile config option parsing in autoc/config_manager.cc: weight file path for eval mode, default nn_weights.dat
- [ ] T071 [US4] Add NNInitMethod config option parsing in autoc/config_manager.cc: xavier or uniform, default xavier
- [ ] T072 [US4] Update autoc/autoc.ini with all NN config parameters (commented out, ControllerType defaults to GP)
- [ ] T073 [US4] Update autoc/autoc-eval.ini with NN config parameters for eval mode (NNWeightFile, ControllerType)
- [ ] T074 [US4] Wire NNPopulation into autoc main loop in autoc/autoc.cc: when ControllerType=NN, use nn_evolve_generation() instead of GP evolution
- [ ] T075 [US4] Wire NNControllerBackend into evalTask() dispatch: create backend from NNGenome weights when evaluating NN individual
- [ ] T076 [US4] Integrate with elite store in autoc/autoc.cc: re-evaluate elite NN individuals, track best across generations (reuse existing GP elite store)
- [ ] T077 [US4] Add EvalLogger::logNNWeightStats() in autoc/eval_logger.cc: per-layer weight statistics (mean, stdev, min, max) to data.dat each generation
- [ ] T078 [US4] Update data.dat output: when ControllerType=NN, call logNNWeightStats instead of GP S-expression dump
- [ ] T079 [US4] Update data.stc output: log NN-specific metrics (weight magnitude per layer, mutation sigma stats) via EvalLogger
- [ ] T080 [US4] Update S3 archive output in autoc/autoc.cc: store NNPopulation with topology metadata per generation
- [ ] T081 [US4] Verify build and all tests: `cd ~/GP && make && cd build && ctest --output-on-failure`
- [ ] T082 [US4] Integration test: run autoc with ControllerType=NN for 5 generations, verify fitness values are computed and decrease (improve) across generations

**Checkpoint**: Full NN evolution pipeline working. Can run `./build/autoc` with ControllerType=NN and see NN weights evolving, fitness improving, results stored in S3 with `nn-` prefix.

---

## Phase 6: Renderer & CRRCSim Integration — US6

**Goal**: Update renderer and CRRCSim to handle NN archives and ControllerType awareness.

**Independent Test**: Renderer can browse both GP (`autoc-*/`) and NN (`nn-*/`) S3 archives. CRRCSim can load NN weight files for live visualization.

**Functional Requirements**: FR-016, FR-017

### Implementation for US6

- [ ] T083 [US6] Update renderer S3 regex in autoc/renderer.cc: match both `autoc-.*/gen(\d+)\.dmp` (GP) and `nn-.*/gen(\d+)\.dmp` (NN) prefixes
- [ ] T084 [US6] Add ControllerType-aware controller metadata display in autoc/renderer.cc: show NN weight summary instead of GP tree dump when format detected as NN
- [ ] T085 [US6] Update CRRCSim inputdev_autoc.cpp: add NN format detection alongside GP/bytecode for live visualization in ~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp
- [ ] T086 [US6] Verify CRRCSim build: `cd ~/crsim/crrcsim-0.9.13/build && make`
- [ ] T087 [US6] Verify renderer build: `cd ~/GP && make`

**Checkpoint**: Renderer browses NN and GP archives. CRRCSim handles NN format in live visualization.

---

## Phase 7: Embedded Code Generation & Deployment — US5

**Goal**: Generate standalone C++ inference code for xiao-gp embedded deployment. Verify bit-exact parity with desktop inference.

**Independent Test**: Generated nn_program_generated.cpp compiles on both desktop and PlatformIO. Control outputs match desktop nn_forward() for same weights and inputs.

**Functional Requirements**: FR-008, FR-013

### Implementation for US5

- [ ] T088 [US5] Implement nn2cpp tool in autoc/nn2cpp.cc: read weight file, generate nn_program_generated.cpp with embedded weights as static const float[], unrolled layer loops, fast_tanh() calls
- [ ] T089 [US5] Add nn2cpp build target in autoc/CMakeLists.txt
- [ ] T090 [US5] Generate generatedNNProgram() function with same signature as generatedGPProgram(PathProvider&, AircraftState&, gp_scalar)
- [ ] T091 [US5] Write desktop parity test in autoc/tests/nn_evaluator_tests.cc: nn2cpp-generated code produces identical outputs to nn_forward() for same weights/inputs
- [ ] T092 [US5] Copy nn_evaluator_portable.cc/h to ~/xiao-gp/include/GP/autoc/ (shared-source pattern matching gp_evaluator_portable)
- [ ] T093 [US5] Generate test nn_program_generated.cpp and place in ~/xiao-gp/generated/
- [ ] T094 [US5] Update ~/xiao-gp/src/msplink.cpp: add build-time selection between generatedGPProgram() and generatedNNProgram()
- [ ] T095 [US5] Update ~/xiao-gp/platformio.ini: add build flag to select GP vs NN generated program
- [ ] T096 [US5] Verify xiao-gp PlatformIO build: `cd ~/xiao-gp && pio run`
- [ ] T097 [US5] Verify all three repo builds pass: GP (`cd ~/GP && make`), CRRCSim (`cd ~/crsim/crrcsim-0.9.13/build && make`), xiao-gp (`cd ~/xiao-gp && pio run`)

**Checkpoint**: End-to-end deployment pipeline: evolve NN → nnextractor → nn2cpp → xiao-gp firmware. Bit-exact inference on desktop and embedded.

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Training validation, documentation, final integration

- [ ] T098 [P] Run first full training comparison: NN (14-16-8-3, pop=500, 50 gen) vs GP on same scenario set, compare convergence curves and control smoothness
- [ ] T099 [P] Update CLAUDE.md: add NN-related build commands, config options, new executables (nnextractor, nn2cpp), ControllerType param
- [ ] T100 [P] Update specs/BACKLOG.md: mark 013 as in-progress/complete, update related deferred items (4D fitness, smoothness)
- [ ] T101 Run release checklist from constitution: all tests pass, all three repo builds pass, CLAUDE.md updated
- [ ] T102 Validate success criteria SC-001 through SC-005: lower fitness than GP, consistent convergence, <1ms embedded inference, smooth control outputs, cross-platform parity

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
- Constitution requires Build Stability: verify `cd ~/GP && make` after every phase
- Constitution grants Dual-Mode Parity exemption: NN is a third mode, not a reimplementation of GP
- NN desktop ↔ embedded parity IS required: same weights must produce bit-exact same outputs
- All sensor functions already exist in gp_evaluator_portable.cc — reuse, don't rewrite
- No backward compatibility with GP/bytecode archives — intentional fork in the road
- GP-only params silently ignored when ControllerType=NN (backward compatible autoc.ini)
- S3 archive prefix: `nn-{timestamp}/` for NN, `autoc-{timestamp}/` for GP
- **Architecture risk**: US1 must eliminate eval pipeline duplication before adding NN. Key deliverables: ControllerType enum in EvalData, shared scenario builder, single evalTask(), ControllerFactory in minisim. See plan.md "Architectural Risk" section.

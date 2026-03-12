# Tasks: 013-neuroevolution

**Input**: Design documents from `/specs/013-neuroevolution/`
**Prerequisites**: plan.md, spec.md, research.md, data-model.md

**Tests**: Required by project constitution (Testing-First principle). Tests written before implementation.

**Organization**: Tasks grouped by spec implementation phase. Each phase is an independently testable increment.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which implementation phase (US1=Unify Eval, US2=NN Core, US3=Serialization, US4=Evolution, US5=Embedded Codegen)

---

## Phase 1: Setup (Build System & File Scaffolding)

**Purpose**: Create new source files and update CMakeLists.txt for NN support

- [ ] T001 Add NN source files to CMakeLists.txt: nn_evaluator_portable.cc/h, fitness_computer.cc/h, eval_logger.cc/h, eval_backend.h, nn_serialization.cc/h in autoc/CMakeLists.txt
- [ ] T002 [P] Create stub header autoc/eval_backend.h with ControllerBackend interface declaration
- [ ] T003 [P] Create stub files autoc/fitness_computer.h and autoc/fitness_computer.cc
- [ ] T004 [P] Create stub files autoc/eval_logger.h and autoc/eval_logger.cc
- [ ] T005 [P] Create stub files autoc/nn_evaluator_portable.h and autoc/nn_evaluator_portable.cc
- [ ] T006 [P] Create stub files autoc/nn_serialization.h and autoc/nn_serialization.cc
- [ ] T007 Add test targets to CMakeLists.txt: nn_evaluator_tests.cc, nn_serialization_tests.cc, nn_evolution_tests.cc in autoc/tests/
- [ ] T008 Verify build compiles with stubs: `cd ~/GP && make`

---

## Phase 2: Unify Evaluation Pipeline — US1 (Absorbs 007-unify-eval)

**Goal**: Extract ~350 lines of duplicated eval logic into shared components. Define ControllerBackend interface so GP tree, bytecode, and NN all use the same fitness/logging/scenario paths.

**Independent Test**: Both GP tree mode (EvaluateMode=0) and bytecode mode (EvaluateMode=1) produce identical fitness values after refactor.

**Functional Requirements**: FR-008b, FR-012

### Tests for US1

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T009 [US1] Write parity test in autoc/tests/gp_evaluator_tests.cc: verify FitnessComputer::computeStepPenalty matches existing inline fitness calc for known distance/attitude values
- [ ] T010 [P] [US1] Write parity test: verify FitnessComputer::computeCrashPenalty matches existing crash penalty for known fraction_completed values
- [ ] T011 [P] [US1] Write parity test: verify FitnessComputer::computeAttitudeScale matches existing attitude scale computation

### Implementation for US1

- [ ] T012 [US1] Define ControllerBackend interface in autoc/eval_backend.h: virtual evaluate(AircraftState&, PathProvider&) and virtual getName() methods
- [ ] T013 [US1] Implement FitnessComputer in autoc/fitness_computer.cc: extract computeStepPenalty, computeCrashPenalty, computeAttitudeScale from autoc.cc lines ~1278-1350
- [ ] T014 [US1] Implement EvalLogger in autoc/eval_logger.cc: extract logStepHeader, logStep, logGenerationStats, logBestController from autoc.cc and autoc-eval.cc
- [ ] T015 [US1] Implement GPTreeBackend in autoc/autoc-eval.cc: wrap existing MyGene::evaluate() behind ControllerBackend interface
- [ ] T016 [US1] Implement BytecodeBackend in autoc/autoc-eval.cc: wrap existing GPBytecodeInterpreter::evaluate() behind ControllerBackend interface
- [ ] T017 [US1] Refactor evalTask() in autoc/autoc-eval.cc to use FitnessComputer, EvalLogger, and ControllerBackend instead of inline logic
- [ ] T018 [US1] Refactor BytecodeEvaluationGP::evalTask() in autoc/autoc-eval.cc to use same shared components
- [ ] T019 [US1] Verify parity: run GP tree mode and bytecode mode with same inputs, confirm identical fitness values
- [ ] T020 [US1] Verify build: `cd ~/GP && make` and `cd ~/GP/build && ctest --output-on-failure`

**Checkpoint**: Evaluation pipeline unified. GP tree and bytecode modes use shared FitnessComputer, EvalLogger, and ControllerBackend interface. All existing tests pass. Ready for NN backend.

---

## Phase 3: NN Evaluator Core — US2

**Goal**: Implement portable NN forward pass with tanh LUT. Create NNControllerBackend conforming to Phase 2 interface.

**Independent Test**: Forward pass produces correct outputs for known weight matrices. Output values always in [-1, 1]. Deterministic for same inputs.

**Functional Requirements**: FR-001, FR-004, FR-005, FR-011, FR-014

### Tests for US2

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T021 [P] [US2] Write forward pass test in autoc/tests/nn_evaluator_tests.cc: single-layer identity weights produce expected output
- [ ] T022 [P] [US2] Write forward pass test: multi-layer 14-16-8-3 topology with known weights produces expected outputs (hand-computed reference)
- [ ] T023 [P] [US2] Write output range test: random weights always produce outputs in [-1, 1] via tanh
- [ ] T024 [P] [US2] Write tanh LUT test: fast_tanh() matches std::tanh() within 1e-3 tolerance across [-5, 5]
- [ ] T025 [P] [US2] Write Xavier init test: initialized weights have zero mean and expected variance (1/fan_in) within statistical bounds
- [ ] T026 [P] [US2] Write determinism test: same weights + same inputs produce bit-exact same outputs across repeated calls
- [ ] T027 [P] [US2] Write weight count test: topology {14,16,8,3} produces exactly 403 weights

### Implementation for US2

- [ ] T028 [US2] Define NNGenome struct in autoc/nn_evaluator_portable.h: weights vector, topology vector, fitness, generation, mutation_sigma per data-model.md
- [ ] T029 [US2] Implement nn_forward() in autoc/nn_evaluator_portable.cc: feedforward pass with layer-major weight layout, tanh activation
- [ ] T030 [US2] Implement fast_tanh() LUT in autoc/nn_evaluator_portable.cc: 512-entry table with linear interpolation (same pattern as existing sin LUT in gp_evaluator_portable.cc)
- [ ] T031 [US2] Implement nn_weight_count() utility: compute total weights+biases from topology vector
- [ ] T032 [US2] Implement nn_xavier_init() in autoc/nn_evaluator_portable.cc: Xavier/Glorot initialization per topology
- [ ] T033 [US2] Implement nn_gather_inputs() in autoc/nn_evaluator_portable.cc: call 14 sensor functions (executeGetDPhi, executeGetDist, etc.) to build input vector
- [ ] T034 [US2] Implement NNControllerBackend in autoc/nn_evaluator_portable.cc: wraps nn_gather_inputs + nn_forward + setPitch/setRoll/setThrottle, conforms to ControllerBackend interface from eval_backend.h
- [ ] T035 [US2] Verify build and tests: `cd ~/GP && make && cd build && ctest --output-on-failure`

**Checkpoint**: NN forward pass works, tanh LUT accurate, NNControllerBackend plugs into unified eval pipeline. Can manually test NN eval with hardcoded weights.

---

## Phase 4: NN Serialization & Archive — US3

**Goal**: Binary serialization format for NN genomes (RPC transport, S3 storage, minisim detection). Incompatible with GP/bytecode formats by design.

**Independent Test**: Round-trip serialize/deserialize preserves all NNGenome fields exactly. Minisim detects NN format by magic number.

**Functional Requirements**: FR-006, FR-007, FR-008a

### Tests for US3

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T036 [P] [US3] Write round-trip test in autoc/tests/nn_serialization_tests.cc: serialize NNGenome → deserialize → all fields match exactly (weights, topology, fitness, generation)
- [ ] T037 [P] [US3] Write magic number test: serialized data starts with "NN01" magic bytes
- [ ] T038 [P] [US3] Write format detection test: NN data detected as NN, GP tree data not detected as NN, bytecode data not detected as NN
- [ ] T039 [P] [US3] Write corrupt data test: truncated/invalid data returns error, does not crash

### Implementation for US3

- [ ] T040 [US3] Implement nn_serialize() in autoc/nn_serialization.cc: write NNGenome to binary format per data-model NNSerializationFormat (magic "NN01", topology, weights, metadata)
- [ ] T041 [US3] Implement nn_deserialize() in autoc/nn_serialization.cc: read binary format back to NNGenome with validation
- [ ] T042 [US3] Implement nn_detect_format() in autoc/nn_serialization.cc: check first 4 bytes for "NN01" magic
- [ ] T043 [US3] Add Boost binary serialization adapter for NNGenome in autoc/nn_serialization.cc (for RPC transport via EvalData.gp blob)
- [ ] T044 [US3] Integrate NN format detection in autoc/minisim.cc: extend existing format detection (GP tree / bytecode / NN) using magic number check
- [ ] T045 [US3] Implement NNInterpreter in autoc/minisim.cc: load NN weights from file, evaluate via nn_forward() + sensor gathering, same interface as GPBytecodeInterpreter
- [ ] T046 [US3] Implement nnextractor tool in autoc/nnextractor.cc: extract best NNGenome from S3 archive → standalone weight file
- [ ] T047 [US3] Add nnextractor build target in autoc/CMakeLists.txt
- [ ] T048 [US3] Verify build and tests: `cd ~/GP && make && cd build && ctest --output-on-failure`

**Checkpoint**: NN genomes can be serialized/deserialized for RPC and S3. Minisim auto-detects NN format. nnextractor extracts best weights.

---

## Phase 5: Evolution Integration — US4

**Goal**: Wire NNPopulation into existing evolution engine. Arithmetic crossover, Gaussian mutation, tournament selection on weight vectors. Configurable via autoc.ini.

**Independent Test**: Run NN evolution for 5 generations on a trivial problem (constant target) and verify fitness improves.

**Functional Requirements**: FR-002, FR-003, FR-009, FR-010, FR-014

### Tests for US4

> **Write tests FIRST, ensure they FAIL before implementation**

- [ ] T049 [P] [US4] Write crossover test in autoc/tests/nn_evolution_tests.cc: arithmetic crossover of two known genomes produces expected blended weights
- [ ] T050 [P] [US4] Write mutation test: Gaussian mutation changes weights, mean change is near zero over many trials, sigma controls spread
- [ ] T051 [P] [US4] Write self-adaptive sigma test: mutation_sigma itself mutates and stays positive
- [ ] T052 [P] [US4] Write population init test: all individuals have correct weight count, all weights finite, no NaN/Inf
- [ ] T053 [P] [US4] Write tournament selection test: higher fitness individuals selected more often over many tournaments

### Implementation for US4

- [ ] T054 [US4] Implement NNPopulation struct in autoc/autoc.h: vector of NNGenome, shared topology, generation counter, best_fitness/best_index per data-model.md
- [ ] T055 [US4] Implement nn_arithmetic_crossover() in autoc/autoc.cc: BLX-alpha blend of two parent weight vectors
- [ ] T056 [US4] Implement nn_gaussian_mutation() in autoc/autoc.cc: per-weight Gaussian perturbation with self-adaptive sigma
- [ ] T057 [US4] Implement nn_tournament_select() in autoc/autoc.cc: reuse existing tournament selection logic adapted for NNGenome fitness
- [ ] T058 [US4] Implement nn_init_population() in autoc/autoc.cc: create NNPopulation with Xavier-initialized individuals
- [ ] T059 [US4] Implement nn_evolve_generation() in autoc/autoc.cc: selection → crossover → mutation → evaluation loop for one generation
- [ ] T060 [US4] Add ControllerType config option in autoc/autoc.ini and parsing in autoc/config_manager.cc: GP | NN
- [ ] T061 [US4] Add NNTopology config option in autoc/autoc.ini: comma-separated layer sizes (e.g., 14,16,8,3)
- [ ] T062 [US4] Add NNMutationSigma config option in autoc/autoc.ini: initial mutation sigma (default 0.1)
- [ ] T063 [US4] Wire NNPopulation into autoc main loop in autoc/autoc.cc: when ControllerType=NN, use nn_evolve_generation() instead of GP evolution
- [ ] T064 [US4] Wire NNControllerBackend into evalTask() dispatch: create backend from NNGenome weights when evaluating NN individual
- [ ] T065 [US4] Integrate with elite store in autoc/autoc.cc: re-evaluate elite NN individuals, track best across generations
- [ ] T066 [US4] Update data.stc output in autoc/autoc.cc: log NN-specific metrics (weight magnitude per layer, mutation sigma stats) via EvalLogger
- [ ] T067 [US4] Update data.dat output: replace GP S-expression dump with NN weight summary (mean/stdev per layer)
- [ ] T068 [US4] Update S3 archive output: store NNPopulation with topology metadata per generation
- [ ] T069 [US4] Verify build and all tests: `cd ~/GP && make && cd build && ctest --output-on-failure`
- [ ] T070 [US4] Integration test: run autoc with ControllerType=NN for 5 generations, verify fitness values are computed and decrease (improve) across generations

**Checkpoint**: Full NN evolution pipeline working. Can run `./build/autoc` with ControllerType=NN and see NN weights evolving, fitness improving, results stored in S3.

---

## Phase 6: Embedded Code Generation & Deployment — US5

**Goal**: Generate standalone C++ inference code for xiao-gp embedded deployment. Verify bit-exact parity with desktop inference.

**Independent Test**: Generated nn_program_generated.cpp compiles on both desktop and PlatformIO. Control outputs match desktop nn_forward() for same weights and inputs.

**Functional Requirements**: FR-008, FR-013

### Implementation for US5

- [ ] T071 [US5] Implement nn2cpp tool in autoc/nn2cpp.cc: read weight file, generate nn_program_generated.cpp with embedded weights, unrolled layer loops, fast_tanh() calls
- [ ] T072 [US5] Add nn2cpp build target in autoc/CMakeLists.txt
- [ ] T073 [US5] Generate generatedNNProgram() function with same signature as generatedGPProgram(PathProvider&, AircraftState&, gp_scalar)
- [ ] T074 [US5] Copy nn_evaluator_portable.cc/h to ~/xiao-gp/include/GP/autoc/ (shared-source pattern matching gp_evaluator_portable)
- [ ] T075 [US5] Generate test nn_program_generated.cpp and place in ~/xiao-gp/generated/
- [ ] T076 [US5] Update ~/xiao-gp/src/msplink.cpp: add build-time selection between generatedGPProgram() and generatedNNProgram()
- [ ] T077 [US5] Verify xiao-gp PlatformIO build: `cd ~/xiao-gp && pio run`
- [ ] T078 [US5] Verify desktop bit-exact parity: nn2cpp-generated code produces identical outputs to nn_forward() for same weights/inputs (test in autoc/tests/nn_evaluator_tests.cc)
- [ ] T079 [US5] Verify all three repo builds pass: GP (`cd ~/GP && make`), CRRCSim (`cd ~/crsim/crrcsim-0.9.13/build && make`), xiao-gp (`cd ~/xiao-gp && pio run`)

**Checkpoint**: End-to-end deployment pipeline: evolve NN → nnextractor → nn2cpp → xiao-gp firmware. Bit-exact inference on desktop and embedded.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Architecture search, documentation, final integration validation

- [ ] T080 [P] Run first full training comparison: NN (14-16-8-3, pop=500, 50 gen) vs GP on same scenario set, compare convergence
- [ ] T081 [P] Update CLAUDE.md: add NN-related build commands, config options, new executables (nnextractor, nn2cpp)
- [ ] T082 [P] Update specs/BACKLOG.md: mark 013 as in-progress/complete, update related items
- [ ] T083 Run release checklist from constitution: all tests pass, all three repo builds pass, CLAUDE.md updated

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Unify Eval / US1 (Phase 2)**: Depends on Setup — BLOCKS all subsequent phases
- **NN Core / US2 (Phase 3)**: Depends on US1 (needs ControllerBackend interface)
- **Serialization / US3 (Phase 4)**: Depends on US2 (needs NNGenome struct and nn_forward)
- **Evolution / US4 (Phase 5)**: Depends on US2 + US3 (needs NNGenome, forward pass, serialization)
- **Embedded / US5 (Phase 6)**: Depends on US4 (needs trained weights to generate code from)
- **Polish (Phase 7)**: Depends on US4 minimum; US5 for full validation

### Critical Path

```
Setup → US1 (unify eval) → US2 (NN core) → US3 (serialization) → US4 (evolution) → US5 (embedded) → Polish
```

Phases are strictly sequential because each builds on the previous phase's artifacts. Within each phase, test tasks marked [P] can run in parallel, and some implementation tasks marked [P] can run in parallel.

### Parallel Opportunities Within Phases

- **Phase 1**: T002-T006 all create independent stub files — fully parallel
- **Phase 2**: T009-T011 test tasks are parallel; T015-T016 backend wrappers are parallel
- **Phase 3**: T021-T027 test tasks are parallel; T028+T030 (struct + LUT) are parallel
- **Phase 4**: T036-T039 test tasks are parallel; T040-T042 (serialize/deserialize/detect) are parallel
- **Phase 5**: T049-T053 test tasks are parallel; T054-T058 (population + operators) partially parallel
- **Phase 6**: T071-T073 sequential (nn2cpp → function → build)

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
5. + US5 → Embedded deployment (value: evolved NN flies on real hardware)

---

## Notes

- Constitution requires Testing-First: write tests before implementation, verify they fail
- Constitution requires Build Stability: verify `cd ~/GP && make` after every phase
- Constitution grants Dual-Mode Parity exemption: NN is a third mode, not a reimplementation of GP
- NN desktop ↔ embedded parity IS required: same weights must produce bit-exact same outputs
- All sensor functions already exist in gp_evaluator_portable.cc — reuse, don't rewrite
- No backward compatibility with GP/bytecode archives — intentional fork in the road

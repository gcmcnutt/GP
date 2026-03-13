# Implementation Plan: 013-neuroevolution

**Branch**: `013-neuroevolution` | **Date**: 2026-03-12 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/013-neuroevolution/spec.md`

## Summary

Replace tree-based GP controller with an evolved fixed-topology neural network (14→16→8→3,
403 parameters, tanh activation). The NN addresses the GP representation ceiling: bang-bang
control dominance, combinatorial node competition, and inability to compute cross-terms from
existing sensors. Evolution uses simple GA with Gaussian mutation on flat weight vectors,
reusing existing population management and fitness infrastructure. Deployment targets both
desktop (training) and nRF52840 Cortex-M4F (embedded flight).

Phase 1 unifies the evaluation pipeline (absorbs 007-unify-eval) to create a clean
`ControllerBackend` interface before adding NN as a third controller type.

## Technical Context

**Language/Version**: C++17 (g++, CMake 3.10+)
**Primary Dependencies**: Eigen3 (vectors, quaternions), Boost (serialization, logging, threads), GoogleTest 1.14.0
**Storage**: S3 for evolution archives, local files for diagnostics (data.dat, data.stc)
**Testing**: GoogleTest 1.14.0 (unit + integration), manual evolution runs (convergence validation)
**Target Platform**: aarch64 Linux desktop (training), ARM Cortex-M4F @ 64 MHz (embedded deployment)
**Project Type**: Library + application (core GP lib, autoc evolution system, xiao-gp embedded)
**Performance Goals**: NN inference <1 ms on nRF52840 (baseline topology); fitness eval identical to GP (~160s/gen at pop=500)
**Constraints**: 256 KB RAM, 1 MB Flash on embedded; no ML framework dependencies; FP32 only

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

### I. Testing-First ✅

All phases include unit tests. Phase 2 (NN core) has known-weight regression tests and
output range validation. Phase 4 (evolution) has convergence smoke tests. Phase 5 (embedded)
has cross-platform parity tests.

### II. Build Stability ✅

All three repos (GP, CRRCSim, xiao-gp) must compile at each phase boundary. Incremental
changes with CI verification. New files added to CMakeLists.txt as they're created.

### III. Dual-Mode Parity — JUSTIFIED EXEMPTION

The NN controller is intentionally **not** parity-constrained with GP tree or bytecode
evaluation. It is a new controller type with its own serialization format and evaluation
path. However, the NN evaluation itself must produce identical results across:
- Desktop training (nn_evaluator_portable.cc)
- Desktop deployment (NN weight file loaded by minisim)
- Embedded deployment (nn2cpp generated code)

This is **NN self-parity** across platforms, not GP/NN cross-parity.

## Project Structure

### Documentation (this feature)

```text
specs/013-neuroevolution/
├── plan.md              # This file
├── research.md          # Phase 0 output (12 decisions)
├── data-model.md        # Phase 1 output (6 entities)
├── spec.md              # Feature specification
└── tasks.md             # Phase 2 output (/speckit.tasks)
```

### Source Code (across 3 repositories)

```text
# GP repository (~/GP)
autoc/
├── nn_evaluator_portable.cc   # [NEW] Portable NN forward pass (desktop + embedded)
├── nn_evaluator_portable.h    # [NEW] NN forward pass header + normalization constants
├── nn_serialization.cc        # [NEW] NNGenome Boost serialization + "NN01" binary format
├── nn_serialization.h         # [NEW] Serialization header
├── nn2cpp.cc                  # [NEW] NN weight → standalone C++ code generator
├── nnextractor.cc             # [NEW] Extract best NN from S3 → weight file
├── gp_math_utils.h            # [NEW] Shared GPrandGaussian() utility (Box-Muller over GPrand)
├── eval_backend.h             # [NEW] ControllerBackend interface
├── fitness_computer.cc        # [NEW] Extracted shared fitness computation
├── fitness_computer.h         # [NEW] Fitness computation header
├── eval_logger.cc             # [NEW] Extracted shared logging (data.dat, data.stc)
├── eval_logger.h              # [NEW] Logging header
├── nn_population.cc           # [NEW] NNPopulation: genomes, selection, crossover, mutation
├── nn_population.h            # [NEW] Population header
├── autoc.cc                   # [MOD] Add ControllerType config, NN population management
├── autoc.h                    # [MOD] Add NN config options, NNGenome struct
├── autoc-eval.cc              # [MOD] Refactor to use ControllerBackend interface
├── autoc.ini                  # [MOD] Add ControllerType, NNTopology, NNMutationSigma
├── autoc-eval.ini             # [MOD] Add NNWeightFile, ControllerType for eval mode
├── config_manager.cc          # [MOD] Parse ControllerType + 6 NN config params
├── minisim.cc                 # [MOD] Add NN format detection alongside GP/bytecode
├── renderer.cc                # [MOD] ControllerType awareness, NN S3 prefix regex
├── CMakeLists.txt             # [MOD] Add new source files
├── tests/
│   ├── nn_evaluator_tests.cc  # [NEW] Forward pass unit tests
│   ├── nn_population_tests.cc # [NEW] Evolution operator tests
│   ├── nn_serialization_tests.cc # [NEW] Serialization round-trip tests
│   └── fitness_computer_tests.cc # [NEW] Extracted fitness computation tests

# CRRCSim repository (~/crsim/crrcsim-0.9.13)
src/mod_inputdev/inputdev_autoc/
├── inputdev_autoc.cpp         # [MOD] Add NN format detection for live visualization

# xiao-gp repository (~/xiao-gp)
├── generated/
│   └── nn_program_generated.cpp  # [NEW] nn2cpp output (drop-in for gp_program_generated.cpp)
├── include/
│   └── nn_evaluator_portable.h   # [SYMLINK/COPY] Shared header from GP
├── src/
│   └── nn_evaluator_portable.cc  # [SYMLINK/COPY] Shared portable evaluator from GP
└── platformio.ini                # [MOD] Build flag to select GP vs NN generated program
```

**Structure Decision**: Follows existing patterns — portable evaluator shared across repos,
code generator produces standalone .cpp files, build system selects controller type.

## Architectural Risk: Eval Pipeline Spaghetti

The current eval pipeline has structural debt that US1 must address before adding NN.

### Current Problems

1. **God class**: `autoc.cc` is 2,630 lines — scenario building, RPC, fitness, S3, config
   dispatch all in one file. `BytecodeEvaluationGP` is an inline nested class.

2. **~100 lines duplicated**: `MyGP::evalTask()` and `BytecodeEvaluationGP::evalTask()`
   share identical scenario metadata building (wind/entry variations, demetic grouping,
   path timing). Adding NN as a third copy would triple this.

3. **Implicit format detection**: `minisim.cc` guesses GP vs bytecode by inspecting the
   first bytes of the payload (`firstBytes[0] >= 0x16`). No type tag in `EvalData`.
   Adding a third format to this heuristic is fragile.

4. **No common evaluator interface**: GP tree uses `MyGene::evaluate(vector<Path>&, MyGP&, gp_scalar)`,
   bytecode uses `GPBytecodeInterpreter::evaluate(AircraftState&, vector<Path>&, gp_scalar)`.
   Different signatures, dispatched via if/else in minisim.

5. **No controller type declaration**: `EvaluateMode` is a bare int (0/1). Config lacks
   an explicit controller type enum or factory pattern.

### US1 Refactoring Strategy

US1 must clean this up **before** NN is added. The approach:

1. **Add `ControllerType` enum** to `EvalData` struct — explicit type tag replaces
   magic byte heuristics. Values: `GP_TREE`, `BYTECODE`, `NEURAL_NET`.

2. **Extract scenario builder** — the 100-line metadata construction block becomes a
   shared function: `buildScenarioMetadata(EvalData&, scenario, variations, ...)`.
   Called once from a single `evalTask()` implementation.

3. **Define `ControllerBackend` interface** — unified `evaluate(AircraftState&, PathProvider&)`
   signature. GPTreeBackend, BytecodeBackend, NNControllerBackend all implement it.

4. **Single evalTask()** — parameterized by controller type. Builds scenario metadata
   (shared), serializes controller payload (type-specific), sends RPC, receives results.
   No subclass duplication.

5. **Factory in minisim** — `ControllerFactory::create(controllerType, payload)` returns
   the right backend. No if/else chains, no byte inspection.

### What NOT to Over-Engineer

- No template metaprogramming — simple virtual dispatch is fine for 3 types
- No plugin system — controller types are compile-time known
- No separate scenario_builder.h — a free function in the existing eval code is sufficient
- The `ControllerBackend` interface is 2 virtual methods, not an abstract framework

### Risk Mitigation

| Risk | Mitigation |
|------|-----------|
| US1 scope creep (refactoring too much) | Strict scope: extract duplication, add type tag, define interface. No other changes to autoc.cc |
| Breaking GP/bytecode parity during refactor | Parity test (T023) runs both modes with identical inputs, must produce identical fitness |
| New interface too rigid for future controllers | Interface is minimal (2 methods). Future controller types add implementations, not interface changes |

## Complexity Tracking

| Violation | Why Needed | Simpler Alternative Rejected Because |
|-----------|------------|-------------------------------------|
| Dual-mode parity exemption | NN is a new controller type, not a GP variant | Forcing NN/GP parity is meaningless — different representations |
| New serialization format ("NN01") | Flat weight arrays are fundamentally different from GP trees | Reusing GP/bytecode format would require wasteful encoding |

# Research: 014-nn-training-signal

**Date**: 2026-03-14

## R1: libgp.a Dependency Surface

**Decision**: Drop libgp.a entirely. No symbols retained.

**Findings**: autoc uses 46 of 103 exported symbols from libgp.a, across 10 files. Key dependencies:

| Category | Symbols | Replacement |
|----------|---------|-------------|
| PRNG | `GPrand()`, `GPsrand()`, `GPRandomPercent()` | `std::mt19937` (already used in NN code) |
| Config | `GPConfiguration`, `GPVariables`, `GPConfigVarInformation` | inih library (see R2) |
| GP Tree | `GP`, `GPGene`, `GPNode`, `GPNodeSet`, `GPAdfNodeSet`, `GPContainer` | Remove entirely (NN-only fork) |
| Population | `GPPopulation`, selection/crossover/mutation | NNPopulation already self-contained |
| Init/Misc | `GPInit()`, `GPRegisterClass()`, `GPExitSystem()` | Remove (GP serialization infrastructure) |
| Load/Save | `GPGene::load/save`, `GP::load/save` | Remove (NN uses "NN01" format) |

**data.dat / data.stc**: Both written by autoc.cc, NOT by libgp. Safe to keep as-is during refactoring.

**Transitive include analysis**: `gp.h` is self-contained (only pulls `<iostream>`). No deep transitive chains. However, `gp.h` defines the *entire* class hierarchy that autoc inherits from (GPObject → GPContainer → GP → GPPopulation, GPGene, GPNode, GPNodeSet, GPAdfNodeSet, GPVariables). The dependency is inheritance-deep, not include-deep. `gpconfig.h` is similarly flat (defines GPConfiguration, GPConfigVarInformation).

**Files requiring GP header removal** (with what each uses from gp.h):

| File | `#include` | What it uses |
|------|-----------|--------------|
| autoc.h | gp.h | MyGP : GP, MyGene : GPGene, MyPopulation : GPPopulation |
| autoc.cc | gp.h, gpconfig.h | GPInit, GPRegisterClass, GPrand, GPsrand, GPVariables, GPConfiguration |
| autoc-eval.cc | gp.h | MyGene::evaluate(), GPNode::value(), tree traversal |
| config_manager.h | gp.h, gpconfig.h | GPVariables (struct), GPConfigVarInformation, GPConfiguration |
| config_manager.cc | (via .h) | GPConfiguration::read(), DATAINT/DATADOUBLE/DATASTRING |
| minisim.cc | gp.h | GP::load(), GPGene, GPRegisterClass, tree deserialization |
| pathgen.cc | gp.h | GPrand() only |
| gpextractor.cc | gp.h | GP::load(), GPGene::NthChild(), tree walking |
| nnextractor.cc | gp.h | GPConfiguration (for ini parsing only) |
| renderer.cc | gp.h | GP::load(), tree deserialization for S3 replay |
| gp_math_utils.h | gp.h | GPrand() only |

**CMakeLists.txt changes**: Remove `include_directories(../include)`, `link_directories(../lib)`, and `gp` from `target_link_libraries`.

**Alternatives considered**: Keeping a minimal GP subset inline — rejected because the entire GP tree evaluation, bytecode, and serialization infrastructure is unused in the NN-only fork.

---

## R2: Config Parser Replacement

**Decision**: Use **inih** (benhoyt/inih) — BSD-licensed, 3 files (<400 LOC), zero config file changes needed.

**Rationale**:
- `autoc.ini` is flat key=value with `#` comments — inih handles this natively with empty section name `""`
- ~45 config keys, all flat (no nesting, no arrays, no sections)
- Types needed: int, double, string — inih provides `GetInteger()`, `GetReal()`, `GetString()`
- Returns `std::string` (solves current `char*`/`strdup` problem)
- Used in Google Fuchsia, systemd — battle-tested

**Current config usage** (config_manager.cc): `GPConfiguration` parses autoc.ini into a table of `{key, type, pointer}` entries. Types are `DATAINT`, `DATADOUBLE`, `DATASTRING`. Simple mechanical translation to inih API.

**Alternatives considered**:
- toml++ (28k LOC, format change required) — overpowered
- SimpleIni (2800 LOC, MIT) — good but heavier than needed
- mINI (900 LOC, requires sections) — would need `[General]` section added
- inipp (300 LOC, MIT) — viable runner-up but less battle-tested

---

## R3: sep-CMA-ES Implementation

**Decision**: Custom implementation (~200-300 LOC C++ with Eigen). No external library dependency.

**Algorithm** (Ros & Hansen 2008): Replaces full N×N covariance matrix with diagonal vector. O(N) per generation instead of O(N³).

**Hyperparameters for N=531**:

| Parameter | Value |
|-----------|-------|
| lambda (population) | 50 |
| mu (parents) | 25 |
| mu_eff | ~13.2 |
| c_sigma | 0.028 |
| d_sigma | 1.03 |
| c_c | 0.0075 |
| c_1 (sep-scaled) | 0.00125 |
| c_mu (sep-scaled) | 0.014 |
| sigma_0 (initial) | 0.3-0.5 |
| E[||N(0,I)||] | ~22.99 |
| min lambda | 22 (formula: 4 + floor(3*ln(531))) |

**Key insight**: sep-CMA-ES learning rates must be scaled by (N+2)/3 = 177.7 to compensate for diagonal approximation.

**Checkpoint state**: ~21 KB total — mean (531), diagonal covariance (531), p_sigma (531), p_c (531), sigma (1), generation (1), best_x (531), best_fitness (1). Population is regenerated from state, not saved.

**Parallelism**: Natural ask-tell pattern. Generate 50 candidates (instant), evaluate all in parallel across 49 scenarios (2,450 sims per generation), update (instant). Fits existing EvalThreads infrastructure.

**Alternatives considered**:
- libcmaes (C++, Eigen-based, LGPL) — full library but LGPL license and heavier than needed
- CMA-ESpp — no sep-CMA-ES variant
- Hansen's c-cmaes — ANSI C reference, no sep variant

---

## R4: Boost Dependency Surface

**Decision**: Remove Boost in phases. Phase 1 (trivial) + Phase 2 (medium) first, Phase 3 (serialization) during transport simplification.

**Component inventory (9 components across 18 files)**:

| Component | Replacement | Effort |
|-----------|-------------|--------|
| boost::thread | std::thread | TRIVIAL (drop-in) |
| boost::mutex/cv | std::mutex/cv | TRIVIAL (drop-in) |
| boost::format | sprintf/std::format | TRIVIAL |
| boost::date_time | std::chrono | TRIVIAL |
| boost::process | fork()+execv() | LOW-MEDIUM |
| boost::log | Custom logger (~30 LOC) or spdlog | MEDIUM |
| boost::asio | POSIX sockets (~150 LOC wrapper) | MEDIUM |
| boost::iostreams | std::stringstream | MEDIUM (coupled to serialization) |
| boost::serialization | Custom binary format | HIGH (9 versioned structs, RPC transport) |

**Key serialization structures**: EvalData (v6), EvalResults (v9), AircraftState (v3), ScenarioMetadata (v6), DebugSample, PhysicsTraceEntry, Path (v2), GPBytecode/Header.

**Phased approach**:
- Phase 1 (1-2h): Threading, format, date_time → `CMakeLists: Boost REQUIRED COMPONENTS serialization` only
- Phase 2 (4-6h): Sockets, logging, iostreams
- Phase 3 (8-10h): Custom serialization replacing Boost archives

**Note**: Since this is a fork with no backward compat requirement, Phase 3 can be simplified — no need to read old versioned formats. Write a single clean NN-only binary format.

---

## R5: Source Layout Convention

**Decision**: Standard C++ project layout for standalone autoc.

```
autoc/
├── CMakeLists.txt
├── include/autoc/
│   ├── nn/          # NNGenome, NNPopulation, nn_forward, nn_topology
│   ├── eval/        # Evaluator, simulation interface, aircraft_state
│   └── util/        # Config, S3, logging, data output
├── src/
│   ├── nn/
│   ├── eval/
│   ├── util/
│   ├── autoc.cc     # Main evolution executable
│   └── minisim.cc   # Simulation worker executable
├── tools/
│   ├── nnextractor.cc
│   ├── nn2cpp.cc
│   └── renderer.cc
├── tests/
│   ├── nn_evaluator_tests.cc
│   ├── nn_population_tests.cc
│   ├── nn_serialization_tests.cc
│   └── ...
├── third_party/
│   └── inih/        # Vendored config parser
└── autoc.ini
```

**All file moves via `git mv`** to preserve history.

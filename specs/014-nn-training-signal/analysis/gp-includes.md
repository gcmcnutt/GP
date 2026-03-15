# GP Header Include Analysis for autoc/

Analysis of `gp.h` and `gpconfig.h` usage across all autoc source files.
Generated for spec 014-nn-training-signal to guide libgp decoupling.

## Library Chain

```
autoc_common (static lib) --> links libgp.a
  All executables (autoc, minisim, renderer, gpextractor, nnextractor, etc.) link autoc_common
```

## Direct Includes

### Files that `#include "gp.h"` directly

| File | Also includes gpconfig.h | GP Symbols Used | Replacement Strategy |
|------|--------------------------|-----------------|---------------------|
| `autoc.h` | No | `GPGene`, `GPGene::GPContainer::Nth()`, `GPNode`, `GP`, `GPObject`, `GPUserID`, `GPAdfNodeSet` (extern), class inheritance: `MyGene : GPGene`, `MyGP : GP`, `MyGP::stdFitness`, `MyGP::fitnessValid` | **GP-only** -- core GP tree classes. NN mode does not use these classes. Keep for GP mode, no change needed for NN-only refactor. |
| `autoc-eval.cc` | No | `GPAdfNodeSet` (global instance), `GPRegisterClass()`, `GPContainer`, `GPNode`, `GPNodeSet`, `GPGene`, `GP`, `MyGene`, `MyGP`, `createNodeSet()` builds `GPNodeSet` with `GPNode` objects | **GP-only** -- GP tree node set creation and registration. NN mode does not use GP trees. Keep for GP mode. |
| `autoc.cc` | Yes (`gpconfig.h`) | `GPInit()`, `GPrand()` (extensively for PRNG), `GPExitSystem()`, `GPVariables` (config), `GPPopulation` (class inheritance: `MyPopulation : GPPopulation`), `GPContainer::Nth()`, `GPAdfNodeSet`, full GP evolution loop | **Mixed** -- GP evolution + config + PRNG. See detailed breakdown below. |
| `pathgen.cc` | No | `GPrand()` only (3 calls for random point generation) | **Replace**: `GPrand()` -> `std::mt19937` or passed-in RNG |
| `minisim.cc` | No | `MyGP` (via `autoc.h`), `gp.load()`, `gp.resolveNodeValues()`, `gp.NthMyGene()`, `gp.printOn()`, `GPAdfNodeSet` (via `initializeSimGP()`) | **GP-only** -- GP tree deserialization and evaluation in worker. NN path uses `NNGenome` directly. Keep for GP mode. |
| `renderer.cc` | No | `GP gp; gp.load(inStream);` -- single use to extract fitness from serialized GP data (line 1238) | **GP-only** -- used for GP fitness extraction. NN renderer path could use NN deserialization instead. Low priority. |
| `gp_math_utils.h` | No | `GPrand()` (2 calls: `GPrandGaussian()`, `GPrandDouble()`) | **Replace**: `GPrand()` -> `std::mt19937` or inject RNG |
| `config_manager.h` | Yes (`gpconfig.h`) | `GPVariables` (static member), `GPConfiguration` (static member), `GPConfigVarInformation` (in createConfigArray) | **Replace**: `GPVariables` -> custom config struct, `GPConfiguration` -> inih parser |

### Files that `#include "gpconfig.h"` directly

| File | GP Symbols Used | Replacement Strategy |
|------|-----------------|---------------------|
| `autoc.cc` | `GPConfiguration` constructor (via `ConfigManager`) | **Replace**: inih-based config parser |
| `config_manager.h` | `GPConfiguration`, `GPConfigVarInformation`, `GPVariables` | **Replace**: inih-based config parser + custom config struct |

## Transitive Includes (via autoc.h, config_manager.h, gp_math_utils.h)

These files do NOT directly include GP headers but pull them in transitively.

| File | Includes | GP Symbols Actually Used | Replacement Strategy |
|------|----------|--------------------------|---------------------|
| `config_manager.cc` | `config_manager.h` | `GPVariables`, `GPConfiguration`, `GPConfigVarInformation` (config array definition) | **Replace**: inih + custom struct |
| `gpextractor.cc` | `config_manager.h`, `autoc.h` | `GPInit()`, `GPRegisterKernelClasses()`, `GPRegisterClass()`, `GPAdfNodeSet`, `MyGene`, `MyGP`, `gp.load()`, `gp.resolveNodeValues()` | **GP-only** -- GP tree extraction tool. No change needed. |
| `nnextractor.cc` | `config_manager.h`, `autoc.h` | None directly (dummy `MyGP::evaluate()` stubs for linker) | **Decouple**: should not need GP headers at all. Needs autoc_common refactoring to separate GP and NN link deps. |
| `gp_evaluator_portable.h` | `autoc.h` | None from gp.h (uses `Operators` enum and `AircraftState` from autoc.h) | **Already clean** -- only uses autoc-defined types |
| `threadpool.h` | `autoc.h` | None from gp.h directly | **Already clean** -- transitive only |
| `gp_bytecode.cc` | `autoc.h` | None from gp.h directly (uses `Operators` enum) | **Already clean** -- transitive only |
| `bytecode2cpp.cc` | `autoc.h` | None from gp.h directly | **Already clean** -- transitive only |
| `nn_population.cc` | `gp_math_utils.h` | `GPrand()`, `GPrandGaussian()`, `GPrandDouble()` (via gp_math_utils.h wrappers) | **Replace**: inject `std::mt19937` or abstract RNG interface |
| `nn_evaluator_portable.cc` | `gp_math_utils.h` | `GPrandGaussian()` (Xavier init) | **Replace**: inject RNG for weight initialization |
| `variation_generator.h` | (included by autoc.cc) | `GPrand()` (extensively for variation generation, ~20+ calls) | **Replace**: `GPrand()` -> passed-in RNG or `std::mt19937` |

## Detailed Symbol Breakdown for autoc.cc

autoc.cc is the most complex consumer. Symbols used from GP headers:

### From gp.h:
| Symbol | Usage | NN Mode Needs It? |
|--------|-------|-------------------|
| `GPInit(1, gpSeed)` | PRNG seeding + copyright | Yes (for `GPrand()` PRNG) |
| `GPrand()` | ~15 calls for scenario variation pre-fetch, wind seeds, migration | Yes (single PRNG architecture) |
| `GPExitSystem()` | 1 call (error handling) | Replace with `std::abort()` or throw |
| `GPVariables` | Config struct (PopulationSize, NumberOfGenerations, etc.) | Yes (shared config) |
| `GPPopulation` | `MyPopulation : GPPopulation` -- GP evolution only | No (NN uses `NNPopulation`) |
| `GPContainer::Nth()` | Population member access | No (NN uses `std::vector`) |
| `GPAdfNodeSet` | Node set for GP tree creation | No |
| `GPRegisterClass()` | Class registration for load/save | No |
| `GPObject` | Base class for `MyPopulation` | No |

### From gpconfig.h:
| Symbol | Usage | NN Mode Needs It? |
|--------|-------|-------------------|
| `GPConfiguration` | Config file parser (via ConfigManager) | Replace with inih |
| `GPConfigVarInformation` | Config variable descriptors | Replace with inih |
| `GPDataType` | Enum for config types | Replace with inih |

## Summary: Replacement Strategy by Category

### 1. PRNG (`GPrand`, `GPsrand`, `GPInit`)
**Files affected**: autoc.cc, pathgen.cc, gp_math_utils.h, variation_generator.h, nn_population.cc, nn_evaluator_portable.cc
**Strategy**: Replace with `std::mt19937` seeded at startup. The "single PRNG architecture" currently relies on `GPrand()` consuming values in a deterministic order. A replacement must preserve this property (single global `std::mt19937` instance).

### 2. Config System (`GPConfiguration`, `GPConfigVarInformation`, `GPVariables`)
**Files affected**: config_manager.h, config_manager.cc, autoc.cc
**Strategy**: Replace with inih parser. `GPVariables` fields (PopulationSize, NumberOfGenerations, TournamentSize, etc.) move to a new config struct or merge into `ExtraConfig`.

### 3. GP Tree Classes (`GP`, `GPGene`, `GPNode`, `GPNodeSet`, `GPAdfNodeSet`, `GPPopulation`, `GPContainer`, `GPObject`)
**Files affected**: autoc.h, autoc-eval.cc, autoc.cc, minisim.cc, renderer.cc, gpextractor.cc
**Strategy**: **Keep for GP mode**. These are intrinsic to GP tree evolution and cannot be removed without removing GP support. NN-only builds could ifdef or separate into GP-specific translation units.

### 4. Registration/Serialization (`GPRegisterClass`, `GPRegisterKernelClasses`, `GPCreateRegisteredClassObject`)
**Files affected**: autoc-eval.cc, gpextractor.cc
**Strategy**: **Keep for GP mode**. Only used for GP tree load/save.

### 5. Error Handling (`GPExitSystem`)
**Files affected**: autoc.cc (1 call)
**Strategy**: Replace with `throw` or `std::abort()`.

# Implementation Plan: Variations Redux

**Branch**: `003-variations-redux` | **Date**: 2026-03-09 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/003-variations-redux/spec.md`

## Summary

Enable robust training variations (entry conditions, wind direction, variable rabbit speed) for zero-shot sim-to-real transfer. The key finding from research is that **crrcsim integration is already complete** — all variation offsets are unpacked from ScenarioMetadata and applied in aircraft launch and wind initialization. The remaining work is: (1) tighten fitness distance parameters, (2) progressively enable variation configs, and (3) validate training produces recovery-capable controllers.

## Technical Context

**Language/Version**: C++17 (g++, CMake 3.10+)
**Primary Dependencies**: Eigen3 (vectors, quaternions), Boost (serialization, logging, threads)
**Storage**: N/A (in-memory state, S3 for evolution artifacts)
**Testing**: GoogleTest 1.14.0 (82 existing evaluator tests)
**Target Platform**: Linux ARM64 (NVIDIA DGX Spark), also x86 WSL2
**Project Type**: Multi-repo system (GP + crrcsim)
**Performance Goals**: No hard target — scale scenarios incrementally
**Constraints**: Both repos must use same Boost version for serialization compatibility
**Scale/Scope**: Population 20K, up to 36 scenarios/generation, 200 generations

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Testing-First | PASS | Existing 82 evaluator tests cover sensor/interpolation. Fitness param change is config-only, validated by training runs. No new code paths to unit test. |
| II. Build Stability | PASS | Changes are config (#define) and .ini only for GP repo. crrcsim has no code changes. Both repos will be built and verified. |
| III. Dual-Mode Parity | PASS | Fitness constants are shared between GP tree and bytecode eval paths. No new operators or evaluation logic changes. |

## Project Structure

### Documentation (this feature)

```text
specs/003-variations-redux/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0: crrcsim already integrated finding
├── data-model.md        # Phase 1: entities and config progression
├── quickstart.md        # Phase 1: step-by-step validation guide
└── checklists/
    └── requirements.md  # Spec quality checklist
```

### Source Code (changes)

```text
# GP repo (autoc/)
autoc/autoc.h            # DISTANCE_NORM 5.0→2.0, DISTANCE_POWER 1.5→2.0
autoc/autoc.ini          # Enable variations progressively

# crrcsim repo — NO CODE CHANGES
# All variation application code already exists:
#   src/crrc_main.cpp:251-254 (entry offsets)
#   src/mod_windfield/windfield.cpp:569 (wind offset)
#   src/mod_landscape/crrc_builtin_scenery.cpp:741 (wind offset)
#   src/mod_landscape/model_based_scenery.cpp:780 (wind offset)
#   src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp:517-521 (metadata unpacking)
```

## Complexity Tracking

No constitution violations. The feature is primarily config changes and training validation.

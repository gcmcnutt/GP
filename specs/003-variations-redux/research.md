# Research: Variations Redux

**Date**: 2026-03-09
**Feature**: 003-variations-redux

## Key Finding: crrcsim Integration Already Complete

### Decision: No new crrcsim code needed for entry/wind variations
**Rationale**: Investigation of the crrcsim codebase reveals that VARIATIONS1 integration is fully implemented:

- **ScenarioMetadata unpacking** (`inputdev_autoc.cpp:517-521`): Offset values are extracted from ScenarioMetadata and stored in Global variables when pathSelector changes.
- **Entry variations** (`crrc_main.cpp:251-254`): `initialize_flight_model()` applies heading, roll, pitch, and speed offsets before calling `initAirplaneState()`.
- **Wind direction** (`windfield.cpp:569`, `crrc_builtin_scenery.cpp:741`, `model_based_scenery.cpp:780`): All three scenery implementations add `Global::windDirectionOffset` to the base wind direction.
- **Debug logging** (`inputdev_autoc.cpp:524-528`): Conditional `DETAILED_LOGGING` prints received variation values.
- **Global storage** (`global.h:84-91`, `global.cpp:57-62`): Five double fields with proper defaults (0.0 offsets, 1.0 speed factor).

**Alternatives considered**: None needed — the code exists and follows the VARIATIONS1 spec exactly.

### Decision: FR-005 logging is behind DETAILED_LOGGING define
**Rationale**: The debug logging in `inputdev_autoc.cpp:524-528` is gated by `#ifdef DETAILED_LOGGING`. For production use, this may need to be unconditional or controlled by a runtime flag rather than a compile-time define.
**Action**: Evaluate during P2 validation whether the existing logging is sufficient or needs to be made always-on.

## Fitness Distance Tuning

### Decision: DISTANCE_NORM=2.0, DISTANCE_POWER=2.0
**Rationale**: User chose aggressive tightening. The current NORM=5.0 creates a "good enough" zone at ~15-20m. With NORM=2.0 and POWER=2.0:
- Below 2m: cost compresses (sub-linear penalty)
- At 2m: linear crossover
- Above 2m: quadratic penalty (much steeper than current 1.5-power)
- At 5m (old norm): penalty is (5/2)^2 = 6.25 vs old (5/5)^1.5 = 1.0

This should push tracking distance down significantly, but risks making crashes more costly relative to loose tracking. Monitor crash rate in baseline run.

**Alternatives considered**:
- NORM=3.0, POWER=1.5 (moderate, conservative)
- NORM=3.0, POWER=2.0 (moderate norm, steeper curve)
- NORM=4.0, POWER=2.0 (gentle norm, steeper curve)

## Scenario Scaling Strategy

### Decision: Incremental scaling, no hard throughput target
**Rationale**: Current 20-thread setup on ARM64 DGX Spark runs ~1200 sim/sec. With 36 scenarios per individual, generation time scales roughly linearly. Rather than setting arbitrary throughput targets, scale scenarios incrementally and observe:
1. Start with WindScenarios=1 (baseline with tightened fitness only)
2. Scale to 9 (wind variations enabled)
3. Scale to 36 (full entry+wind variations)

If any step produces unacceptable generation times, reduce scenario count.

## ScenarioMetadata Serialization

### Decision: Trust current compatibility
**Rationale**: Both repos currently use the same Boost version on the same platform (ARM64 Ubuntu). The Boost serialization cross-version issue is tracked in the BACKLOG as a separate migration item. No verification testing needed for this feature.

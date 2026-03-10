# Data Model: Variations Redux

**Date**: 2026-03-09
**Feature**: 003-variations-redux

## Entities

### ScenarioMetadata (existing, no changes)

Already version 5 with variation offset fields. Serialized via Boost binary archives between autoc and crrcsim.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| entryHeadingOffset | double | 0.0 | Radians, offset from config heading |
| entryRollOffset | double | 0.0 | Radians, initial roll attitude |
| entryPitchOffset | double | 0.0 | Radians, initial pitch attitude |
| entrySpeedFactor | double | 1.0 | Multiplier on config launch velocity |
| windDirectionOffset | double | 0.0 | Radians, offset from base wind dir |

### VariationSigmas (existing, no changes)

Config-driven Gaussian distribution widths. Defined in `variation_generator.h`.

| Field | Type | Config Key | Config Unit | Internal Unit |
|-------|------|------------|-------------|---------------|
| headingSigma | double | EntryHeadingSigma | degrees | radians |
| rollSigma | double | EntryRollSigma | degrees | radians |
| pitchSigma | double | EntryPitchSigma | degrees | radians |
| speedSigma | double | EntrySpeedSigma | fraction | fraction |
| windDirectionSigma | double | WindDirectionSigma | degrees | radians |

### Fitness Constants (modified)

| Constant | Current | New | File |
|----------|---------|-----|------|
| DISTANCE_NORM | 5.0 | 2.0 | autoc/autoc.h |
| DISTANCE_POWER | 1.5 | 2.0 | autoc/autoc.h |
| ATTITUDE_NORM | 0.349 | 0.349 (unchanged) | autoc/autoc.h |
| ATTITUDE_POWER | 1.5 | 1.5 (unchanged) | autoc/autoc.h |

### Config Parameters (modified)

| Parameter | Current | Phase 1 (P1) | Phase 2 (P2) | Phase 3 (P3) | File |
|-----------|---------|--------------|--------------|--------------|------|
| EnableEntryVariations | 0 | 0 | 0 | 1 | autoc.ini |
| EnableWindVariations | 0 | 0 | 1 | 1 | autoc.ini |
| WindScenarios | 1 | 1 | 9 | 36 | autoc.ini |
| RabbitSpeedSigma | 0.0 | 0.0 | 0.0 | 0.0 (P4: 2.0) | autoc.ini |
| DISTANCE_NORM | 5.0 | 2.0 | 2.0 | 2.0 | autoc.h |
| DISTANCE_POWER | 1.5 | 2.0 | 2.0 | 2.0 | autoc.h |

## State Transitions

No new state transitions. The existing scenario lifecycle is:
1. autoc generates VariationOffsets from seed + sigmas
2. Offsets stored in ScenarioMetadata
3. ScenarioMetadata serialized to crrcsim via Boost RPC
4. crrcsim unpacks to Global:: variables on path change
5. `Simulation->reset()` triggers `initialize_flight_model()` (reads entry offsets) and wind init (reads wind offset)

## Data Flow (unchanged)

```
autoc.ini sigmas → VariationSigmas::fromDegrees() → generateVariations(seed, sigmas)
    → VariationOffsets → ScenarioMetadata fields → Boost serialize → crrcsim
    → Global:: variables → initialize_flight_model() / getWindComponents()
```

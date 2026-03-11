# Data Model: Entry Position Variations with Intercept-Budget Fitness Scaling

**Feature**: 005-entry-fitness-ramp
**Date**: 2026-03-10

## Modified Entities

### VariationSigmas (variation_generator.h)

| Field | Type | Description | Change |
|-------|------|-------------|--------|
| headingSigma | double | radians | existing |
| rollSigma | double | radians | existing |
| pitchSigma | double | radians | existing |
| speedSigma | double | fraction | existing |
| windDirectionSigma | double | radians | existing |
| **positionRadiusSigma** | **double** | **meters (horizontal radius)** | **NEW** |
| **positionAltSigma** | **double** | **meters (vertical/Down)** | **NEW** |

### VariationOffsets (variation_generator.h)

| Field | Type | Description | Change |
|-------|------|-------------|--------|
| entryHeadingOffset | double | radians | existing |
| entryRollOffset | double | radians | existing |
| entryPitchOffset | double | radians | existing |
| entrySpeedFactor | double | multiplier | existing |
| windDirectionOffset | double | radians | existing |
| **entryNorthOffset** | **double** | **meters (NED North)** | **NEW** |
| **entryEastOffset** | **double** | **meters (NED East)** | **NEW** |
| **entryAltOffset** | **double** | **meters (NED Down, negative=up)** | **NEW** |

### ScenarioMetadata (minisim.h)

| Field | Type | Description | Change |
|-------|------|-------------|--------|
| pathVariantIndex | int | -1 = unset | existing |
| windVariantIndex | int | -1 = unset | existing |
| windSeed | unsigned int | PRNG seed | existing |
| scenarioSequence | uint64_t | ordering | existing |
| bakeoffSequence | uint64_t | ordering | existing |
| enableDeterministicLogging | bool | debug flag | existing |
| entryHeadingOffset | double | radians | existing |
| entryRollOffset | double | radians | existing |
| entryPitchOffset | double | radians | existing |
| entrySpeedFactor | double | multiplier | existing |
| windDirectionOffset | double | radians | existing |
| **entryNorthOffset** | **double** | **meters** | **NEW (version 6)** |
| **entryEastOffset** | **double** | **meters** | **NEW (version 6)** |
| **entryAltOffset** | **double** | **meters (NED Down)** | **NEW (version 6)** |

**Serialization**: Version bump from 5 to 6. Loading version ≤5 sets position offsets to 0.0 (backward compatible).

### ExtraConfig (autoc.h)

| Field | Type | Default | Description | Change |
|-------|------|---------|-------------|--------|
| entryHeadingSigma | double | 45.0 | degrees | existing |
| entryRollSigma | double | 22.5 | degrees | existing |
| entryPitchSigma | double | 7.5 | degrees | existing |
| entrySpeedSigma | double | 0.1 | fraction | existing |
| windDirectionSigma | double | 45.0 | degrees | existing |
| **entryPositionRadiusSigma** | **double** | **0.0** | **meters, horizontal radius (0=disabled)** | **NEW** |
| **entryPositionAltSigma** | **double** | **0.0** | **meters, vertical offset (0=disabled)** | **NEW** |

## New Entities

### Intercept Budget Constants (autoc.h)

| Constant | Type | Default | Description |
|----------|------|---------|-------------|
| INTERCEPT_SCALE_FLOOR | double | 0.1 | Minimum scale factor at t=0 |
| INTERCEPT_SCALE_CEILING | double | 1.0 | Maximum scale factor (full penalty) |
| INTERCEPT_BUDGET_MAX | double | 15.0 | Maximum budget cap in seconds |
| INTERCEPT_TURN_RATE | double | π/4 | Estimated turn rate (rad/s) for budget calc |
| ENTRY_SAFE_RADIUS | gp_scalar | ~55.0 | Max entry radius (SIM_PATH_RADIUS_LIMIT - 15m margin) |
| ENTRY_SAFE_ALT_MIN | gp_scalar | ~-22.0 | Shallowest safe altitude (SIM_MIN_ELEVATION - 15m NED) |
| ENTRY_SAFE_ALT_MAX | gp_scalar | ~-105.0 | Deepest safe altitude (SIM_MAX_ELEVATION + 15m NED) |

### Intercept Budget (computed, not stored)

Computed once per scenario at start of fitness evaluation. Not serialized.

| Field | Type | Description |
|-------|------|-------------|
| budget_seconds | double | Estimated time to intercept |
| Derived from | — | displacement, heading offset, aircraft speed, rabbit speed |

### Intercept Scale Factor (computed per-step, not stored)

| Field | Type | Range | Description |
|-------|------|-------|-------------|
| scale | gp_scalar | [FLOOR, CEILING] | `floor + (ceiling - floor) * min(1, (step_time / budget))²` |

## Data Flow

```
startup:
  config_manager reads EntryPositionSigma from autoc.ini
  → stored in ExtraConfig.entryPositionSigma

prefetch (per wind variant):
  generateVariationsFromGPrand(sigmas)
  → generates cylindrical position: Gaussian radius + uniform angle → N/E, Gaussian altitude
  → clamps to safe arena bounds (ENTRY_SAFE_RADIUS, ENTRY_SAFE_ALT_MIN/MAX)
  → produces VariationOffsets with entryNorthOffset, entryEastOffset, entryAltOffset
  → stored in gScenarioVariations table

per-generation, per-scenario:
  populateVariationOffsets(meta)
  → copies offsets into ScenarioMetadata (with RAMP_LANDSCAPE scaling)
  → ScenarioMetadata serialized to minisim/crrcsim via RPC

simulator (minisim or crrcsim):
  → applies position offset to initial aircraft position
  → runs simulation, returns AircraftState sequence

fitness computation (autoc.cc):
  → computes interceptBudget from |aircraftState[0].pos - path[0].start| + heading offset
  → per-step: interceptScale = f(step_time / interceptBudget)
  → distance_sum += pow(interceptScale * distance / NORM, POWER)
  → attitude_sum += pow(interceptScale * attitudeDelta / NORM, POWER)
```

## CRRCSim Globals (global.h)

| Field | Type | Default | Description | Change |
|-------|------|---------|-------------|--------|
| entryHeadingOffset | double | 0.0 | radians | existing |
| entryRollOffset | double | 0.0 | radians | existing |
| entryPitchOffset | double | 0.0 | radians | existing |
| entrySpeedFactor | double | 1.0 | multiplier | existing |
| **entryNorthOffset** | **double** | **0.0** | **meters** | **NEW** |
| **entryEastOffset** | **double** | **0.0** | **meters** | **NEW** |
| **entryAltOffset** | **double** | **0.0** | **meters (NED Down)** | **NEW** |

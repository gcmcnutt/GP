# Research: Entry Position Variations with Intercept-Budget Fitness Scaling

**Feature**: 005-entry-fitness-ramp
**Date**: 2026-03-10

## R1: Intercept Budget Estimation Approach

**Decision**: Simple geometric estimation using distance, heading offset, aircraft speed, and rabbit speed.

**Formula**:
```
turn_time = |heading_offset| / turn_rate_estimate  (turn_rate_estimate ~ 45°/s = π/4 rad/s)
closure_distance = sqrt(northOffset² + eastOffset²)
closure_time = closure_distance / aircraft_speed
rabbit_compensation = closure_time * (rabbit_speed / aircraft_speed)
budget = turn_time + closure_time + rabbit_compensation
budget = clamp(budget, 0, INTERCEPT_BUDGET_MAX)
```

**Rationale**: The user specified "1-digit accuracy" and "crude hacktor" estimation. This geometric approach captures the main contributors (turn time + straight-line closure + rabbit motion) without simulating. Overestimates are benign (gentler ramp), underestimates are unlikely since we're adding components conservatively.

**Alternatives considered**:
- Energy-based estimation (too complex for marginal accuracy gain)
- Look-up table from pre-simulated intercepts (over-engineered)
- Fixed budget per sigma level (loses per-scenario adaptation)

## R2: Scaling Function Shape

**Decision**: Quadratic ramp with configurable floor and ceiling: `scale = floor + (ceiling - floor) * min(1, (t / budget))²`

**Rationale**: Quadratic `t²` provides:
- At t=0: scale = floor (small but nonzero signal)
- At t=budget/2: scale = floor + 0.25 * (ceiling - floor) (~0.325 at floor=0.1)
- At t=budget: scale = ceiling (full penalty)
- Smooth, monotonically increasing, easy to reason about
- Before the power function: `pow(scale * distance / DISTANCE_NORM, DISTANCE_POWER)` — the power function amplifies the scaling effect further

**Alternatives considered**:
- Sigmoid: smooth but has a plateau early that reduces gradient signal
- Linear: too much penalty too early in intercept
- Cubic: too aggressive suppression early, not enough middle-range signal
- Square root: wrong direction — too much penalty early

## R3: Position Offset Implementation Pattern

**Decision**: Follow exact pattern of existing entry attitude offsets in variation_generator.h → ScenarioMetadata → crrcsim Global → crrc_main.cpp. Use cylindrical coordinate generation to match the cylindrical arena.

**Rationale**: The arena is a cylinder (radius 70m, altitude -7m to -120m NED, centered at origin at SIM_INITIAL_ALTITUDE). Cylindrical generation (Gaussian radius + uniform angle + Gaussian altitude) is natural for this geometry and avoids corner effects that a Cartesian Gaussian would have.

**Arena geometry** (from aircraft_state.h):
```
SIM_PATH_RADIUS_LIMIT = 70m    (crash if exceeded)
SIM_MIN_ELEVATION = -7m        (crash: too low/ground in NED)
SIM_MAX_ELEVATION = -120m      (crash: too high in NED)
SIM_INITIAL_ALTITUDE = -25m    (path origin Z)
SIM_INITIAL_LOCATION_DITHER = 30m  (path waypoint spread)
```

**Safe entry envelope** (~15m margin from crash bounds):
- Radius: 0 to ~55m from origin
- Altitude: ~-22m to ~-105m NED

**Key implementation points**:
- `VariationSigmas`: Add `positionRadiusSigma` (meters), `positionAltSigma` (meters)
- `VariationOffsets`: Add `entryNorthOffset`, `entryEastOffset`, `entryAltOffset` (meters)
- `ScenarioMetadata`: Add same three fields, bump to version 6
- `generateVariationsFromGPrand()`: Generate Gaussian radius + uniform angle → N/E, plus Gaussian altitude offset. Clamp to safe bounds.
- crrcsim: `posX += entryNorthOffset`, `posY += entryEastOffset`, altitude offset in `crrc_main.cpp`
- minisim: Add to `initialPosition` before creating AircraftState

**Position generation**:
```
radius = |gaussian(positionRadiusSigma)|    // half-normal: always positive
angle = uniform(0, 2π)
entryNorthOffset = radius * cos(angle)
entryEastOffset = radius * sin(angle)
entryAltOffset = gaussian(positionAltSigma)

// Clamp to safe bounds
total_radius = sqrt(north² + east²)
if total_radius > ENTRY_SAFE_RADIUS: scale down
clamp altitude to [ENTRY_SAFE_ALT_MAX, ENTRY_SAFE_ALT_MIN]
```

**Future expansion**: Existing attitude sigmas (heading=45°, roll=22.5°, pitch=7.5°) can be widened toward all-attitude (heading=180°, roll=180°, pitch=90°) as the training ramp matures. The RAMP_LANDSCAPE infrastructure already scales these from 0→full over generations.

## R4: Scaling Application in Fitness Loop

**Decision**: Apply scaling to both `distance` and `attitude_delta` before their respective `pow()` calls in the per-step fitness accumulation loop.

**Current code** (autoc.cc:1192-1194):
```cpp
distance_sum += pow(distance / DISTANCE_NORM, DISTANCE_POWER);
attitude_sum += pow(attitude_delta / ATTITUDE_NORM, ATTITUDE_POWER);
```

**New code**:
```cpp
gp_scalar intercept_scale = computeInterceptScale(stepTime, interceptBudget);
distance_sum += pow(intercept_scale * distance / DISTANCE_NORM, DISTANCE_POWER);
attitude_sum += pow(intercept_scale * attitude_delta / ATTITUDE_NORM, ATTITUDE_POWER);
```

**Rationale**: Scaling before the power function means the effect is amplified — at scale=0.1, a 50m distance becomes effectively 5m, which is at the DISTANCE_NORM crossover. This is the desired behavior: during intercept, even large distances produce moderate penalties rather than crushing ones.

## R5: Backward Compatibility

**Decision**: When `EntryPositionSigma=0` (default), no position offsets are generated and intercept budget computes to near-zero (degenerate case), resulting in immediate full penalty from step 1 — identical to current behavior.

**Verification**: The scaling function with budget=0 produces `scale = floor + (ceiling - floor) * min(1, (t/0.001))² = ceiling` for all t > 0, which is 1.0 (full penalty). Test this explicitly.

## R6: Intercept Budget in Fitness Computation Context

**Decision**: Intercept budget is computed at the start of each scenario's fitness evaluation in `autoc.cc`, using the ScenarioMetadata entry offsets and initial path point position.

**Note**: The hardcoded initial position in minisim.cc `(-2.19, 5.49, -36.93)` is a legacy artifact matching crrcsim's field start position. The canonical test origin is `(0, 0, SIM_INITIAL_ALTITUDE)` = `(0, 0, -25)`. Entry position offsets should be relative to this test origin. This cleanup should be part of this feature.

**crrcsim start position**: Comes from `autoc_config.xml` `<start position="field" />` which resolves to the scenery start position via `Global::scenery->getStartPosition()`. Position is computed as `posX`/`posY`/`Altitude` in `crrc_main.cpp:210-245`. Entry position offsets should be added after the VARIATIONS1 attitude block (line 254):
```cpp
posX += Global::entryNorthOffset;
posY += Global::entryEastOffset;
Altitude += Global::entryAltOffset;
```
The NED→crrcsim coordinate mapping needs verification but appears to be: North→posX (with sign flip via `-player_pos(2)`), East→posY.

**Available inputs at fitness computation time**:
- `path[0].start` — path start position (known)
- `aircraftState[0].getPosition()` — initial aircraft position (includes offset)
- `meta.entryHeadingOffset` — heading offset (known)
- `SIM_INITIAL_VELOCITY` — aircraft speed (constant)
- `gRabbitSpeedConfig.nominal` — nominal rabbit speed (known)
- Position offset is implied by `|aircraftState[0].getPosition() - path[0].start|`

**Note**: The budget doesn't need the raw position sigma — it uses the actual displacement for this specific scenario, which is deterministic.

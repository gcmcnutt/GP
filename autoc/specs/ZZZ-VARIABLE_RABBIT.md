# VARIABLE RABBIT SPEED

## Overview

The "rabbit" is the reference target point the aircraft follows along the generated path. Currently, the rabbit moves at a **fixed constant velocity** of 16.0 m/s (`SIM_RABBIT_VELOCITY`). This constant speed has trickled into several simplifications throughout the codebase.

This specification proposes **variable rabbit speed** to improve GP training robustness. The rabbit's speed would vary pseudo-randomly with:
- Gaussian distribution around a nominal speed (e.g., 10-20 m/s around 15 m/s)
- Smooth variations over random time intervals (0.5-5 second cycles)
- Deterministic generation from the existing `WindSeedBase` PRNG

This trains the GP to handle a following target that accelerates, decelerates, and varies its pace - more realistic than a constant-speed reference.

## Related Documentation

- [ZZZ-VARIATIONS1.md](ZZZ-VARIATIONS1.md) - Entry/wind variations (similar PRNG pattern)
- [TEMPORAL_STATE.md](TEMPORAL_STATE.md) - Error history nodes (related: velocity history)
- [variation_generator.h](../variation_generator.h) - Existing PRNG implementation

## Current Implementation

### Path Generation (pathgen.h)

All path points compute `simTimeMsec` from distance traveled at fixed velocity:

```cpp
// pathgen.h line 90 (GenerateRandom), similar in all generators
gp_scalar simTimeMsec = (odometer / SIM_RABBIT_VELOCITY) * static_cast<gp_scalar>(1000.0f);
```

The constant is defined in `aircraft_state.h:36`:
```cpp
#define SIM_RABBIT_VELOCITY static_cast<gp_scalar>(16.0f)  // m/s
```

### Path Structure (aircraft_state.h:54-70)

```cpp
class Path {
public:
  gp_vec3 start;              // 3D waypoint position (NED meters)
  gp_vec3 orientation;        // Direction hint (always UnitX currently)
  gp_scalar distanceFromStart; // Accumulated distance along path
  gp_scalar radiansFromStart;  // Accumulated turn angle
  gp_scalar simTimeMsec;       // PRECOMPUTED timestamp at this waypoint
};
```

The `simTimeMsec` field is **baked during path generation**, not computed at runtime.

### Path Consumption

All three consumers use **time-based lookup**:

**minisim.cc (lines 264-267):**
```cpp
while (aircraftState.getThisPathIndex() < path.size() - 2 &&
       (path.at(aircraftState.getThisPathIndex()).simTimeMsec < duration_msec)) {
  aircraftState.setThisPathIndex(aircraftState.getThisPathIndex() + 1);
}
```

**inputdev_autoc.cpp (lines 738-742):**
```cpp
while (pathIndex < path.size() - 2 &&
       (path.at(pathIndex).simTimeMsec < simTimeMsec)) {
  pathIndex++;
}
```

**xiao-gp/msplink.cpp (getRabbitPathIndex, lines 857-874):**
```cpp
int getRabbitPathIndex(unsigned long elapsed_msec) {
  for (size_t i = current_path_index; i < flight_path.size(); i++) {
    if (flight_path[i].simTimeMsec >= elapsed_msec) {
      return (int)i;
    }
  }
  return (int)(flight_path.size() - 1);
}
```

### Lookahead Calculation (aircraft_state.h:187)

The `getPathIndex()` function uses fixed velocity for speculative lookahead:
```cpp
const gp_scalar distanceGoal = distanceSoFar + steps * SIM_RABBIT_VELOCITY * (SIM_TIME_STEP_MSEC / 1000.0f);
```

At 10Hz (100ms steps), this advances ~1.6m per step (16.0 * 0.1).

## Proposed Design

### Speed Variation Model

The rabbit speed varies smoothly using **layered Perlin-style noise** (or similar):

| Parameter | Value | Description |
|-----------|-------|-------------|
| `RabbitSpeedNominal` | 15.0 m/s | Center of Gaussian distribution |
| `RabbitSpeedSigma` | 3.0 m/s | 1σ deviation (68% within 12-18 m/s) |
| `RabbitSpeedMin` | 8.0 m/s | Hard floor (prevent stall-like following) |
| `RabbitSpeedMax` | 25.0 m/s | Hard ceiling (prevent unreachable target) |
| `RabbitSpeedCycleMin` | 0.5 sec | Minimum variation cycle period |
| `RabbitSpeedCycleMax` | 5.0 sec | Maximum variation cycle period |

**NOT physics-based**: The speed variation is random, not derived from energy or flight dynamics. This is intentional - we want to train robustness to arbitrary speed changes.

### Speed Profile Generation

Generate a **speed profile** alongside the path. The profile consists of time-tagged speed values that vary smoothly:

```cpp
struct RabbitSpeedPoint {
    gp_scalar timeMsec;       // Time along simulation (0 to SIM_TOTAL_TIME_MSEC)
    gp_scalar speed;          // Rabbit speed at this time (m/s)
    gp_scalar acceleration;   // Rate of speed change (m/s^2, for interpolation)
};
```

**Generation algorithm:**
1. Start at `t=0` with speed drawn from `Gaussian(nominal, sigma)`, clamped to [min, max]
2. Pick a random cycle duration from uniform `[cycleMin, cycleMax]`
3. At end of cycle, draw new target speed from Gaussian
4. Interpolate smoothly (cosine or cubic) between speeds
5. Repeat until `t > SIM_TOTAL_TIME_MSEC`

This produces ~20-200 speed waypoints per 100-second path, depending on cycle lengths.

### Path Structure (Unchanged)

The Path class stays the same - no new fields needed:

```cpp
class Path {
public:
  gp_vec3 start;
  gp_vec3 orientation;
  gp_scalar distanceFromStart;
  gp_scalar radiansFromStart;
  gp_scalar simTimeMsec;       // Now computed from variable speed profile
};
```

The variable speed is **baked into `simTimeMsec`** during path generation:

```cpp
// During path generation
gp_scalar simTimeMsec = 0.0f;
gp_scalar lastDistance = 0.0f;

for (each path point at distance d) {
    gp_scalar speed = getSpeedAtTime(simTimeMsec);  // From speed profile
    gp_scalar dt = (d - lastDistance) / speed * 1000.0f;
    simTimeMsec += dt;
    lastDistance = d;

    pathPoint.simTimeMsec = simTimeMsec;  // Variable speed encoded in timing
}
```

No `rabbitSpeed` field needed - the timing implicitly encodes speed. If we later need explicit speed access (e.g., `GETRABBITSPEED` terminal), we can compute it from adjacent points: `speed = Δdistance / Δtime`.

### Speed Profile Generator

Add to `variation_generator.h` or new `rabbit_speed_generator.h`:

```cpp
struct RabbitSpeedConfig {
    double nominal;      // Center speed (m/s)
    double sigma;        // 1σ deviation (m/s)
    double minSpeed;     // Floor (m/s)
    double maxSpeed;     // Ceiling (m/s)
    double cycleMin;     // Min cycle duration (seconds)
    double cycleMax;     // Max cycle duration (seconds)
};

// Generate speed profile for given duration
// Continues PRNG sequence from windSeed (call after generateVariations())
std::vector<RabbitSpeedPoint> generateSpeedProfile(
    unsigned int& seed,  // reference - continues existing PRNG sequence
    const RabbitSpeedConfig& config,
    double totalDurationSec
);

// Query speed at arbitrary time (interpolates between waypoints)
double getSpeedAtTime(
    const std::vector<RabbitSpeedPoint>& profile,
    double timeSec
);
```

### Configuration (autoc.ini)

Add rabbit speed parameters to existing variation section:

```ini
# =============================================================================
# SCENARIO VARIATIONS (extends existing section)
# =============================================================================
# ... existing WindScenarios, WindSeedBase, entry/wind sigmas ...

# Rabbit speed distribution (Gaussian around nominal, clamped to [min, max])
# Set Sigma=0 for constant speed at Nominal
RabbitSpeedNominal           = 15.0    # m/s - center of distribution
RabbitSpeedSigma             = 3.0     # m/s - 1σ deviation (0 = constant speed)
RabbitSpeedMin               = 8.0     # m/s - hard floor
RabbitSpeedMax               = 25.0    # m/s - hard ceiling

# Variation timing (random cycle lengths within this range)
# Ignored if Sigma=0
RabbitSpeedCycleMin          = 0.5     # seconds
RabbitSpeedCycleMax          = 5.0     # seconds
```

**Uses existing WindSeedBase PRNG** - no separate seed. The rabbit speed profile is generated from the same `windSeed` that drives entry/wind variations, extending the existing `generateVariations()` pattern in `variation_generator.h`.

**Constant speed mode**: Set `RabbitSpeedSigma = 0` and `RabbitSpeedNominal = 16.0` to replicate current fixed-speed behavior.

## Implementation Plan

### Phase 1: Speed Profile Generator (autoc only)

**Files to modify:**
| File | Changes |
|------|---------|
| `variation_generator.h` | Add speed profile generation (extend existing PRNG) |
| `autoc.h` | Add RabbitSpeedConfig to ExtraConfig |
| `config_manager.cc` | Parse new config parameters |

**Extend variation_generator.h** - add speed profile generation using the same PRNG pattern:

```cpp
// Add to variation_generator.h

struct RabbitSpeedConfig {
    double nominal;    // m/s
    double sigma;      // m/s (0 = constant speed)
    double minSpeed;   // m/s
    double maxSpeed;   // m/s
    double cycleMin;   // seconds
    double cycleMax;   // seconds
};

struct RabbitSpeedPoint {
    double timeSec;
    double speed;
};

// Generate speed profile using same LCG PRNG as generateVariations()
// Called after generateVariations() consumes some PRNG state - this continues the sequence
inline std::vector<RabbitSpeedPoint> generateSpeedProfile(
    unsigned int& seed,  // NOTE: reference - continues PRNG sequence
    const RabbitSpeedConfig& cfg,
    double totalDurationSec
) {
    std::vector<RabbitSpeedPoint> profile;

    // Constant speed mode
    if (cfg.sigma <= 0.0) {
        profile.push_back({0.0, cfg.nominal});
        profile.push_back({totalDurationSec, cfg.nominal});
        return profile;
    }

    // Reuse same LCG PRNG pattern
    auto nextDouble = [&seed]() -> double {
        seed = seed * 1103515245 + 12345;
        return static_cast<double>((seed >> 16) & 0x7FFF) / 32768.0;
    };

    auto gaussian = [&nextDouble](double mean, double sigma) -> double {
        double u1 = nextDouble() * 0.999 + 0.001;
        double u2 = nextDouble();
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        return mean + z * sigma;
    };

    auto clampSpeed = [&cfg](double s) -> double {
        return (s < cfg.minSpeed) ? cfg.minSpeed :
               (s > cfg.maxSpeed) ? cfg.maxSpeed : s;
    };

    double t = 0.0;
    double currentSpeed = clampSpeed(gaussian(cfg.nominal, cfg.sigma));
    profile.push_back({t, currentSpeed});

    while (t < totalDurationSec) {
        double cycleDuration = cfg.cycleMin + nextDouble() * (cfg.cycleMax - cfg.cycleMin);
        double targetSpeed = clampSpeed(gaussian(cfg.nominal, cfg.sigma));

        // Cosine-eased interpolation points (~10Hz)
        int numSteps = static_cast<int>(cycleDuration / 0.1);
        if (numSteps < 2) numSteps = 2;

        for (int i = 1; i <= numSteps && (t + i * cycleDuration / numSteps) <= totalDurationSec; i++) {
            double frac = static_cast<double>(i) / numSteps;
            double easedFrac = 0.5 * (1.0 - cos(M_PI * frac));
            profile.push_back({t + frac * cycleDuration,
                               currentSpeed + easedFrac * (targetSpeed - currentSpeed)});
        }

        t += cycleDuration;
        currentSpeed = targetSpeed;
    }

    return profile;
}

inline double getSpeedAtTime(const std::vector<RabbitSpeedPoint>& profile, double timeSec) {
    if (profile.empty()) return 16.0;
    if (timeSec <= profile.front().timeSec) return profile.front().speed;
    if (timeSec >= profile.back().timeSec) return profile.back().speed;

    // Binary search + linear interpolation
    size_t lo = 0, hi = profile.size() - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (profile[mid].timeSec <= timeSec) lo = mid;
        else hi = mid;
    }
    double frac = (timeSec - profile[lo].timeSec) / (profile[hi].timeSec - profile[lo].timeSec);
    return profile[lo].speed + frac * (profile[hi].speed - profile[lo].speed);
}
```

### Phase 2: Path Generation Integration

**Files to modify:**
| File | Changes |
|------|---------|
| `pathgen.h` | Integrate speed profile into path generation |
| `minisim.h` | Add speed config to ScenarioMetadata |

**Key change in pathgen.h** - replace fixed velocity calculation:

```cpp
// OLD (fixed velocity):
gp_scalar simTimeMsec = (odometer / SIM_RABBIT_VELOCITY) * 1000.0f;

// NEW (variable velocity):
// speedProfile passed to generator, or generated from scenario seed
gp_scalar speed = getSpeedAtTime(speedProfile, simTimeMsec / 1000.0);
gp_scalar dt = segmentDistance / speed * 1000.0f;
simTimeMsec += dt;
// No rabbitSpeed field - timing encodes speed implicitly
```

### Phase 3: Time-Domain Lookahead (Critical Refactor)

> **NOTE**: This refactor should be done regardless of variable rabbit speed. The current hybrid approach (time→distance→search) is unnecessarily complex. Pure time-based lookahead is cleaner and semantically correct: "1 step ahead" means "0.1 seconds ahead", not "1.6 meters ahead".

**Important insight**: The current `getPathIndex()` converts time→distance using fixed velocity, then searches by distance. This is a hybrid approach that breaks with variable rabbit speed.

**Current implementation (aircraft_state.h:187):**
```cpp
const gp_scalar distanceGoal = distanceSoFar + steps * SIM_RABBIT_VELOCITY * (SIM_TIME_STEP_MSEC / 1000.0f);
// Then search path by distance...
```

**Problem**: With variable speed, the time→distance mapping is wrong.

**Solution**: Pure time-based lookahead. Path points already have `simTimeMsec` - search by time directly:

```cpp
// NEW (pure time-domain lookahead):
inline int getPathIndex(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {
    if (std::isnan(arg)) {
        return pathProvider.getCurrentIndex();
    }

    int steps = CLAMP_DEF((int)arg, -5, 5);
    int currentStep = CLAMP_DEF(pathProvider.getCurrentIndex(), 0, pathProvider.getPathSize() - 1);

    // Time-based lookahead: N steps = N * 100ms
    const gp_scalar currentTimeMsec = pathProvider.getPath(currentStep).simTimeMsec;
    const gp_scalar timeGoalMsec = currentTimeMsec + steps * SIM_TIME_STEP_MSEC;

    constexpr gp_scalar kEps = static_cast<gp_scalar>(0.1f);  // 0.1ms epsilon

    if (steps > 0) {
        while (currentStep < pathProvider.getPathSize() - 1 &&
               pathProvider.getPath(currentStep).simTimeMsec + kEps < timeGoalMsec) {
            currentStep++;
        }
    } else if (steps < 0) {
        while (currentStep > 0 &&
               pathProvider.getPath(currentStep).simTimeMsec - kEps > timeGoalMsec) {
            currentStep--;
        }
    }

    return currentStep;
}
```

**This is cleaner and works correctly regardless of variable speed mode.**

### Phase 3a: Path Resolution (Nyquist Requirement)

For time-based lookahead to work, path points must have sufficient **temporal resolution**.

**Current state:**
- Straight segments: 1.0m spacing
- At 16 m/s → 62.5ms/point (16 Hz) - acceptable for 10 Hz control
- At 8 m/s → 125ms/point (8 Hz) - **below Nyquist for 10 Hz control!**

**Requirement**: For 10 Hz control, need ≥20 Hz path resolution (Nyquist).

**At minimum rabbit speed (8 m/s):**
- Need 50ms/point max
- 50ms × 8 m/s = 0.4m spacing

**Solution**: Reduce path point spacing from 1.0m to **0.4m**:

```cpp
// pathgen.h - addStraightSegment
const gp_scalar step = 0.4f;  // 0.4m spacing (ensures ≥20Hz at min rabbit speed)

// pathgen.h - addHorizontalTurn
const gp_scalar step = 0.02f;  // ~1° spacing (finer angular resolution)
```

**Memory impact:**
- 100 seconds × 20 m/s average = 2000m path
- At 0.4m spacing: 5000 points
- At ~32 bytes/point: 160KB (was 64KB at 1.0m spacing)
- Still acceptable for desktop; embedded needs guard update

**Path segment limits to update:**

| Location | Current | New | Notes |
|----------|---------|-----|-------|
| `xiao-gp/include/embedded_pathgen_selector.h:15` | `MAX_EMBEDDED_PATH_SEGMENTS = 400` | `1000` | Path 4 needs ~357 at 1.0m → ~892 at 0.4m |
| `minisim.cc` | Dynamic vector | No change | Uses `std::vector`, no hard limit |
| `inputdev_autoc.cpp` | Dynamic vector | No change | Uses `std::vector`, no hard limit |

**Alternative**: Keep 1.0m spacing but interpolate between points during lookahead. This trades runtime computation for memory.

### Phase 3b: Consumer Updates

With time-domain lookahead in place, consumers need minimal changes:

**minisim.cc** - Path index lookup unchanged (already time-based). `getPathIndex()` now works correctly.

**inputdev_autoc.cpp** - Same, no changes needed beyond using updated `getPathIndex()`.

**xiao-gp/msplink.cpp** - `getRabbitPathIndex()` unchanged. `getPathIndex()` works correctly.

### Phase 4: Embedded Support

**embedded_pathgen_selector.h** - Update path generation for embedded use:
- Speed profile generation must work on XIAO (limited RAM)
- Speed is baked into `simTimeMsec` during generation - no extra memory per point
- Speed profile itself is ~100-200 floats (~400-800 bytes) - acceptable

### Phase 5: Fitness Considerations

With variable rabbit speed, the GP may need to:
1. Manage throttle to match varying target speed
2. Handle acceleration/deceleration smoothly
3. Not over-penalize when rabbit speeds up (aircraft temporarily falls behind)

**Potential fitness adjustments:**
```cpp
// Current asymmetric penalty (harsher for overshooting):
#define WAYPOINT_AHEAD_POWER 1.8
#define WAYPOINT_BEHIND_POWER 1.2

// With variable rabbit, may want to relax "behind" penalty when rabbit accelerating
// This can be revisited after initial testing
```

## Files Changed Summary

### Phase 0: Time-Domain Lookahead (do first, independent of variable speed)

| File | Changes |
|------|---------|
| `aircraft_state.h` | Refactor `getPathIndex()` to search by `simTimeMsec` instead of distance |
| `pathgen.h` | Reduce spacing from 1.0m to 0.4m for Nyquist compliance |
| `xiao-gp/.../embedded_pathgen_selector.h` | Increase `MAX_EMBEDDED_PATH_SEGMENTS` from 400 to 1000 |

### Phase 1-5: Variable Rabbit Speed

| File | Changes |
|------|---------|
| `autoc.ini` | Add rabbit speed parameters to variation section |
| `autoc.h` | Add `RabbitSpeedConfig` to `ExtraConfig` |
| `config_manager.cc` | Parse new parameters |
| `variation_generator.h` | Add `generateSpeedProfile()` using existing PRNG |
| `pathgen.h` | Integrate variable speed into `simTimeMsec` calculation |
| `embedded_pathgen_selector.h` | Update for variable speed support |

## Backward Compatibility

- **Constant speed**: Set `RabbitSpeedSigma = 0`, `RabbitSpeedNominal = 16.0` for legacy behavior
- **No Path struct changes** - `simTimeMsec` field usage unchanged, just computed differently
- **Missing config params**: Default to `Sigma=0, Nominal=16.0` (legacy constant speed)

## Testing Strategy

1. **Unit tests**: Verify speed profile generation stays within bounds
2. **Path generation test**: Confirm simTimeMsec accumulates correctly with variable speed
3. **Visual inspection**: Render paths with speed colormap, verify smooth variation
4. **Regression test**: Run with `RabbitSpeedSigma=0, RabbitSpeedNominal=16.0`, verify identical to current behavior
5. **Evolution test**: Run GP training, verify fitness improves and throttle management emerges

## Success Criteria

1. Speed profiles are smooth (no sudden jumps)
2. Path timestamps are consistent with integrated variable speed
3. GP evolves throttle management behaviors
4. Flight test shows improved robustness to varying approach speeds

## Risks

1. **Computation cost**: Speed profile generation adds ~1ms per scenario
   - Mitigation: Generate once per scenario, reuse across population

2. **Memory footprint**: ~4 bytes/pathpoint increase for rabbitSpeed field
   - Mitigation: Acceptable for 1000-point paths (~4KB increase)

3. **Training difficulty**: Variable speed may make training harder
   - Mitigation: Start with small sigma (1.0), increase gradually

4. **Embedded constraints**: XIAO has limited RAM
   - Mitigation: Path already fits; 4 bytes/point is acceptable

## Future Extensions

- **GETRABBITSPEED terminal**: GP can observe current rabbit speed
- **GETRABBITACCEL terminal**: GP can observe rabbit acceleration
- **Physics-based speed**: Vary speed based on path curvature (slower in turns)
- **Adaptive sigma**: Start training with small variations, increase over generations

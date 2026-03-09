# Spec: Path Interpolation & Evaluator Improvements

**Feature ID**: 002-path-interpolation
**Status**: Clarified
**Created**: 2026-03-06
**Clarified**: 2026-03-06

## Overview

This feature fixes the path indexing jitter bug discovered in 001-gp-eval-tests and verifies related evaluator improvements (FITNESS_SIMPLIFY, RAMP_LANDSCAPE) are working correctly.

## Problem Statement

### Primary: Discrete Waypoint Snapping

The current `getPathIndex()` function snaps to discrete waypoints, causing:

1. **Timing jitter overshoot**: Waypoint at t=98ms is skipped when looking for t=100ms because `98 < 100`, jumping to t=200ms instead
2. **Discontinuous sensors**: GETDPHI/GETDTHETA jump discretely between waypoints
3. **Real-time sensitivity**: On xiao-gp, eval loop jitter causes erratic sensor values

**Evidence**: Tests `NavigationOps.EarlyTimestampOvershoot` and `NavigationOps.NonUniformTimestamps` in `gp_evaluator_tests.cc` document this behavior.

### Secondary: Verify Existing Implementations

The following are implemented but need verification:
- **FITNESS_SIMPLIFY**: 2-objective fitness (distance + attitude) with path-relative scaling
- **RAMP_LANDSCAPE**: Gradual variation scaling during training

## Solution

### Path Interpolation

Replace discrete waypoint lookup with time-based position interpolation:

```cpp
gp_vec3 getInterpolatedTargetPosition(
    PathProvider& pathProvider,
    gp_scalar currentTimeMsec,
    gp_scalar offsetSteps
) {
    // 1. Calculate goal time (clamped to ±1 second = ±10 steps)
    gp_scalar clampedSteps = std::clamp(offsetSteps, -10.0f, 10.0f);
    gp_scalar goalTimeMsec = currentTimeMsec + clampedSteps * SIM_TIME_STEP_MSEC;

    // 2. Binary search for bracketing waypoints
    // Find i where path[i].simTimeMsec <= goalTime < path[i+1].simTimeMsec

    // 3. Linear interpolate position
    gp_scalar frac = (goalTimeMsec - p0.simTimeMsec) / (p1.simTimeMsec - p0.simTimeMsec);
    return p0.start + frac * (p1.start - p0.start);
}
```

**Boundary Handling**: Clamp offset to ±1 second (±10 steps at 100ms/step). Future estimates beyond this range will be fuzzed in a later feature (see BACKLOG.md "Error Cone for Future Path Points").

### Updated Navigation Functions

- `executeGetDPhi()` - use interpolated position
- `executeGetDTheta()` - use interpolated position
- `executeGetDTarget()` - use interpolated position

### Remove Discrete Function

- Remove `getPathIndex()` from `aircraft_state.h`

## FITNESS_SIMPLIFY Design

### Objectives

Two-objective fitness replacing previous 4-objective approach:

| Objective | Units | Formula |
|-----------|-------|---------|
| Distance to rabbit | meters | `pow(distance, DISTANCE_POWER)` |
| Attitude delta | scaled meters | `attitude_change * attitude_scale` |

### Key Parameters

- **DISTANCE_POWER = 1.2**: Nonlinear penalty keeps small excursions expensive
- **attitude_scale**: Path-relative scaling normalizes radians to effective meters

### Path-Relative Scaling

```cpp
// Per-path: compute scale from THIS path's geometry
gp_scalar path_distance = path.back().distanceFromStart;  // total meters
gp_scalar path_turn_rad = path.back().radiansFromStart;   // total radians
gp_scalar attitude_scale = path_distance / std::max(path_turn_rad, 1.0f);

// Per-step fitness
gp_scalar distance = (aircraftPosition - rabbitPosition).norm();
gp_scalar attitude_change = fabs(roll - prev_roll) + fabs(pitch - prev_pitch);
gp_scalar step_fitness = pow(distance, DISTANCE_POWER) + attitude_change * attitude_scale;
```

**Rationale**: "1 radian of excess turning costs same as `attitude_scale` meters of distance error" - makes objectives comparable.

### Removed Metrics

- Direction alignment (redundant - bad direction increases distance)
- Energy deviation (loophole: low+fast = same energy as high+slow)
- Per-objective weight constants (MOVEMENT_DIRECTION_WEIGHT, etc.)

### Kept

- CRASH_COMPLETION_WEIGHT - orthogonal concern for soft lexicographic ordering

### Known Issue: Straight-Line Path Oscillation

When `path_turn_rad` approaches 0 (straight-line paths), the craft tends to oscillate around the path axis. This is a fundamental sensor/control issue beyond just scaling:

- **Symptom**: Craft "wants to turn around the path" even when on course
- **Impact**: Energy waste, unstable tracking on straight segments
- **Root cause**: GETDPHI/GETDTHETA sensors give non-zero values when offset from path, encouraging corrections even when heading is correct
- **Future work**: Consider damping or clamping sensor values when heading alignment is good, or add "on-course" detection to fitness

For this feature, use `max(path_turn_rad, 1.0f)` as a reasonable clamp.

## RAMP_LANDSCAPE Design

### Purpose

Entry and rabbit speed variations corrupt the fitness landscape early in training. Ramping variations from 0 to full allows:

1. Early generations establish good tracking behavior (easy problem)
2. Progressive introduction of harder scenarios
3. Existing skills adapt incrementally to entry recovery

### Configuration

```ini
# Generations per step (0 = disabled, use full sigmas immediately)
VariationRampStep = 5
```

### Scale Calculation

```cpp
double computeVariationScale(int currentGen, int totalGens, int rampStep) {
    if (rampStep <= 0) return 1.0f;  // Disabled

    int stepIndex = currentGen / rampStep;
    int totalSteps = totalGens / rampStep;
    if (totalSteps == 0) return 1.0f;

    return std::min(1.0f, static_cast<double>(stepIndex) / static_cast<double>(totalSteps));
}
```

### Applied Sigmas

Scale multiplies **all** variation sigmas before scenario generation:

- EntryHeadingSigma * scale
- EntryRollSigma * scale
- EntryPitchSigma * scale
- EntrySpeedSigma * scale
- RabbitSpeedSigma * scale
- WindSpeedSigma * scale (if added to scenario payload)
- WindDirectionSigma * scale (if added to scenario payload)

**Note**: Wind parameters may need to be added to the GP→sim scenario payload. Currently entry/rabbit are included but wind may not be. Progressive difficulty also interacts with scenario count.

### Example Behavior

With `VariationRampStep=5` over 100 generations:
- Gen 0-4: scale = 0.00 (nominal entries)
- Gen 5-9: scale = 0.05
- Gen 50-54: scale = 0.50
- Gen 95-99: scale = 0.95

### Logging

```
<info> Gen 10: VariationScale=0.10 (step 2/20)
```

## Scope

### In Scope

| Item | Status | Priority | Work Required |
|------|--------|----------|---------------|
| Path Interpolation | New | **1** | Implement, test (foundation for others) |
| Update gp_evaluator_tests | Existing | **1** | Update for interpolation |
| FITNESS_SIMPLIFY | Implemented | 2 | Verify working (after interpolation) |
| RAMP_LANDSCAPE | Implemented | 2 | Verify working (after interpolation) |

**Implementation Order**: Path interpolation must be completed first - it fixes the jitter bug that would otherwise corrupt fitness signals during verification of the other features.

### Out of Scope

- New GP operators
- Fitness function redesign (FITNESS_SIMPLIFY already done)
- Embedded (xiao-gp) changes (uses same portable evaluator)
- Future estimate fuzzing (separate backlog item)
- Target pose/orientation interpolation (future backlog item - may be needed for pose estimation sensors)

## Files Affected

| File | Changes |
|------|---------|
| `gp_evaluator_portable.h` | Add `getInterpolatedTargetPosition()` |
| `gp_evaluator_portable.cc` | Implement interpolation, update `executeGetD*()` |
| `aircraft_state.h` | Remove `getPathIndex()`, add `getMaxTimeMsec()` to PathProvider |
| `tests/gp_evaluator_tests.cc` | Update navigation tests for interpolated positions |
| `autoc.cc` | Verify FITNESS_SIMPLIFY and RAMP_LANDSCAPE implementation |

## Verification

### Path Interpolation Tests

1. `TEST(NavigationOps, InterpolationMidpoint)` - verify lerp at 50% between waypoints
2. `TEST(NavigationOps, InterpolationContinuity)` - verify no jumps at waypoint boundaries
3. `TEST(NavigationOps, InterpolationJitterRobust)` - verify jittery times work correctly
4. Update existing `EarlyTimestampOvershoot` test to verify fix

### FITNESS_SIMPLIFY Verification

**Phase 1: Code Inspection + Unit Tests**
- Confirm `DISTANCE_POWER = 1.2` is used in fitness calculation
- Confirm `attitude_scale = path_distance / path_turn_rad` is computed per-path
- Verify removed metrics (direction, energy) are not computed
- Add unit test for `computeAttitudeScale()` with various path geometries

**Phase 2: Short Evolution Run**
- Run ~10 generations with small population
- Verify fitness values are reasonable (not NaN, not exploding)
- Confirm 2-objective formula produces monotonic improvement

### RAMP_LANDSCAPE Verification

**Phase 1: Code Inspection + Unit Tests**
- Verify `computeVariationScale()` implementation matches spec
- Add unit test for scale computation at various generation/step values
- Confirm scale is applied to all variation sigmas

**Phase 2: Short Evolution Run**
- Set `VariationRampStep = 5`, run ~20 generations
- Verify log shows "VariationScale" increasing over generations
- Verify gen-0 has scale=0 (no variations)
- Test with `VariationRampStep=0` produces full sigmas immediately

## Success Criteria

1. **Smoothness**: Sensor values change continuously, no discrete jumps
2. **Jitter robustness**: 10ms timing jitter causes <1% sensor value change
3. **Existing tests pass**: All 77 evaluator tests still pass
4. **Build passes**: `cd ~/GP && make` succeeds
5. **FITNESS_SIMPLIFY**: 2-objective fitness produces monotonic improvement
6. **RAMP_LANDSCAPE**: Variation scale ramps correctly per generation

## Clarifications

The following design decisions were clarified during spec review:

1. **Boundary handling**: Clamp interpolation offset to ±1 second (±10 steps). Future fuzzing is a separate backlog item.

2. **Straight-line paths**: Known issue - craft oscillates around path axis due to sensor design. Use `max(path_turn_rad, 1.0f)` as clamp for now; damping/on-course detection is future work.

3. **Wind ramp**: RAMP_LANDSCAPE should apply to all variations including wind. Wind parameters may need to be added to GP→sim scenario payload.

4. **Implementation order**: Path interpolation first (priority 1), then verify FITNESS_SIMPLIFY/RAMP_LANDSCAPE (priority 2).

5. **Orientation interpolation**: Position-only for now. Target pose estimation is a future backlog item.

6. **Verification method**: Both code inspection + unit tests AND short evolution runs for FITNESS_SIMPLIFY/RAMP_LANDSCAPE.

## Known Issues (Mar 2026 — aarch64 migration)

### Float Precision in Interpolator (Non-Determinism Risk)

The path interpolation system has 32-bit float precision issues that may cause subtle
non-determinism across re-evaluations, observed as fitness regressions on the elite
individual (fitness jumps 4x worse then recovers next gen, zero divergences detected).

**Root causes identified:**

1. **`unsigned long` → `float` cast loses precision at high sim times.**
   `aircraftState.getSimTimeMsec()` is `unsigned long` but cast to `gp_scalar` (float,
   ~7 decimal digits) at `gp_evaluator_portable.cc:354`. After ~100s of sim time
   (100,000ms), float cannot distinguish consecutive milliseconds. The binary search
   at line 328 uses exact `<=` comparison — a 1-ULP rounding difference selects a
   different waypoint bracket.

2. **Accumulated float odometer in path generation (`pathgen.h:67-94`).**
   Each waypoint's `simTimeMsec` derives from `odometer / velocity * 1000`, where
   `odometer += newDistance` accumulates across ~50+ loop iterations with transcendental
   `norm()` calls. Deterministic per-seed but fragile to optimizer reordering.

3. **Binary search has no epsilon tolerance (`gp_evaluator_portable.cc:322-333`).**
   Exact float comparison means issues (1) and (2) can select different waypoint pairs
   on re-evaluation, producing different interpolated positions and sensor values.

**Constraint:** `double` is too slow on target embedded hardware (Xiao-GP). Must stay
32-bit float (`gp_scalar`).

**Proposed fix direction (revisit in this feature):**

- The rabbit moves forward only — no backtracking needed
- Search area is small: ±MAX_OFFSET_STEPS (±1 second at 10Hz)
- Waypoint spacing is ~1.6m at 16 m/s
- **Forward-scan from last known index** instead of binary search eliminates float
  comparison sensitivity entirely — just walk forward from `getCurrentIndex()`
- **Integer millisecond timestamps** in path would avoid float time comparison altogether
- Path is generated once per scenario — accumulated error is consistent within one eval,
  so the fix only needs to ensure the *lookup* is deterministic, not the path times

### `-ffast-math` Removed on aarch64

Performance builds originally used `-ffast-math` (5000 sims/sec). Removed because:
- `-ffinite-math-only` optimizes away `std::isnan()` guards (GP trees produce NaN routinely)
- Remaining sub-flags (`-funsafe-math-optimizations`, `-fassociative-math`, `-freciprocal-math`)
  interact with NEON FMA to distort the fitness landscape vs x86 AVX

Current perf flags: `-O3 -g -march=native -mtune=native -funroll-loops -flto`

**Throughput regression:** ~200 sims/sec on short paths (longSequential) vs 5000 with
`-ffast-math`. CPU utilization drops to ~5% per worker — not slow computation but
idle/waiting. Longer paths (progressiveDistance) peg CPUs normally. Root cause TBD.

### Temporal History Wipe (crrcsim — fixed)

`aircraftState` was aggregate-reconstructed every eval tick, destroying the temporal
history ring buffer. Fixed with in-place setters + `clearHistory()` on path reset.
Committed to crrcsim main as `f65c8ad`.

## References

- BACKLOG.md: "Path Interpolation for Smooth Target Tracking"
- BACKLOG.md: "Error Cone for Future Path Points" (future work)
- BACKLOG.md: "Target Pose Estimation/Interpolation" (future work)
- autoc/tests/gp_evaluator_tests.cc (EarlyTimestampOvershoot test)
- autoc/specs/UNIFY.md Appendix D — detailed bug fix notes

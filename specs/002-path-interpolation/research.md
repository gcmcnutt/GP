# Research: Path Interpolation & Evaluator Improvements

**Feature**: 002-path-interpolation
**Date**: 2026-03-06

## Research Items

### 1. Current `getPathIndex()` Implementation

**Location**: `autoc/aircraft_state.h:179-209`

**Current behavior**:
```cpp
inline int getPathIndex(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {
    // Clamp steps to -5..+5
    int steps = CLAMP_DEF((int)arg, -5, 5);

    // Calculate goal time
    const gp_scalar currentTimeMsec = pathProvider.getPath(currentStep).simTimeMsec;
    const gp_scalar timeGoalMsec = currentTimeMsec + steps * SIM_TIME_STEP_MSEC;

    // Linear search forward/backward to find first waypoint at-or-after goal time
    if (steps > 0) {
        while (currentStep < pathProvider.getPathSize() - 1 &&
               pathProvider.getPath(currentStep).simTimeMsec + kEps < timeGoalMsec) {
            currentStep++;
        }
    }
    // ... backward search similar
    return currentStep;
}
```

**Bug identified**: Returns discrete waypoint index, not interpolated position. When timing jitter causes goal time to be slightly early (e.g., 98ms instead of 100ms), the algorithm overshoots to next waypoint (200ms) causing sensor value jumps.

**Decision**: Replace with `getInterpolatedTargetPosition()` that returns a position vector (linear interpolation between bracketing waypoints).

**Rationale**: Interpolation eliminates discrete jumps, making sensors robust to real-time timing jitter.

**Alternatives considered**:
- Improve epsilon handling in current function - rejected because still returns discrete index
- Add position caching at discrete indices - rejected because doesn't solve jitter issue

### 2. FITNESS_SIMPLIFY Implementation Status

**Decision**: Already implemented in `autoc/autoc.cc`

**Evidence**:
- `DISTANCE_POWER = 1.2` defined in `autoc.h:19`
- `attitude_scale` computed per-path at lines 1165-1170, 1736-1741
- Two-objective formula at lines 1193, 1764: `distance_sum += pow(distance, DISTANCE_POWER); attitude_sum += attitude_delta`
- Final fitness uses `distance_sum + attitude_sum * attitude_scale`

**Verification needed**: Unit test for `computeAttitudeScale()` behavior with edge cases (straight-line paths where `path_turn_rad → 0`).

### 3. RAMP_LANDSCAPE Implementation Status

**Decision**: Already implemented in `autoc/autoc.cc`

**Evidence**:
- `gVariationRampStep` global at line 108
- `computeVariationScale()` function at lines 116-125
- Applied to all variation sigmas in `populateVariationOffsets()` at lines 282-290
- Applied to rabbit speed variation at lines 309-311
- Config loaded from `VariationRampStep` in `config_manager.cc:119`
- Logging at generation step boundaries at lines 2007-2014

**Verification needed**: Unit test for `computeVariationScale()` at boundary conditions (gen=0, mid-training, end of training).

### 4. Binary Search Algorithm Design

**Decision**: Use `std::lower_bound` for O(log n) lookup instead of linear search.

**Rationale**: Paths can have 100-1000 waypoints. Binary search is more efficient and clearer.

**Implementation approach**:
```cpp
// Binary search for lower bracket: path[i].simTimeMsec <= goalTime
auto it = std::lower_bound(path.begin(), path.end(), goalTimeMsec,
    [](const Path& p, gp_scalar t) { return p.simTimeMsec < t; });

// Get bracketing indices
size_t i1 = std::distance(path.begin(), it);
size_t i0 = (i1 > 0) ? i1 - 1 : 0;

// Linear interpolation
gp_scalar frac = (goalTimeMsec - path[i0].simTimeMsec) /
                 (path[i1].simTimeMsec - path[i0].simTimeMsec);
return path[i0].start + frac * (path[i1].start - path[i0].start);
```

### 5. Step Range Expansion

**Decision**: Expand step clamp from ±5 to ±10 (±1 second at 100ms/step)

**Rationale**: Spec clarification #1 specified ±1 second as the boundary. Current code clamps to ±5 (500ms).

**Change required**: Update `CLAMP_DEF((int)arg, -5, 5)` to `CLAMP_DEF((int)arg, -10, 10)` or use a named constant.

## Unknowns Resolved

| Unknown | Resolution |
|---------|------------|
| How is getPathIndex implemented? | Linear search with overshoot bug - replace with interpolation |
| Is FITNESS_SIMPLIFY implemented? | Yes, verified in autoc.cc |
| Is RAMP_LANDSCAPE implemented? | Yes, verified in autoc.cc |
| What algorithm for interpolation? | Binary search + linear lerp |
| Step range boundary? | ±10 steps (±1 second) per spec clarification |

## Dependencies

| Dependency | Status | Notes |
|------------|--------|-------|
| Eigen3 | Available | Used for `gp_vec3` position vectors |
| GoogleTest | Available | Existing test framework |
| std::lower_bound | Available | C++17 STL algorithm |

## No Outstanding Research

All NEEDS CLARIFICATION items resolved. Ready for Phase 1 design.

# Data Model: Path Interpolation

**Feature**: 002-path-interpolation
**Date**: 2026-03-06

## Entities

### Path (existing)

Represents a single waypoint on the target trajectory.

| Field | Type | Description |
|-------|------|-------------|
| `start` | `gp_vec3` | Position in NED coordinates (meters) |
| `simTimeMsec` | `gp_scalar` | Simulation time at this waypoint (milliseconds) |
| `distanceFromStart` | `gp_scalar` | Cumulative distance from path origin (meters) |
| `radiansFromStart` | `gp_scalar` | Cumulative turn angle from path origin (radians) |

**Location**: `autoc/aircraft_state.h`

### PathProvider (existing)

Container for path waypoints with index tracking.

| Method | Returns | Description |
|--------|---------|-------------|
| `getPath(int idx)` | `Path&` | Get waypoint at index |
| `getPathSize()` | `int` | Number of waypoints |
| `getCurrentIndex()` | `int` | Current rabbit index |
| `getMaxTimeMsec()` | `gp_scalar` | **NEW**: Last waypoint time |

**Location**: `autoc/aircraft_state.h`

### AircraftState (existing)

Current aircraft state for evaluation.

| Field | Type | Description |
|-------|------|-------------|
| `simTimeMsec` | `gp_scalar` | Current simulation time (milliseconds) |
| `position` | `gp_vec3` | Aircraft position in NED (meters) |
| ... | | (other fields unchanged) |

**Location**: `autoc/aircraft_state.h`

## New Functions

### getInterpolatedTargetPosition()

**Signature**:
```cpp
gp_vec3 getInterpolatedTargetPosition(
    PathProvider& pathProvider,
    gp_scalar currentTimeMsec,
    gp_scalar offsetSteps
);
```

**Parameters**:
- `pathProvider`: Path waypoints container
- `currentTimeMsec`: Current simulation time (milliseconds)
- `offsetSteps`: Lookahead/lookbehind steps (positive = future, negative = past)

**Returns**: Interpolated position in NED coordinates (gp_vec3)

**Algorithm**:
1. Clamp `offsetSteps` to ±10 (±1 second)
2. Calculate goal time: `goalTime = currentTime + offsetSteps * 100ms`
3. Binary search for bracketing waypoints: `path[i].time <= goalTime < path[i+1].time`
4. Linear interpolate: `pos = lerp(path[i].start, path[i+1].start, frac)`
5. Return interpolated position

**Edge cases**:
- `goalTime < path[0].time`: Return `path[0].start`
- `goalTime >= path[N-1].time`: Return `path[N-1].start`
- Single waypoint path: Return that waypoint
- NaN offsetSteps: Return position at current rabbit index

**Location**: `autoc/gp_evaluator_portable.h/.cc`

### Removed Functions

| Function | Reason |
|----------|--------|
| `getPathIndex()` | Replaced by `getInterpolatedTargetPosition()` |

## Updated Functions

### executeGetDPhi()

**Before**:
```cpp
int idx = getPathIndex(pathProvider, aircraftState, arg);
gp_vec3 targetPos = pathProvider.getPath(idx).start;
```

**After**:
```cpp
gp_vec3 targetPos = getInterpolatedTargetPosition(
    pathProvider, aircraftState.simTimeMsec, arg);
```

Same update for: `executeGetDTheta()`, `executeGetDTarget()`

## Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `SIM_TIME_STEP_MSEC` | 100 | Milliseconds per step (existing) |
| `MAX_OFFSET_STEPS` | 10 | Maximum lookahead/lookbehind steps |
| `MAX_OFFSET_MSEC` | 1000 | Maximum offset in milliseconds |

## Validation Rules

1. `offsetSteps` MUST be clamped to `[-MAX_OFFSET_STEPS, MAX_OFFSET_STEPS]`
2. NaN input MUST return current rabbit position (safe fallback)
3. Empty path MUST be handled gracefully (return zero vector or assert)
4. Time interpolation fraction MUST be clamped to `[0, 1]` for numerical stability

## State Transitions

Not applicable - pure functional transformation with no state changes.

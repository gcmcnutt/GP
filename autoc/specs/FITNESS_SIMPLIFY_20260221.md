# Fitness Function Simplification

**Date:** 2026-02-21
**Status:** Proposed
**Problem:** Current 4-objective fitness allows degenerate solutions (fly low-and-fast matches energy but wrong altitude)
**Context:** Basic flight controller layer - "dumb controller" that tracks path smoothly. Strategy (below bad, ahead bad) belongs in higher layers (see layered controller specs).

## Evolution (Not Cycling)

| Phase | Metrics | Lesson |
|-------|---------|--------|
| Original | 4 various | Too complex, unclear gradients |
| Jan 2026 | 3 (waypoint, direction, energy) | Removed throttle, added energy |
| Recent | +attitude_delta (4th) | Needed for smoothness |
| **Now** | 2 (distance + attitude) | Direction redundant, energy has loophole |

Each step learned something concrete - this is genuine simplification.

## Current State

Four objectives with different units and power-law scaling:

| Objective | Units | Power | Issue |
|-----------|-------|-------|-------|
| Waypoint Distance | meters | 1.2-1.8 | Core signal, keep |
| Movement Direction | scaled 0-40 | 1.3 | Redundant - bad direction → distance grows |
| Energy Deviation | J/kg | 1.0-1.5 | **Loophole**: low+fast = same energy |
| Attitude Delta | radians | 1.0 | Needed, but wrong scale |

## Proposed: Two Objectives

1. **Distance to rabbit** (meters) - primary tracking signal
2. **Attitude delta** (radians → scaled to meters) - smoothness constraint

### Why This Works

- **Distance alone handles direction**: If you're not aligned with path, rabbit moves away, distance grows
- **Distance alone handles altitude**: Below target = more meters away
- **Attitude delta prevents spiraling**: Constant radius spiral = constant distance, but roll changes explode
- **No energy loophole**: Can't game distance by trading speed for altitude

## Normalization Problem

Raw units don't compare:
- Distance: ~5m average error × 100 steps = 500 meter-steps
- Attitude: ~44 rad optimal (7 revolutions) + some excess = small number

Need to scale attitude to "effective meters" so objectives are comparable.

## Solution: Path-Relative Scaling

Use path geometry to normalize:

```
path_distance = path.back().distanceFromStart   // total meters
path_turn_rad = path.back().radiansFromStart    // total radians (turnmeter)
attitude_scale = path_distance / path_turn_rad  // meters per radian
```

For a 500m path with 44 rad of turns: `attitude_scale ≈ 11.4` meters/radian.

This means: "1 radian of excess turning costs same as 11.4m of distance error"

### Fitness Formula

```c
// Per-path: compute scale from THIS path's geometry (not global)
// Important: training may have multiple paths per scenario
gp_scalar path_distance = path.back().distanceFromStart;
gp_scalar path_turn_rad = path.back().radiansFromStart;
gp_scalar attitude_scale = path_distance / std::max(path_turn_rad, 1.0f);

// Per-step fitness (within this path)
gp_scalar distance = (aircraftPosition - rabbitPosition).norm();
gp_scalar attitude_change = fabs(roll - prev_roll) + fabs(pitch - prev_pitch);

// Simple nonlinear penalty - small excursions still matter
constexpr gp_scalar DISTANCE_POWER = 1.2;

gp_fitness step_fitness = pow(distance, DISTANCE_POWER) +
                          attitude_change * attitude_scale;

step_fitness_sum += step_fitness;  // drop STEP_TIME_WEIGHT (constant multiplier)
```

**Multi-path scenarios**: Each path computes its own `attitude_scale` from its own geometry. This ensures a tight 100m circuit and a long 500m path both get appropriately scaled attitude penalties.

## Properties

1. **Self-normalizing**: No manually tuned constants - scale derived from path
2. **Physically meaningful**: Both terms in effective meters
3. **Gradient-preserving**: Power-law on distance keeps selection pressure
4. **Degenerate-resistant**: No loopholes for low-fast or spiral solutions

## Alternative: Accumulate Then Scale (Per-Path)

Could also accumulate raw values and scale once at path end:

```c
// During this path's flight
constexpr gp_scalar DISTANCE_POWER = 1.2;
distance_sum += pow(distance, DISTANCE_POWER);
attitude_sum += attitude_change;  // raw radians (roll + pitch)

// At this path's end
gp_scalar attitude_scale = path_distance / std::max(path_turn_rad, 1.0f);
localFitness = distance_sum + attitude_sum * attitude_scale;

// Aggregate across paths in scenario
stdFitness += localFitness;
```

Same result, slightly cleaner accumulation. Scale computed fresh for each path.

## Edge Cases

- **path_turn_rad = 0** (straight line): Use fallback scale, e.g., `path_distance / (2*PI)`
- **Very short paths**: Scale still works, just smaller numbers
- **Crash mid-flight**: Scale from partial path? Or use full path scale + completion penalty

## Migration

**Remove:**
- `MOVEMENT_DIRECTION_WEIGHT`, `DIRECTION_SCALE` - direction alignment redundant
- `ALTITUDE_LOW_POWER`, `ALTITUDE_HIGH_POWER` - energy metric removed entirely
- `STEP_TIME_WEIGHT` - physics cruft, dimensionless now
- `ATTITUDE_DELTA_SCALE` - replaced by path-relative scaling
- `WAYPOINT_AHEAD_POWER`, `WAYPOINT_BEHIND_POWER` - asymmetry belongs in strategy layer
- All energy calculation code

**Keep:**
- `CRASH_COMPLETION_WEIGHT` - orthogonal concern (soft lexicographic)

**Add:**
- `DISTANCE_POWER = 1.2` - single nonlinear exponent (small excursions still expensive)

## Design Decisions

1. **No altitude asymmetry** - "below is bad" belongs in higher strategy layer, not dumb controller
2. **Roll + pitch equal** - let evolution find optimal weighting; both contribute to energy/turn rate
3. **Single distance power** - keep it simple; 1.2 provides gradient without over-penalizing

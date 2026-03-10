# [TRANSFERRED to specs/010-safety-layer] LAYERED_CONTROLLER: Evolving Beyond Basic Tracking

## The Current Problem

The GP has one job: output PITCH, ROLL, THROTTLE commands each timestep. It receives:
- Raw state (position, velocity, orientation)
- Navigation helpers (GETDPHI, GETDTHETA, GETDTARGET)

This is enough for basic "go there" but the GP has no awareness of:
- Am I about to crash?
- Am I overshooting the path?
- Am I in an intercept phase vs tracking phase?
- What's my energy state?
- What's coming up on the path?

The result: GP evolves reactive controllers that work "okay" but can't handle edge cases or plan ahead.

## Layered Architecture Concept

```
┌─────────────────────────────────────────────────────────┐
│ STRATEGY LAYER (optional, evolved separately?)          │
│ - Phase detection (intercept → track → recover)         │
│ - Energy management                                      │
│ - Path preview / anticipation                           │
├─────────────────────────────────────────────────────────┤
│ SAFETY LAYER (constraints, always active)               │
│ - Envelope protection (don't exceed limits)             │
│ - Geofence (stay in bounds)                             │
│ - Stall prevention                                       │
│ - "Don't get ahead of rabbit" limiting                  │
├─────────────────────────────────────────────────────────┤
│ TRACKING LAYER (core GP evolution target)               │
│ - Basic attitude control                                │
│ - Path following                                         │
│ - Speed regulation                                       │
└─────────────────────────────────────────────────────────┘
         ↓ PITCH, ROLL, THROTTLE commands ↓
```

## Option A: Single GP with Dynamic State Flags

Instead of separate layers, give the GP more situational awareness via new terminals.

### New State Flags (Boolean/Normalized)

| Terminal | Description | Range |
|----------|-------------|-------|
| `ISINTERCEPTING` | Target far ahead, need to catch up | 0/1 |
| `ISTRACKING` | On path, maintaining | 0/1 |
| `ISRECOVERING` | Off path, returning | 0/1 |
| `ISOVERSHOOTING` | Ahead of rabbit | 0/1 |
| `NEARSTALL` | Airspeed dangerously low | 0/1 |
| `NEARGROUND` | Altitude dangerously low | 0/1 |
| `NEARBOUNDS` | Approaching geofence | 0/1 |
| `HIGHENERGY` | Fast/high, can trade for maneuver | 0/1 |
| `LOWENERGY` | Slow/low, need to conserve | 0/1 |

### New Continuous State

| Terminal | Description | Range |
|----------|-------------|-------|
| `GETTRACKERROR` | Cross-track error (lateral distance to path) | [-1, 1] normalized |
| `GETPATHERROR` | Along-track error (ahead/behind rabbit) | [-1, 1] normalized |
| `GETENERGY` | Specific energy (altitude + speed²) | [-1, 1] normalized |
| `GETPATHCURVE` | Upcoming path curvature | [-1, 1] |
| `GETPATHCLIMB` | Upcoming path climb rate | [-1, 1] |

### Pros
- Single GP, simpler architecture
- GP can learn to use flags or ignore them
- Backward compatible (existing GPs still work)

### Cons
- Even more terminals (search space grows)
- GP may never learn to use the flags effectively
- No guarantee safety constraints are respected

## Option B: Hardcoded Safety + Evolved Tracking

Safety layer is NOT evolved - it's hand-coded constraints that clamp GP outputs.

```cpp
// GP produces raw commands
float gp_pitch = evaluateGP();
float gp_roll = evaluateGP();
float gp_throttle = evaluateGP();

// Safety layer overrides
if (altitude < MIN_ALT) {
    pitch = max(pitch, RECOVERY_PITCH);  // Force climb
}
if (airspeed < STALL_SPEED) {
    pitch = min(pitch, DIVE_PITCH);      // Force dive for speed
    throttle = MAX_THROTTLE;
}
if (nearBoundary()) {
    // Turn toward center
    roll = computeTurnAwayFromBoundary();
}

// Final output
setPitch(clamp(pitch, -1, 1));
setRoll(clamp(roll, -1, 1));
setThrottle(clamp(throttle, 0, 1));
```

### Pros
- Safety guaranteed (no evolution needed)
- GP can focus on tracking only
- Smaller search space for GP

### Cons
- Hardcoded safety may fight with GP intent
- Less elegant
- Safety layer needs manual tuning

## Option C: Hierarchical GP (Two Populations)

Evolve two GPs:
1. **Tracking GP** - basic path following (current)
2. **Supervisor GP** - modifies tracking GP outputs based on situation

```
Supervisor inputs: situation flags + Tracking GP outputs
Supervisor outputs: modified PITCH, ROLL, THROTTLE

// Or alternatively:
Supervisor outputs: gain adjustments, mode switches, limits
```

### Pros
- Clean separation of concerns
- Each GP has simpler job
- Can evolve tracking first, then supervisor

### Cons
- Two populations to manage
- Co-evolution is tricky (moving target)
- More complex infrastructure

## Option D: Phase-Based Controller Selection

Different GPs for different phases, with a simple state machine selecting which is active.

```
State Machine:
  INTERCEPT → (when close enough) → TRACK → (if lost) → RECOVER → INTERCEPT
                                         ↓
                                    (if danger)
                                         ↓
                                      ESCAPE
```

Each phase has its own evolved GP or even hand-coded controller:
- INTERCEPT: Aggressive pursuit, high speed
- TRACK: Smooth following, minimize error
- RECOVER: Return to path
- ESCAPE: Emergency maneuver (stall/ground/boundary)

### Pros
- Each GP has very focused job
- Phase transitions are explicit
- Easy to debug ("why did it crash?" → check which phase, which GP)

### Cons
- Multiple GPs to evolve
- State machine transitions need tuning
- Discontinuities at phase boundaries

## Recommended Approach: Incremental Enhancement

### Phase 1: Add Situational Awareness (Low Risk)

Add these terminals to current GP:

```cpp
// Flight phase indicators
ns.putNode(*new GPNode(GETTRACKERROR, "GETTRACKERROR"));  // Cross-track error
ns.putNode(*new GPNode(GETPATHERROR, "GETPATHERROR"));    // Along-track error (+ = ahead)

// Safety indicators
ns.putNode(*new GPNode(GETALTMARGIN, "GETALTMARGIN"));    // Height above min safe alt
ns.putNode(*new GPNode(GETSPEEDMARGIN, "GETSPEEDMARGIN"));// Speed above stall

// Energy state
ns.putNode(*new GPNode(GETENERGY, "GETENERGY"));          // Normalized specific energy
```

Implementation is straightforward - these are just computed values from existing state.

### Phase 2: Hardcoded Safety Envelope (Medium Risk)

Add post-GP safety clamping:

```cpp
// After GP evaluation, before sending to sim:
if (getAltitude() < CRITICAL_ALT) {
    pitch = std::max(pitch, MIN_PITCH_FOR_CLIMB);
}
if (getAirspeed() < CRITICAL_SPEED) {
    throttle = MAX_THROTTLE;
    pitch = std::min(pitch, MAX_PITCH_FOR_SPEED);
}
```

This catches the GP's worst mistakes without requiring it to learn safety.

### Phase 3: Remove Useless Terminals (Node Diet)

Based on analysis of successful GPs, remove terminals that never appear in good solutions:

Candidates for removal:
- GETVELX, GETVELY, GETVELZ (use GETVEL instead)
- GETALPHA, GETBETA (aerodynamic angles rarely useful at this level)
- GETPITCH, GETROLL, GETTHROTTLE (readback of own commands)

### Phase 4: Evaluate Higher-Level Nodes (Experimental)

If Phase 1-3 don't improve success rate, try collapsing common patterns:

```cpp
// TRACKPITCH = SETPITCH(GETDTHETA(reasonable_lookahead))
ns.putNode(*new GPNode(TRACKPITCH, "TRACKPITCH", 1));  // arg = aggressiveness

// TRACKROLL = SETROLL(GETDPHI(reasonable_lookahead))
ns.putNode(*new GPNode(TRACKROLL, "TRACKROLL", 1));    // arg = aggressiveness
```

These encode the "known good" structure while allowing GP to tune the parameters.

## Computing the New Terminals

### GETTRACKERROR (Cross-Track Error)

```cpp
gp_scalar executeGetTrackError(PathProvider& pathProvider, AircraftState& state) {
    // Find closest point on path
    int currentIdx = state.getThisPathIndex();
    gp_vec3 pathPoint = pathProvider.getPath(currentIdx).start;
    gp_vec3 pathTangent = pathProvider.getPath(currentIdx).direction;

    // Vector from path point to aircraft
    gp_vec3 offset = state.getPosition() - pathPoint;

    // Remove along-track component
    gp_vec3 crossTrack = offset - pathTangent * offset.dot(pathTangent);

    // Signed lateral error (positive = right of path)
    gp_scalar lateralError = crossTrack.dot(state.getRightVector());

    // Normalize to [-1, 1] with 10m = full scale
    return CLAMP_DEF(lateralError / 10.0f, -1.0f, 1.0f);
}
```

### GETPATHERROR (Along-Track Error)

```cpp
gp_scalar executeGetPathError(PathProvider& pathProvider, AircraftState& state) {
    // How far ahead/behind the rabbit are we?
    int currentIdx = state.getThisPathIndex();
    gp_vec3 rabbitPos = pathProvider.getPath(currentIdx).start;
    gp_vec3 pathTangent = pathProvider.getPath(currentIdx).direction;

    gp_vec3 toRabbit = rabbitPos - state.getPosition();
    gp_scalar alongTrack = toRabbit.dot(pathTangent);

    // Positive = behind rabbit (need to catch up)
    // Negative = ahead of rabbit (slow down)
    // Normalize with 20m = full scale
    return CLAMP_DEF(alongTrack / 20.0f, -1.0f, 1.0f);
}
```

### GETENERGY (Specific Energy)

```cpp
gp_scalar executeGetEnergy(AircraftState& state) {
    // Specific energy: E = h + v²/2g
    gp_scalar altitude = -state.getPosition().z();  // NED: -z = up
    gp_scalar speed = state.getRelVel();

    gp_scalar kineticTerm = (speed * speed) / (2.0f * 9.81f);
    gp_scalar totalEnergy = altitude + kineticTerm;

    // Normalize relative to nominal flight (50m alt, 15 m/s)
    gp_scalar nominalEnergy = 50.0f + (15.0f * 15.0f) / (2.0f * 9.81f);

    return CLAMP_DEF((totalEnergy - nominalEnergy) / nominalEnergy, -1.0f, 1.0f);
}
```

## Open Questions

1. **Should safety override be visible to GP?** If GP knows its outputs will be clamped, it might learn to work with that. Or it might learn to depend on it.

2. **How to evaluate strategy?** Fitness currently rewards path completion. Strategy (energy management, anticipation) only helps indirectly. May need explicit fitness terms.

3. **Continuous vs discrete modes?** Phase flags (ISINTERCEPTING) are discrete. Cross-track error is continuous. GP may handle continuous better (no sharp transitions).

4. **Lookahead as a first-class concept?** Many issues trace back to the lookahead parameter for GETDPHI/GETDTHETA. Maybe:
   - GETDPHI_NEAR, GETDPHI_FAR (fixed lookaheads)
   - Or compute optimal lookahead based on speed/curvature

5. **Fitness function alignment?** If we add energy awareness but don't penalize energy waste, GP won't use it. Fitness may need terms for:
   - Smoothness (penalize control oscillation)
   - Efficiency (penalize throttle extremes)
   - Anticipation (reward early response to upcoming curves)

---

## Observations from 003-variations-redux (2026-03-09)

### Lookahead is the Binding Constraint on Tracking Distance

Analysis of a training run with tightened fitness params (DISTANCE_NORM=2.0, DISTANCE_POWER=2.0) showed:

- **Median tracking distance: ~10.5m**, down from ~15-20m baseline
- **Turn overshoots peak at ~16m** with dphi climbing from 22° to 46° through tight turns
- **Pitch saturated at 1.0** during worst overshoots — the aircraft is physically trying as hard as it can
- **Straight segment tracking is ~8-9m** — the fitness tuning helped here

The GP's lookahead window via `GETDPHI(steps)` / `GETDTHETA(steps)` gives roughly 10-30m of path preview at ~16 m/s cruise speed. When a tight turn begins, the aircraft is already committed by the time it sees the curvature ramp up. It can't pre-position or anticipate.

**Key insight**: The ~10m median is likely near the **physical ceiling** for reactive-only control with this craft, path geometry, and lookahead window. Further fitness tuning will yield diminishing returns. The turn overshoots are a planning/anticipation problem, not a control responsiveness problem.

### What a Strategy Layer Would Need to Address

1. **Turn anticipation**: "A sharp turn is coming in 50m — start bleeding speed and banking early." This requires lookahead beyond the GP's current ~10-30m reactive window.
2. **Path feasibility**: "This turn radius is tighter than my min turn radius at current speed — I need to slow down NOW, not when I'm in the turn."
3. **Arena awareness**: "The path curves back toward center after this turn — don't panic-correct, the geometry will bring me back."

These are all examples of reasoning that requires a longer planning horizon than the current GP terminal set provides. The `GETPATHCURVE` terminal proposed in Option A above would be a lightweight first step.

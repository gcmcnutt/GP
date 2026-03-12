# Research: Distance Temporal Sensor Nodes

**Feature**: 012-distance-temporal-nodes
**Date**: 2026-03-12

## No NEEDS CLARIFICATION Items

All technical decisions are resolved by the spec and the existing GETDPHI_PREV/RATE implementation pattern.

## Pattern Analysis: Existing Temporal Nodes

### Decision: Follow GETDPHI_PREV/RATE pattern exactly

**Rationale**: The existing temporal nodes (implemented in 003-variations-redux, documented in ZZZ-TEMPORAL_STATE.md) provide a proven, tested pattern across all 15 touchpoint files. Replicating this pattern for distance minimizes risk and ensures consistency.

**Alternatives considered**:
- Custom ring buffer per node: rejected — would duplicate infrastructure and diverge from the shared-index pattern
- Separate time buffer for distance: rejected — the existing timeHistory_[] is shared across all buffers and indexed by the same historyIndex_

### Key Pattern Details (from codebase analysis)

**Ring buffer structure** (aircraft_state.h):
- `gp_scalar dPhiHistory_[HISTORY_SIZE]` — one float array per signal
- `unsigned long timeHistory_[HISTORY_SIZE]` — shared timestamp array
- `int historyIndex_` — shared write pointer (wraps at HISTORY_SIZE)
- `int historyCount_` — shared valid sample count
- HISTORY_SIZE = 10 (1 second at 10Hz)

**Recording** (minisim.cc line 240):
- Called before GP evaluation each tick
- `aircraftState.recordErrorHistory(dPhi, dTheta, duration_msec)`
- New: extend signature to include distance, OR add separate recordDistanceHistory() call

**Rate computation** (gp_evaluator_portable.cc lines 410-437):
- `dt = (timestamp[0] - timestamp[1]) / 1000.0f` (ms → seconds)
- Default dt = 0.1f when timestamps equal or dt < 0.001f
- Clamp result to [-10, 10]
- Return 0.0f if historyCount < 2

**Accessor** (aircraft_state.h lines 275-287):
- `getHistoricalDPhi(int n)` — lookback n entries in ring buffer
- Clamp n to [0, historyCount_ - 1]
- Index calculation: `(historyIndex_ - 1 - n + HISTORY_SIZE) % HISTORY_SIZE`

### Decision: GETDIST is nullary (no lookahead argument)

**Rationale**: Unlike GETDPHI(n)/GETDTHETA(n) which take a path lookahead step, distance is a scalar — there's no directional component that would benefit from lookahead. The GP can combine GETDIST with other path-relative sensors if needed.

**Alternatives considered**:
- GETDIST(n) with lookahead: rejected — distance to a future path point is less useful than distance to current rabbit position. The GP already has GETDPHI(n)/GETDTHETA(n) for directional lookahead.

### Decision: Buffer raw distance, not GETDTARGET output

**Rationale**: GETDTARGET returns `CLAMP((distance - 10) / relVel, -1, 1)` — a composite signal mixing distance and speed with a hard-coded offset. Buffering raw meters gives the GP cleaner primitives. The GP already has GETVEL for speed, so it can compose distance/speed relationships itself.

**Alternatives considered**:
- Buffer GETDTARGET output: rejected — muddied signal, speed changes affect derivative even when distance is constant
- Buffer both raw and GETDTARGET: rejected — redundancy dilutes selection pressure

### Decision: Deprecate GETDTARGET from TrainingNodes

**Rationale**: GETDTARGET is a "pre-computed shortcut" that bakes in assumptions (10m offset, speed normalization) the GP should discover on its own from primitives. Keeping it alongside GETDIST would create redundancy. Opcode remains for backward compatibility with existing evolved GPs.

## Units Consistency Verification

| Node | Returns | Units | Rate units |
|------|---------|-------|------------|
| GETDPHI(n) | Roll bearing to target | radians | GETDPHI_RATE → rad/s |
| GETDTHETA(n) | Pitch angle to target | radians | GETDTHETA_RATE → rad/s |
| GETDIST | Distance to rabbit | meters | GETDIST_RATE → m/s |
| GETVEL | Airspeed | m/s | (no rate node) |
| GETDTARGET(n) | (distance-10)/speed | dimensionless | (no rate — deprecated) |

All rate nodes use the same dt computation: `(timestamp[0] - timestamp[1]) / 1000.0f` with 0.1s default.
All rate nodes clamp to [-10, 10] in their respective units.
All history nodes share the same ring buffer index and timestamp array.

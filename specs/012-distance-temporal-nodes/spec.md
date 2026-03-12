# Feature Specification: Distance Temporal Sensor Nodes

**Feature Branch**: `012-distance-temporal-nodes`
**Created**: 2026-03-12
**Status**: Draft
**Input**: User description: "Add temporal distance sensor nodes (GETDTARGET_PREV, GETDTARGET_RATE) to give GP controllers derivative and historical distance-to-rabbit awareness for PD/PID-like throttle control"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - GP Evolves Distance Derivative Control (Priority: P1)

The GP evolution engine has access to new temporal distance nodes (GETDTARGET_PREV and GETDTARGET_RATE) alongside existing nodes. During evolution, GP trees can use GETDTARGET_RATE to sense whether the aircraft is closing on or falling behind the rabbit, enabling proportional-derivative throttle control that damps the current 10-55m distance oscillation.

**Why this priority**: This is the core value — without the rate signal, the GP can only do proportional distance control, which produces undamped oscillation around the target following distance. The existing GETDPHI_RATE/GETDTHETA_RATE nodes proved that temporal derivatives unlock PD-like control for attitude; this extends the same pattern to distance.

**Independent Test**: Can be fully tested by running a short evolution (10-20 generations) and verifying that GP trees incorporate GETDTARGET_RATE nodes in their throttle control expressions.

**Acceptance Scenarios**:

1. **Given** a GP configuration with GETDTARGET_RATE in TrainingNodes, **When** evolution runs for 20+ generations, **Then** at least some individuals in the population use GETDTARGET_RATE in their trees
2. **Given** a GP tree using GETDTARGET_RATE for throttle, **When** the aircraft is closing on the rabbit (distance decreasing), **Then** GETDTARGET_RATE returns a negative value (closing rate)
3. **Given** a GP tree using GETDTARGET_RATE for throttle, **When** the aircraft is falling behind (distance increasing), **Then** GETDTARGET_RATE returns a positive value (opening rate)

---

### User Story 2 - Distance History Enables Integral-Like Control (Priority: P2)

The GP evolution engine has access to GETDTARGET_PREV(n), which returns the distance-to-rabbit value from n ticks ago. This enables the GP to compute error accumulation patterns (how long has the aircraft been too far/too close), supporting integral-like control strategies that eliminate steady-state distance offset.

**Why this priority**: Historical distance awareness complements the rate signal. While GETDTARGET_RATE addresses oscillation damping, GETDTARGET_PREV enables the GP to detect persistent distance bias and evolve corrective strategies. This mirrors the proven GETDPHI_PREV/GETDTHETA_PREV pattern.

**Independent Test**: Can be tested by verifying GETDTARGET_PREV(n) returns the correct historical distance value for indices 0 through the buffer depth, using unit tests with known distance sequences.

**Acceptance Scenarios**:

1. **Given** the aircraft has been tracking at 20m distance for 5 ticks, **When** GETDTARGET_PREV(0) through GETDTARGET_PREV(4) are queried, **Then** all return approximately 20m distance
2. **Given** the aircraft distance changed from 30m to 15m over 5 ticks, **When** GETDTARGET_PREV(0) and GETDTARGET_PREV(4) are queried, **Then** GETDTARGET_PREV(0) returns ~15m and GETDTARGET_PREV(4) returns ~30m
3. **Given** fewer than n ticks have elapsed since simulation start, **When** GETDTARGET_PREV(n) is queried, **Then** the system returns a safe default (zero or the oldest available value)

---

### User Story 3 - Cross-Platform Consistency (Priority: P2)

The new distance temporal nodes produce identical results across all evaluation paths: GP tree evaluation during training, bytecode interpretation during deployment, and the portable evaluator for embedded targets. A GP controller evolved during training behaves identically when extracted to bytecode and run on a different platform.

**Why this priority**: Cross-platform consistency is a hard requirement for the dual-mode evaluation architecture. Any divergence between training and deployment would make evolved controllers unreliable.

**Independent Test**: Can be tested by evolving a GP that uses GETDTARGET_PREV/GETDTARGET_RATE, extracting to bytecode, and comparing per-step outputs between GP tree and bytecode evaluation modes.

**Acceptance Scenarios**:

1. **Given** a GP tree using GETDTARGET_RATE, **When** evaluated as GP tree and as bytecode on the same scenario, **Then** all control outputs match to floating-point precision
2. **Given** the portable evaluator with GETDTARGET_PREV support, **When** compiled with GP_BUILD and GP_TEST defines, **Then** unit tests confirm identical output to the full evaluator for known input sequences

---

### Edge Cases

- What happens when GETDTARGET_PREV(n) is called with n larger than the history buffer depth (10)? Clamp n to buffer size, return oldest available value.
- What happens when GETDTARGET_PREV is called with a negative or fractional n? Cast to integer, clamp to [0, buffer_depth-1] — same behavior as existing GETDPHI_PREV.
- What happens when GETDTARGET_RATE is computed on the first simulation tick (no prior value)? Return 0.0 — same as existing GETDTHETA_RATE behavior.
- How does the distance history buffer interact with the intercept budget? The buffer stores raw distance values regardless of intercept scaling — the temporal signal reflects actual aircraft-to-rabbit distance, not fitness-weighted distance.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a GETDTARGET_PREV(n) node that returns the distance-to-rabbit value from n simulation ticks ago, using the same circular buffer pattern as GETDPHI_PREV
- **FR-002**: System MUST provide a GETDTARGET_RATE node (zero arguments) that returns the rate of change of distance-to-rabbit in meters per second
- **FR-003**: GETDTARGET_RATE MUST return negative values when distance is decreasing (closing) and positive values when distance is increasing (opening)
- **FR-004**: The distance history buffer MUST store raw distance-to-rabbit values (meters), not the normalized/powered fitness values
- **FR-005**: Both nodes MUST be implemented in all three evaluation paths: GP tree evaluator, bytecode interpreter, and portable evaluator
- **FR-006**: Both nodes MUST be registered as new opcodes in the opcode enumeration with unique bytecode instruction encodings
- **FR-007**: GETDTARGET_RATE MUST be clamped to a reasonable range to prevent extreme values from dominating GP fitness (consistent with GETDPHI_RATE/GETDTHETA_RATE clamping)
- **FR-008**: The distance history buffer MUST use the same depth (10 samples) and ring buffer implementation pattern as the existing phi/theta history buffers
- **FR-009**: Both nodes MUST be configurable in TrainingNodes in autoc.ini (opt-in, not added to all existing configurations by default)
- **FR-010**: Unit tests MUST verify correctness of both nodes for known input sequences, buffer wraparound, edge cases (first tick, buffer overflow), and cross-evaluator consistency

### Key Entities

- **Distance History Buffer**: Circular buffer storing the last 10 distance-to-rabbit measurements (raw meters). One buffer per evaluation context. Updated each simulation tick after distance computation.
- **GETDTARGET_PREV Node**: Unary GP node (1 argument: history index n). Returns buffered distance value at index n. Argument cast to integer, clamped to [0, buffer_depth-1].
- **GETDTARGET_RATE Node**: Nullary GP node (0 arguments). Computes `(distance[0] - distance[1]) / dt` where dt is the simulation tick interval. Clamped to [-10, 10] m/s.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Evolution runs using the new nodes produce GP trees that incorporate GETDTARGET_RATE in throttle control expressions within 20 generations
- **SC-002**: Mean following distance in a 200-generation training run with the new nodes available is closer to the 7.5m target than the current 20.7m baseline
- **SC-003**: Distance standard deviation (currently 7.1m) decreases, indicating reduced oscillation amplitude
- **SC-004**: All existing unit tests continue to pass; new unit tests cover both nodes across all evaluator paths with 100% pass rate
- **SC-005**: Bytecode evaluation of GP trees using the new nodes produces identical per-step control outputs compared to GP tree evaluation (zero divergence in elite re-eval)

## Assumptions

- The simulation tick rate remains ~10Hz (100ms intervals), consistent with the existing temporal node design
- The distance-to-rabbit value is already computed each tick in the evaluation loop before the GP tree is evaluated, so buffering it requires no additional computation
- The existing ring buffer infrastructure (used by GETDPHI_PREV/GETDTHETA_PREV) can be extended with a third buffer for distance without significant memory impact
- GETDTARGET_RATE clamping range of [-10, 10] m/s is appropriate (aircraft closing at 10 m/s relative to rabbit would be an extreme maneuver)

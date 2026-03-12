# Feature Specification: Distance Temporal Sensor Nodes

**Feature Branch**: `012-distance-temporal-nodes`
**Created**: 2026-03-12
**Status**: Draft
**Input**: User description: "Add raw distance sensor with temporal nodes (GETDIST, GETDIST_PREV, GETDIST_RATE) as clean primitives for PD/PID-like throttle control. Deprecate composite GETDTARGET node."

## Context

The GP controller's distance-to-rabbit oscillates between 10-55m because it only has a composite distance signal (GETDTARGET) that bakes in a 10m offset and speed normalization: `CLAMP((distance - 10) / relVel, -1, 1)`. This pre-computed signal:
- Mixes distance and speed into one dimensionless value, preventing the GP from reasoning about them independently
- Hard-codes a 10m offset assumption that may not be optimal
- Has no temporal derivative, so the GP cannot sense closing/opening rate

The existing temporal nodes (GETDPHI_PREV/RATE, GETDTHETA_PREV/RATE) proved that clean primitives in radians with history and derivatives unlock PD-like control for attitude. This feature extends the same pattern to distance using raw meters.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Raw Distance Sensor with Derivative Control (Priority: P1)

The GP evolution engine has access to GETDIST (raw distance to rabbit in meters), GETDIST_RATE (closing/opening rate in m/s), and GETDIST_PREV(n) (distance n ticks ago). Combined with existing GETVEL (airspeed), the GP can compose its own distance-to-throttle mapping and damp oscillation using the rate signal — without being constrained by GETDTARGET's hard-coded assumptions.

**Why this priority**: This is the core value. Raw distance + derivative gives the GP the building blocks to evolve PD throttle control from primitives, just as GETDPHI + GETDPHI_RATE enabled PD roll/pitch control.

**Independent Test**: Run a short evolution (10-20 generations) with GETDIST, GETDIST_RATE, GETDIST_PREV, and GETVEL in TrainingNodes (without GETDTARGET). Verify GP trees use the new nodes.

**Acceptance Scenarios**:

1. **Given** GETDIST in TrainingNodes, **When** the aircraft is 20m from the rabbit, **Then** GETDIST returns 20.0 (meters)
2. **Given** GETDIST_RATE available, **When** the aircraft is closing on the rabbit at 3 m/s, **Then** GETDIST_RATE returns -3.0 (negative = closing)
3. **Given** GETDIST_RATE available, **When** the aircraft is falling behind at 5 m/s, **Then** GETDIST_RATE returns +5.0 (positive = opening)
4. **Given** GETDIST_PREV available, **When** queried with n=3, **Then** returns the raw distance from 3 ticks (~300ms) ago

---

### User Story 2 - GETDTARGET Deprecation (Priority: P1)

GETDTARGET is removed from the active TrainingNodes configuration with a comment explaining why: it's a composite signal that prevents the GP from reasoning about distance and speed independently. The node remains functional in the codebase (for backward compatibility with existing evolved GPs) but is not offered to new evolution runs.

**Why this priority**: Keeping GETDTARGET in the training set alongside GETDIST would create redundancy and dilute selection pressure. Removing it focuses GP search on the cleaner primitives.

**Independent Test**: Verify GETDTARGET is commented out in autoc.ini TrainingNodes with explanatory comment. Verify existing bytecode files using GETDTARGET still evaluate correctly.

**Acceptance Scenarios**:

1. **Given** updated autoc.ini, **When** TrainingNodes is parsed, **Then** GETDTARGET is not in the active node set
2. **Given** an existing GP tree or bytecode using GETDTARGET, **When** evaluated, **Then** it produces identical results (backward compatibility preserved)

---

### User Story 3 - Cross-Platform Consistency (Priority: P2)

All new nodes produce identical results across all evaluation paths: GP tree evaluator, bytecode interpreter, portable evaluator, and bytecode-to-C++ code generator. History buffers reset identically across all paths when starting a new simulation path.

**Why this priority**: Cross-platform consistency is a hard requirement. The bytecode2cpp code generator and gpextractor must also handle the new opcodes.

**Independent Test**: Evolve a GP using GETDIST_RATE, extract to bytecode, compare per-step outputs between GP tree and bytecode modes.

**Acceptance Scenarios**:

1. **Given** a GP tree using GETDIST_RATE, **When** evaluated as GP tree and as bytecode, **Then** all outputs match to floating-point precision
2. **Given** the bytecode2cpp generator, **When** processing bytecode with GETDIST_PREV/GETDIST_RATE opcodes, **Then** it generates correct C++ code with proper stack depth analysis
3. **Given** a new simulation path starts, **When** the distance history buffer resets, **Then** GETDIST_PREV returns 0.0 and GETDIST_RATE returns 0.0 (same reset behavior as GETDPHI_PREV/GETDPHI_RATE)

---

### Edge Cases

- GETDIST_PREV(n) with n > buffer depth (10): clamp to buffer size, return oldest available value
- GETDIST_PREV with negative or fractional n: cast to integer, clamp to [0, depth-1] — same as GETDPHI_PREV
- GETDIST_RATE on first tick (no prior value): return 0.0 — same as GETDPHI_RATE
- GETDIST when aircraft has crashed or is at origin: returns actual Euclidean distance (may be 0 or very large)
- Distance history buffer reset: cleared via clearHistory() at start of each new path, same as dPhi/dTheta buffers

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide GETDIST, a nullary node returning raw Euclidean distance from aircraft to rabbit's current interpolated position in meters (same distance value used in the fitness calculation)
- **FR-002**: System MUST provide GETDIST_PREV(n), a unary node returning the buffered GETDIST value from n simulation ticks ago
- **FR-003**: System MUST provide GETDIST_RATE, a nullary node returning rate of change of distance in meters per second, computed as `(dist[0] - dist[1]) / dt` using actual timestamp deltas from the history buffer (same time computation as GETDPHI_RATE/GETDTHETA_RATE)
- **FR-004**: GETDIST_RATE MUST use the same dt computation as GETDPHI_RATE: actual timestamp delta from ring buffer divided by 1000 to convert ms to seconds, with 0.1s default when dt < 0.001s
- **FR-005**: GETDIST_RATE MUST be clamped to [-10, 10] m/s, consistent with GETDPHI_RATE/GETDTHETA_RATE clamping range
- **FR-006**: The distance history buffer MUST use the same ring buffer depth (10 samples), indexing, and reset behavior as the existing dPhi/dTheta history buffers in AircraftState
- **FR-007**: All three nodes MUST be implemented in all evaluation paths: GP tree evaluator (autoc-eval.cc), bytecode interpreter (gp_bytecode.cc), portable evaluator (gp_evaluator_portable.cc), bytecode-to-C++ generator (bytecode2cpp.cc), and GP-to-bytecode extractor (gpextractor.cc)
- **FR-008**: All three nodes MUST be registered as new opcodes in the enum (autoc.h AND gp_evaluator_portable.h), appended before the `_END` marker to preserve backward compatibility
- **FR-009**: All three nodes MUST be added to the allNodes[] registration table in autoc-eval.cc with correct name strings and argument counts
- **FR-010**: GETDTARGET MUST be removed from TrainingNodes in autoc.ini and autoc-eval.ini with an explanatory comment noting it as deprecated in favor of the cleaner GETDIST primitives. The opcode and evaluation code MUST remain for backward compatibility.
- **FR-011**: The node set reference comment block in autoc.ini and autoc-eval.ini MUST be updated to document the new nodes and the GETDTARGET deprecation
- **FR-012**: CLAUDE.md GP Operators section MUST be updated with the new node descriptions and GETDTARGET deprecation note
- **FR-013**: Unit tests MUST verify all three nodes for: known input sequences, buffer wraparound, first-tick behavior, reset/clearHistory behavior, cross-evaluator consistency (GP tree vs bytecode), and index clamping
- **FR-014**: Both nodes MUST be opt-in via TrainingNodes configuration (not silently added to existing configurations)

### Key Entities

- **Distance History Buffer**: Third ring buffer in AircraftState alongside dPhiHistory_ and dThetaHistory_. Stores last 10 raw distance-to-rabbit values in meters. Shares the same historyIndex_, historyCount_, and timeHistory_ as existing buffers. Updated each tick by recording distance before GP evaluation. Cleared on new path via clearHistory().
- **GETDIST Node**: Nullary GP node (0 arguments). Returns current Euclidean distance to rabbit in meters.
- **GETDIST_PREV Node**: Unary GP node (1 argument: history index n). Returns buffered distance at index n. Clamped to [0, buffer_depth-1].
- **GETDIST_RATE Node**: Nullary GP node (0 arguments). Computes `(dist[0] - dist[1]) / dt` using timestamp-based dt. Clamped to [-10, 10] m/s. Returns 0.0 when insufficient history.

### Implementation Touchpoints

All files requiring updates when adding new GP nodes (derived from existing GETDPHI_PREV/RATE pattern):

| Category | File | What to update |
|----------|------|---------------|
| Opcode enum | autoc.h | Add GETDIST, GETDIST_PREV, GETDIST_RATE before _END |
| Opcode enum (portable) | gp_evaluator_portable.h | Mirror autoc.h enum additions |
| Node registration | autoc-eval.cc | Add to allNodes[] with name strings and arg counts |
| Evaluation dispatch | gp_evaluator_portable.cc | Add switch cases in evaluateGPOperator() |
| Function implementation | gp_evaluator_portable.cc | Implement executeGetDist(), executeGetDistPrev(), executeGetDistRate() |
| Function declarations | gp_evaluator_portable.h | Declare new execute functions |
| Bytecode dispatch | gp_evaluator_portable.cc | Add cases in evaluateBytecodePortable() |
| Bytecode extraction | gpextractor.cc | Add cases for child processing |
| Code generation | bytecode2cpp.cc | Update analyzeStackDepth(), getOperatorName(), generateInstruction() |
| History storage | aircraft_state.h | Add distHistory_[] buffer, extend recordErrorHistory(), clearHistory(), add getHistoricalDist() |
| History recording | minisim.cc | Compute and record distance before GP evaluation |
| Config (train) | autoc.ini | Update TrainingNodes, comment block, deprecate GETDTARGET |
| Config (eval) | autoc-eval.ini | Same as autoc.ini |
| Documentation | CLAUDE.md | Add node descriptions, deprecation note |
| Unit tests | tests/gp_evaluator_tests.cc | Add test cases for all three nodes |

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Evolution runs with GETDIST, GETDIST_PREV, GETDIST_RATE (and without GETDTARGET) produce GP trees that use the distance nodes in throttle control within 20 generations
- **SC-002**: Mean following distance in a 200-generation run with the new primitives is closer to 7.5m target than the current 20.7m baseline
- **SC-003**: Distance standard deviation (currently 7.1m) decreases, indicating reduced oscillation
- **SC-004**: All existing unit tests pass; new tests cover all three nodes across GP tree and bytecode evaluators with 100% pass rate
- **SC-005**: Bytecode evaluation produces identical per-step outputs to GP tree evaluation (zero divergence)
- **SC-006**: Existing GP trees/bytecode files using GETDTARGET continue to evaluate correctly (backward compatibility)

## Clarifications

### Session 2026-03-12

- Q: What position does GETDIST measure distance to? → A: Distance to rabbit's current interpolated position (same value used in fitness distance calculation). The discrete path interpolation is a simulation implementation detail — in real flight, a high-frequency estimator (~1000Hz) would provide much finer position estimates, so the GP should treat GETDIST as "estimated distance to target" regardless of update rate.

## Assumptions

- The simulation tick rate remains ~10Hz (100ms), consistent with existing temporal nodes
- Raw distance to rabbit is already computed each tick in the evaluation loop (used for fitness), so buffering adds negligible cost
- The AircraftState ring buffer can share timeHistory_ and index tracking with the existing dPhi/dTheta buffers (one index, one timestamp array, three value arrays)
- GETDIST_RATE clamping at [-10, 10] m/s is appropriate — aircraft closing at 10 m/s relative to rabbit is extreme (rabbit speed ~16 m/s, aircraft max ~25 m/s)
- GETDIST returns estimated distance to rabbit's current interpolated position (lookahead=0), not a future path point. In simulation this uses discrete path interpolation; in real flight a high-frequency estimator would provide finer resolution. The GP treats it as "distance to target" regardless of update rate.

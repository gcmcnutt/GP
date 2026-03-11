# Feature Specification: Entry Position Variations with Intercept-Budget Fitness Scaling

**Feature Branch**: `005-entry-fitness-ramp`
**Created**: 2026-03-10
**Status**: Draft
**Input**: User description: "Add entry condition variations including random position offsets, and scale distance fitness per-step based on estimated time-to-intercept budget so early large errors don't overwhelm the GP signal during intercept phase."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Intercept-Budget Distance Scaling (Priority: P1)

The researcher runs a GP evolution where the aircraft starts at a random position offset from the path start. The fitness function estimates a time-to-intercept budget from the initial displacement and orientation, then scales the per-step distance penalty proportionally to how far through that budget the current timestep is. Early in intercept (step is 1/7 of budget), distance penalty is scaled down to ~1/7. As the craft approaches the intercept point, scaling ramps toward full penalty. After intercept, full tracking penalty applies.

This lets early-generation GP individuals receive gradient signal for "heading toward the path" rather than being crushed by large absolute distance errors during the approach phase.

**Why this priority**: Without this scaling, entry position variations produce uniformly terrible fitness for all individuals in early generations — no selection gradient, no learning signal. This is the core mechanism that makes entry variations trainable.

**Independent Test**: Run evolution with position-offset entries. Verify that fitness values differentiate between individuals that turn toward the path vs those that fly away, even when both are far from the rabbit during intercept.

**Acceptance Scenarios**:

1. **Given** an aircraft starting 60m from the path with 90-degree heading offset, **When** the intercept budget is estimated at 7 seconds, **Then** at t=1s the distance penalty is scaled to approximately 1/7 of the full penalty value.
2. **Given** an aircraft that has completed its intercept budget, **When** it reaches the tracking phase, **Then** the distance penalty is at full (unscaled) value, identical to the current system.
3. **Given** two GP individuals where one turns toward the path and one flies away, **When** both start from the same offset entry, **Then** the one turning toward the path has meaningfully lower fitness (better).

---

### User Story 2 - Entry Position Variations (Priority: P2)

The researcher configures position offset variations in `autoc.ini`. The system generates random starting positions displaced from the path start point, in addition to existing heading/roll/pitch/speed offsets. The displacement is Gaussian-distributed and configurable via sigma parameters, following the same variation_generator pattern used for attitude offsets.

**Why this priority**: Position offsets create realistic intercept scenarios. Without them, the aircraft always starts at the path — attitude offsets alone aren't enough to train intercept behavior.

**Independent Test**: Run with position variation enabled. Verify aircraft start positions are displaced from path start according to configured sigma, and that ScenarioMetadata carries the offset to the simulator.

**Acceptance Scenarios**:

1. **Given** `EntryPositionSigma=30` configured, **When** scenarios are generated, **Then** starting positions are Gaussian-distributed around the path start with sigma of 30 meters (horizontal).
2. **Given** position offset variations, **When** sent to minisim/crrcsim, **Then** the simulator applies the position offset to the aircraft initial state before simulation begins.
3. **Given** position sigma of 0, **When** scenarios are generated, **Then** behavior is identical to current system (no position offset).

---

### User Story 3 - Intercept Budget Estimation (Priority: P2)

The system computes an approximate time-to-intercept from the initial conditions: distance to path start, heading offset relative to path tangent, and aircraft speed. The estimate is rough (1-digit accuracy is sufficient) — it accounts for turn time plus straight-line closure time plus additional time for the rabbit's motion.

**Why this priority**: The intercept budget drives the scaling function. An overestimate is acceptable (gentle ramp), an underestimate would premature apply full penalty. A simple geometric estimate is preferred over a complex one.

**Independent Test**: Compute intercept budget for known configurations (e.g., 60m offset, 90-degree heading, 20 m/s speed). Verify the estimate is within 50% of a reasonable value (geometric turn + closure + rabbit motion).

**Acceptance Scenarios**:

1. **Given** 60m displacement, 90-degree heading offset, 20 m/s speed, 15 m/s rabbit, **When** intercept budget is computed, **Then** the estimate is approximately 5-10 seconds (turn time ~2s + closure ~3s + rabbit compensation ~2-3s).
2. **Given** 0m displacement and 0-degree heading offset, **When** intercept budget is computed, **Then** the estimate is near zero (already intercepted), and full tracking penalty applies from the start.

---

### Edge Cases

- What happens when the intercept budget estimate is very large (>15s)? Cap at a maximum to prevent near-zero scaling for the entire run.
- What happens when the aircraft never intercepts (flies away for the whole scenario)? The ramp still applies — the distance penalty at budget fraction > 1.0 is full penalty, so the crash/completion penalty dominates as it does today.
- What happens with zero position offset but large heading offset? Budget estimation handles this — it's purely a turn-time intercept with small closure distance.
- What if the aircraft intercepts faster than estimated? The scaling only helps — steps after the budget point get full penalty, which matches actual tracking.

## Clarifications

### Session 2026-03-10

- Q: Should attitude delta penalty also be scaled during intercept window? → A: Yes, scale both distance AND attitude during intercept. Aggressive maneuvering (hard turns, attitude changes) is expected and desired during intercept phase.
- Q: How precise should the intercept budget estimation be? → A: Crude "hacktor" — e.g., geometric distance/speed + ~3s for a 90-degree turn + rabbit motion compensation. Path does NOT publish optimal approach attitude, so any direction is valid for intercept. Simple fudge-factor estimation is fine.
- Q: Per-path or per-scenario intercept budget? → A: Per-scenario. Each of the ~49 scenarios has its own path + entry conditions and gets its own intercept budget. Aggregate fitness is the sum across all scenarios, each with its own ramp. All inputs are deterministic.
- Q: Should there be a floor on the scale factor at t=0? → A: Yes, configurable floor (e.g., ~0.1) and ceiling (likely 1.0) as compile-time constants in header declarations, not autoc.ini params. Non-linear function ensures some signal always gets through but importance grows with estimated time progression.
- Q: Arena geometry and position offset bounds? → A: Arena is a cylinder centered on origin at SIM_INITIAL_ALTITUDE (-25m NED). Crash bounds: radius=SIM_PATH_RADIUS_LIMIT (70m), altitude SIM_MIN_ELEVATION (-7m) to SIM_MAX_ELEVATION (-120m). Entry positions must stay well inside (15m margin → safe radius ~55m, safe altitude -22m to -105m). Position offsets are relative to the test origin (with its -Z bias), generated via Gaussian params. Cylindrical coordinates (radius + angle + altitude offset) are natural for the cylindrical arena. Existing attitude variations can be expanded to all-attitude as training ramp matures — compile-time knobs for entry position envelope.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST estimate a per-scenario intercept time budget from the initial displacement, heading offset, aircraft speed, and rabbit speed.
- **FR-002**: System MUST scale BOTH the per-step distance penalty AND the attitude delta penalty by a function of `(current_step_time / intercept_budget)`, where the scaling is near zero at the start and ramps to 1.0 at the intercept budget time. Aggressive maneuvering during intercept is expected and should not be penalized.
- **FR-003**: The scaling function MUST be applied per-step BEFORE the power function (e.g., `pow(scale * distance / DISTANCE_NORM, DISTANCE_POWER)`), so it modulates the raw values.
- **FR-004**: Steps beyond the intercept budget time MUST receive full (unscaled) distance penalty, identical to current behavior.
- **FR-005**: System MUST add entry position offsets (North, East, Down) to ScenarioMetadata and VariationOffsets, following the existing Gaussian variation_generator pattern. Offsets are relative to the test origin (path start at SIM_INITIAL_ALTITUDE).
- **FR-006**: Position offset MUST be configurable via `autoc.ini` sigma parameters (e.g., `EntryPositionRadiusSigma`, `EntryPositionAltSigma`), with 0 meaning disabled. Cylindrical generation is natural for the cylindrical arena: Gaussian radius + uniform angle + Gaussian altitude offset.
- **FR-011**: Generated entry positions MUST stay within safe arena bounds (compile-time constants with ~15m margin from crash boundaries: max radius ~55m, altitude range approximately -22m to -105m NED). Positions exceeding bounds MUST be clamped.
- **FR-012**: Existing attitude variation sigmas (heading, roll, pitch) MUST be expandable toward all-attitude entries as the training ramp matures. The variation_generator and RAMP_LANDSCAPE infrastructure already support this — sigmas can be increased in autoc.ini.
- **FR-007**: The intercept budget estimation MUST handle the degenerate case of zero displacement/zero heading offset by producing a near-zero budget (immediate full tracking penalty).
- **FR-008**: The intercept budget MUST be capped at a configurable maximum to prevent near-zero scaling for entire runs.
- **FR-009**: The scaling function MUST be non-linear (e.g., quadratic or sigmoid) to provide some gradient signal early while ramping quickly toward full penalty near intercept.
- **FR-010**: The scaling function MUST have a configurable floor (minimum scale, e.g., ~0.1) and ceiling (maximum scale, e.g., 1.0), defined as compile-time constants in a header file alongside other fitness constants (DISTANCE_NORM, etc.), not as autoc.ini runtime parameters.

### Key Entities

- **Intercept Budget**: Estimated time in seconds for the aircraft to reach the path from its initial offset position and heading. Computed once per scenario at the start of each scenario's evaluation. Each of the ~49 scenarios in a generation gets its own budget based on its specific entry conditions.
- **Intercept Scale Factor**: Per-step multiplier on distance and attitude (floor to ceiling, e.g., 0.1 to 1.0) derived from `f(step_time / intercept_budget)`. Applied to raw values before normalization and power function. Floor/ceiling are compile-time constants.
- **Entry Position Offset**: North/East/Down displacement in meters from test origin (path start at SIM_INITIAL_ALTITUDE). Generated in cylindrical coordinates (Gaussian radius + uniform angle + Gaussian altitude), clamped to safe arena bounds. Added to ScenarioMetadata alongside existing heading/roll/pitch/speed offsets.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: GP evolution with 30m position-offset entries produces controllers that successfully intercept and track (fitness below crash threshold) within 50 generations.
- **SC-002**: Evolved controllers handle any entry position within the trained sigma envelope without crashes on 90%+ of test scenarios.
- **SC-003**: A single controller handles both intercept and tracking phases — no separate controller or mode switch needed.
- **SC-004**: Fitness values during intercept-phase training differentiate between good and bad intercept behavior (fitness variance across population > 10% in generation 1).
- **SC-005**: With position sigma of 0 and budget estimation disabled, fitness values are identical to the current system (backward compatibility).

## Assumptions

- The intercept budget is estimated geometrically, not simulated. A rough estimate (±50%) is acceptable because the scaling function is smooth.
- Position offsets include altitude (Down axis in NED) in addition to horizontal. Cylindrical generation (radius sigma + uniform angle + altitude sigma) fits the cylindrical arena naturally.
- Entry positions are clamped to safe arena bounds (~15m inside crash boundaries). Arena: cylinder radius 70m, altitude -7m to -120m NED.
- The scaling function applies to BOTH distance and attitude delta components during the intercept budget window.
- The path does NOT publish an optimal approach attitude — any intercept direction is valid, so the budget estimation does not assume a specific approach vector.
- The rabbit speed for budget estimation uses the nominal configured speed, not per-step variable speed.
- Minisim already supports entry heading/roll/pitch offsets via ScenarioMetadata; adding position offset follows the same pattern.

# Feature Specification: Variations Redux

**Feature Branch**: `003-variations-redux`
**Created**: 2026-03-09
**Status**: Draft
**Input**: Enable and validate entry/wind variations for robust sim-to-real transfer

## Clarifications

### Session 2026-03-09

- Q: What DISTANCE_NORM/POWER values for tightened fitness? → A: NORM=2.0, POWER=2.0 (aggressive tightening)
- Q: What throughput target when scaling to 36 scenarios? → A: No hard target; scale incrementally (1→9→36), adjust if too slow
- Q: Verify ScenarioMetadata serialization compatibility before enabling? → A: Trust it works; both repos on same Boost version. Compatibility is tracked in the Boost migration backlog item, not this feature.

## User Scenarios & Testing

### User Story 1 - Tighter Distance Tracking (Priority: P1)

The trainer produces controllers that track the rabbit more closely than the current ~15-20m steady-state distance. Tightening the fitness distance parameters makes the GP push harder toward the target before diminishing returns kick in.

**Why this priority**: Without tighter tracking, adding variations will make an already-loose controller diverge further. This is the baseline quality gate.

**Independent Test**: Run a 200-gen training with tightened DISTANCE_NORM/POWER and no variations enabled. Compare best-of-generation fitness and rendered tracking distance against the 002 baseline.

**Acceptance Scenarios**:

1. **Given** DISTANCE_NORM reduced and DISTANCE_POWER increased, **When** training completes 200 generations with no variations, **Then** the evolved controller tracks within 10m of the rabbit on straight segments (improved from ~15-20m baseline).
2. **Given** tightened distance params, **When** training completes, **Then** fitness converges monotonically and no regression in crash rate vs baseline.

---

### User Story 2 - Wind Direction Variations (Priority: P2)

The simulator applies wind direction offsets from ScenarioMetadata so the GP trains against varied wind conditions. The crrcsim integration is already complete — this story validates it end-to-end by enabling the config flag and running training.

**Why this priority**: Wind variation is the simplest variation to enable — it only requires a wind direction offset in crrcsim, no changes to aircraft launch. It exercises the full variation pipeline end-to-end.

**Independent Test**: Enable `EnableWindVariations=1` with a small number of scenarios (e.g., 4-9). Run 10 generations and verify via crrcsim console output that different wind directions are applied per scenario.

**Acceptance Scenarios**:

1. **Given** EnableWindVariations=1 and WindDirectionSigma=45, **When** a scenario is sent to crrcsim, **Then** crrcsim applies the wind direction offset from ScenarioMetadata to its base wind direction.
2. **Given** varied wind scenarios, **When** training completes 200 generations, **Then** fitness improves (downward trend) despite harder conditions.

---

### User Story 3 - Entry Condition Variations (Priority: P3)

The simulator applies entry heading, roll, pitch, and speed offsets from ScenarioMetadata at aircraft launch. The crrcsim integration is already complete — this story validates it by enabling the config flag and training with off-nominal entry conditions to evolve recovery behaviors.

**Why this priority**: Entry variations are the core robustness feature. Build on wind-only validation first to isolate issues.

**Independent Test**: Enable `EnableEntryVariations=1` alongside wind variations. Run 10 generations and verify via crrcsim logging that aircraft launch with varied heading/roll/pitch/speed per scenario.

**Acceptance Scenarios**:

1. **Given** EnableEntryVariations=1 with sigma values (heading=45deg, roll=22.5deg, pitch=7.5deg, speed=10%), **When** a scenario launches in crrcsim, **Then** the aircraft starts with the specified offsets applied to heading, roll, pitch, and speed.
2. **Given** entry+wind variations enabled, **When** training completes 200 generations, **Then** evolved controllers recover from 45deg roll entries without crashing in >80% of scenarios.

**Open Design Issues for Entry Variations**:

- **Entry position variation**: Current entry is always near path start. Real-world entry would be anywhere in the arena. Consider adding entry position offsets to ScenarioMetadata so the controller must intercept from arbitrary positions, not just attitude recovery from a fixed launch point.
- **Intercept phase fitness ramp**: Early simulation steps penalize the controller for being far from the rabbit, but with off-nominal entry (especially varied position), the aircraft needs an intercept phase before tracking begins. The fitness penalty for distance should ramp from 0 to full over an estimated intercept window rather than counting early distant steps at full penalty. This window could be dynamic — computed from initial distance to path at launch. Without this, the GP is incentivized to minimize early-step error (which it can't control) rather than developing a good intercept-to-track transition.

---

### User Story 4 - Variable Rabbit Speed (Priority: P4)

Enable non-zero RabbitSpeedSigma so the target moves at varying speeds. This is already fully implemented in autoc — just needs config enablement and training validation.

**Why this priority**: Lowest risk since code is complete. Enable after entry/wind are validated.

**Independent Test**: Set RabbitSpeedSigma>0, run training, verify path timestamps vary per scenario.

**Acceptance Scenarios**:

1. **Given** RabbitSpeedSigma=2.0 with min=8/max=25 m/s, **When** training runs, **Then** rabbit speed varies smoothly between scenarios and the GP evolves speed-adaptive throttle management.

---

### User Story 5 - Progressive Variation Ramp (Priority: P5)

The variation ramp (VariationRampStep) scales variation intensity from 0 to full over the course of training, easing the GP into harder conditions.

**Why this priority**: Tuning knob for training convergence. Validate after individual variations work.

**Independent Test**: Run with VariationRampStep=5 and verify via autoc logs that variation scale increases from 0% to 100% over training generations.

**Acceptance Scenarios**:

1. **Given** VariationRampStep=5 and 200 generations, **When** training runs, **Then** early generations see near-zero variations and later generations see full sigma variations.
2. **Given** ramped variations, **When** compared to non-ramped (full sigma from start), **Then** ramped training converges to equal or better fitness.

---

### Edge Cases

- What happens when entry variations produce extreme attitudes (e.g., 3-sigma = 67.5deg roll)? Expected: aircraft may crash, crash penalty provides gradient.
- What happens when wind offset wraps around 360deg? Expected: crrcsim normalizes to [0, 2pi).
- What happens when entry speed factor is very low (e.g., 0.8x = near stall)? Expected: GP learns throttle-up response or crashes with completion penalty.
- What happens when all variations are enabled simultaneously from gen 0 (no ramp)? Expected: slower convergence but eventual fitness improvement.

## Requirements

### Functional Requirements

- **FR-001**: crrcsim MUST apply wind direction offset from ScenarioMetadata at scenario initialization, adding the offset to the base wind direction from autoc_config.xml.
- **FR-002**: crrcsim MUST apply entry heading offset from ScenarioMetadata at aircraft launch, rotating the launch heading.
- **FR-003**: crrcsim MUST apply entry roll and pitch offsets from ScenarioMetadata as initial aircraft attitude at launch.
- **FR-004**: crrcsim MUST apply entry speed factor from ScenarioMetadata as a multiplier on launch velocity.
- **FR-005**: crrcsim MUST log applied variation values to console for debugging when variations are non-zero.
- **FR-006**: The fitness function distance parameters MUST be set to DISTANCE_NORM=2.0, DISTANCE_POWER=2.0 (from current 5.0/1.5) to produce tighter tracking.
- **FR-007**: The variation ramp MUST scale all variation offsets from 0 to full sigma over training (already implemented in autoc).
- **FR-008**: Training MUST complete 200 generations with all variations enabled without crashes or hangs.
- **FR-009**: WindScenarios count MUST be configurable up to 36 scenarios per generation. No hard throughput target — scale incrementally (1→9→36) and adjust if generations become too slow.

### Key Entities

- **ScenarioMetadata**: Carries variation offsets (entryHeadingOffset, entryRollOffset, entryPitchOffset, entrySpeedFactor, windDirectionOffset) from autoc to crrcsim. Already version 5 with these fields.
- **VariationSigmas**: Configuration struct defining Gaussian distribution widths for each variation type. Already implemented in variation_generator.h.
- **Fitness Constants**: DISTANCE_NORM, DISTANCE_POWER in autoc.h controlling the nonlinear tracking penalty.

## Scope Boundaries

### In Scope
- Validation of existing crrcsim entry and wind variation offset integration
- Fitness distance parameter tuning
- Progressive training validation (baseline -> wind -> entry+wind -> variable rabbit -> ramp)
- Scenario count scaling (1 -> 9 -> 36)

### Out of Scope (Deferred)
- Layered controller / safety layer (separate feature)
- Craft parameter variations: CG, wing loading, drag, power/thrust (future VARIATIONS2)
- Code unification of GP/bytecode eval paths (UNIFY.md)
- LongSequential path Immelman fix
- Control command smoothness penalty (GP evolves smooth naturally with pop=20K + temporal nodes)
- Target pose estimation/interpolation

## Assumptions

- The autoc side variation generation is complete and correct (validated in VARIATIONS1 work).
- ScenarioMetadata serialization (Boost binary, version 5) is compatible between current autoc and crrcsim builds. Both repos currently use the same Boost version; cross-version compatibility is tracked separately in the Boost migration backlog item.
- crrcsim's aircraft launch and wind initialization code has accessible points to apply offsets.
- Current population size (20K) and node set (23 nodes) are sufficient for the expanded training scenarios.

## Success Criteria

### Measurable Outcomes

- **SC-001**: Baseline (no variations) tracking distance improves from ~15-20m to <10m on straight path segments with tightened fitness parameters.
- **SC-002**: Training with wind variations enabled converges to fitness below 25K within 200 generations.
- **SC-003**: Training with entry+wind variations produces controllers that survive >80% of 2-sigma entry conditions (45deg heading, 22.5deg roll).
- **SC-004**: Full variation training (entry+wind+variable rabbit) completes 200 generations without system crashes or hangs.
- **SC-005**: Evolved controllers from variation training show qualitatively different behavior than baseline — specifically, use of GETROLL_RAD/GETPITCH_RAD/GETVEL nodes for recovery from off-nominal entries.

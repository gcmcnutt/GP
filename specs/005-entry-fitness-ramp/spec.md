# Spec: Entry Variations with Ramped Fitness Scaling

**Feature ID**: 005-entry-fitness-ramp
**Status**: Draft
**Created**: 2026-03-10

## Overview

Add entry condition variations (heading, roll, pitch, speed offsets) to training scenarios with a ramped fitness scale that progressively increases difficulty during evolution. This builds robustness for the intercept phase of flight.

## Problem Statement

Current training starts every scenario from a well-aligned entry. Real flights (evidence from Flight 25a/b) show the aircraft often arrives at the tracking phase with off-nominal attitude, speed, or heading. The GP controller must handle these gracefully.

Additionally, throwing full-difficulty variations at generation 0 overwhelms the population before it can learn basic tracking. A ramp schedule lets the population build competence before facing harder conditions.

## Scope

### Entry variations to implement:
- Heading offset (yaw error at intercept)
- Roll offset (banked entry)
- Pitch offset (climbing/diving entry)
- Speed offset (fast/slow arrival)

### Fitness ramp:
- Generation-dependent scaling of variation magnitude
- Start with near-nominal entries, progressively widen distribution
- Configurable ramp schedule (linear, sigmoid, or step)
- Ramp applies to entry variations independently from wind/speed variations

## Success Criteria

- GP controllers handle off-nominal entry conditions
- Fitness improves monotonically despite increasing difficulty (ramp is gradual enough)
- Configurable ramp schedule via autoc.ini
- No regression on nominal-entry scenarios

## Prior Art

- **Reference**: `autoc/specs/ZZZ-VARIATIONS1.md` (entry variation architecture, ScenarioMetadata v5)
- 003-variations-redux validated wind+speed variations approach
- RAMP_LANDSCAPE concept from 002-path-interpolation

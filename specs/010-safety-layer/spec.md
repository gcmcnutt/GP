# Spec: Layered Controller Safety Layer

**Feature ID**: 010-safety-layer
**Status**: Draft
**Created**: 2026-03-10

## Overview

Add hardcoded safety constraints that override GP control outputs when danger is detected, enabling safe exploration of aggressive all-attitude maneuvers.

## Problem Statement

The current tracking layer follows the rabbit but has no awareness of:
- Imminent crash (altitude too low)
- Out-of-bounds flight (geofence violation)
- Stall conditions (speed too low)
- Recovery vs tracking mode distinction

Without safety overrides, GP controllers crash when exploring inverted flight or aggressive maneuvers during evolution.

## Safety Constraints

- **Altitude protection**: if too low, force climb
- **Geofence**: if approaching boundary, break off and head toward origin
- **Stall prevention**: if too slow, force dive + throttle
- **Re-engagement**: when near origin/safe, resume tracking

## Implementation

Post-GP safety clamping in `evaluateOutput()` or similar — overrides GP outputs only when constraints are violated.

## Testing

- Create intentional OOB paths to verify safety layer behavior
- Paths that go toward boundaries, into the ground, etc.
- Verify break-off and re-engagement logic works

## Success Criteria

- GP can explore aggressive maneuvers without crashing
- Safety overrides activate only when needed
- Normal tracking behavior unaffected when constraints are satisfied

## Prior Art

- **Reference**: `autoc/specs/LAYERED_CONTROLLER.md` (full design, Phase 2)

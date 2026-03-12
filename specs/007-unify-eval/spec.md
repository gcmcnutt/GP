# Spec: Unify Evaluation Pipelines

**Feature ID**: 007-unify-eval
**Status**: Draft
**Created**: 2026-03-10

## Overview

Extract and unify the ~350 lines of duplicated evaluation code between `MyGP::evalTask()` (GP tree mode) and `BytecodeEvaluationGP::evalTask()` (bytecode mode) into shared components.

## Problem Statement

Two parallel evaluation pipelines maintain identical fitness computation, scenario metadata building, and logging logic. Any change to fitness, path handling, or eval data must be made in both places — a maintenance burden and divergence risk. This must be unified before adding new fitness logic (entry ramp), safety layers, or GPU offload.

## Proposed Phases

1. Extract `fitness_computer` — shared fitness calculation
2. Extract `eval_data_builder` — shared scenario metadata construction
3. Extract `eval_logger` — shared logging/diagnostics
4. Unify core eval loop — single `evalTask()` parameterized by evaluation backend

## Success Criteria

- Single fitness computation path for both GP tree and bytecode modes
- Both modes produce identical results (existing bytecode verification tests pass)
- No duplicated eval logic between the two code paths
- Clean extension point for future GPU evaluation backend

## Prior Art

- **Reference**: `autoc/specs/UNIFY.md` (detailed analysis, phase plan)
- Bug fixes documented: temporal history wipe, -ffast-math NaN guard, path interpolator float precision

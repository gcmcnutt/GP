# Spec: GP Library Fitness Serialization Precision

**Feature ID**: 006-fitness-precision
**Status**: Draft (GP-specific, not applicable to NN fork)
**Created**: 2026-03-10

## Overview

Fix fitness value serialization to preserve full double precision, so renderer and analysis tools show accurate fitness values.

## Problem Statement

Renderer shows fitness as `164073.000` instead of `164073.136503`. Root cause: `GP::save()` in `src/gp.cc:148-157` uses default ostream precision (6 digits), truncating double-precision fitness values.

## Solution

Add `os << std::setprecision(17) << stdFitness` (17 = `std::numeric_limits<double>::max_digits10`) in the GP save path.

## Scope

- Small, focused fix in `src/gp.cc`
- May also need to audit other serialization points for similar truncation
- Verify renderer displays full precision after fix

## Success Criteria

- Fitness values round-trip through save/load with full double precision
- Renderer shows accurate fractional fitness values

# Spec: Fix Immelman Path Segment

**Feature ID**: 009-immelman-path
**Status**: Draft
**Created**: 2026-03-10

## Overview

Fix the last segment of `longSequentialPath` generation so the Immelman turn connects smoothly, enabling all-attitude training paths.

## Problem Statement

The final path section of `AeroStandard::longSequentialPath` has discontinuities due to coordinate convention issues. The Immelman segment doesn't connect smoothly to the preceding path, creating unrealistic training scenarios.

## Success Criteria

- Immelman segment connects smoothly to preceding path (no position/attitude discontinuities)
- Path is flyable in simulation without crashes at transition points
- Correct NED coordinate conventions throughout

## Location

- `autoc/pathgen.h` (`AeroStandard::longSequentialPath`)

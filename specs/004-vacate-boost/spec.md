# Spec: Vacate Boost Serialization

**Feature ID**: 004-vacate-boost
**Status**: Draft
**Created**: 2026-03-10

## Overview

Replace Boost serialization with a portable, version-stable alternative for all RPC and artifact serialization in the autoc system.

## Problem Statement

Boost serialization is NOT portable across different Boost versions. Archive version headers change between releases (v18 vs v19), breaking cross-architecture artifact sharing via S3. Current issue: Boost 1.74 (x86) cannot read archives from Boost 1.83 (ARM64). This blocks the DGX-to-xiao workflow and any multi-platform development.

Profiling shows serialization is only ~1.5% of runtime, so this is a **portability fix**, not a performance fix.

## Scope

### Data structures requiring migration:
- EvalData (path lists, scenario metadata)
- EvalResults (aircraft state sequences, fitness)
- Path / AircraftState
- ScenarioMetadata
- GPBytecodeInterpreter (binary serialization for RPC)

### Alternatives to evaluate:
- cereal (header-only, similar API to Boost)
- Custom binary protocol (full control, minimal dependencies)
- FlatBuffers (zero-copy, schema-evolution)
- JSON/MessagePack (human-readable option)

## Success Criteria

- All serialized data readable across Boost versions and architectures
- No Boost.Serialization dependency in autoc
- RPC communication between autoc and minisim works cross-platform
- S3 artifacts readable from any build environment
- No regression in serialization performance

## Prior Art

- **Reference**: `autoc/specs/REPLACE_BOOST.md` (detailed analysis, protocol design)
- Path caching optimization identified as high-value (saves 90%+ of request size after first eval)

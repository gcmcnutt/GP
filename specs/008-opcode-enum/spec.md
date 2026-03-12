# Spec: Refactor Opcode Enum

**Feature ID**: 008-opcode-enum
**Status**: Draft
**Created**: 2026-03-10

## Overview

Consolidate the opcode enum to a single authoritative definition and eliminate numeric opcode usage in the bytecode2cpp generator.

## Problem Statement

Opcodes are currently defined/referenced in multiple places, with bytecode2cpp using numeric constants instead of the enum. This creates maintenance risk and makes it easy to introduce opcode mismatches.

## Solution

- Single `enum class Opcode` definition used by all components
- bytecode2cpp generates code using symbolic enum names
- Compile-time guarantee that all components agree on opcode values

## Success Criteria

- One definition of opcodes, used everywhere
- bytecode2cpp output uses symbolic names
- No numeric opcode literals in generated code

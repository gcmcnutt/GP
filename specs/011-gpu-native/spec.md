# Spec: GB10 GPU Native Evaluation

**Feature ID**: 011-gpu-native
**Status**: Draft
**Created**: 2026-03-10

## Overview

Move physics simulation evaluation to the GB10 GPU for 10-100x throughput improvement, enabling massively parallel scenario evaluation.

## Problem Statement

Current evaluation is CPU-bound, with physics simulation dominating runtime. The DGX workstation has a GB10 GPU that could run thousands of simulations in parallel. At GPU scale, the bottleneck shifts from simulation to coordination (Amdahl's law).

## Key Considerations

- Batch evaluations (64-256 per GPU batch)
- Move computation to worker side, reduce round-trips
- Depends on: unified eval pipeline (007), clean serialization (004)
- Physics model (EOM integration, aerodynamics) must be ported to CUDA/compute shaders
- GP tree evaluation could also run on GPU (uniform stack-based bytecode is GPU-friendly)

## Success Criteria

- 10x+ throughput improvement over CPU-only evaluation
- Correct physics results matching CPU implementation
- Scalable to full scenario space (paths x winds x craft x entry x speeds)

## Prior Art

- **Reference**: `autoc/specs/SCALEUP.md` (architecture, profiling results, batching strategy)
- Profiling: physics simulation dominates, serialization only ~1.5%, fitness ~25-30%

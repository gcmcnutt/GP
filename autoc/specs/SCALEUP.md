# SCALEUP: Scaling AutoC for GPU-Accelerated Simulation

## Executive Summary

Current architecture processes ~1200 sims/sec with 8 CPU workers (PERFORMANCE_BUILD with `-ffast-math`).
Determinism issues have been resolved. To scale for GPU-accelerated simulation (potentially 10-100x
throughput), we need to address:

1. **Task Granularity**: Each task = 1 GP individual × 6 paths (too fine-grained for GPU batching)
2. **Serialization Overhead**: Boost binary archives for every RPC call
3. **Physics Engine**: CPU-bound simulation could be ported to SIMD/CUDA

---

## Current Performance (January 2025)

| Build Config | Flags | Sims/sec | Deterministic |
|--------------|-------|----------|---------------|
| Default (Debug) | `-g -O0` | ~240 | Yes |
| PERFORMANCE_BUILD | `-O3 -march=native -flto -ffast-math` | ~1200 | Yes |

**Key achievement**: Determinism is now preserved with full optimizations enabled, thanks to the
headless mode fix in `SimStateHandler.cpp` that bypasses the wall-clock-dependent `dDeltaT`
accumulator when `video.enabled=0`.

---

## Resolved Issues

### Simulation Time Non-Determinism (Fixed)

**Root cause**: The `idle()` function in `SimStateHandler.cpp` used a static `dDeltaT` accumulator
for frame timing correction that persisted across simulation resets. This caused `simTimeMsec` to
vary between evaluations of the same GP, leading to different GP commands at the same physics step.

**Fix**: In headless mode (`video.enabled=0`), bypass the `dDeltaT` correction entirely:
```cpp
bool headlessMode = (cfgfile->getInt("video.enabled", 1) == 0);
if (headlessMode) {
    // Deterministic: fixed multiloop based on nDeltaTicks and dt
    multiloop = (int)(nDeltaTicks/1000.0/Global::dt + 0.5);
} else {
    // Video mode: apply frame timing correction for smooth display
    multiloop = (int)((nDeltaTicks/1000.0 - dDeltaT)/Global::dt + 0.5);
}
```

### `-ffast-math` and Determinism

With the simTime fix in place, `-ffast-math` no longer causes observable divergence in elite
re-evaluations. The floating-point reordering it enables is now safe because:
1. The same physics steps produce the same simTime
2. GP commands are therefore deterministic
3. Physics state evolution is deterministic

---

## Future Scaling Directions

### Option A: SIMD Vectorization of Physics Engine

The crrcsim physics engine (LARCsim FDM) is CPU-bound. Opportunities:
- Vectorize force/moment calculations across multiple scenarios
- SIMD quaternion math for attitude integration
- Batch evaluation of GP bytecode across multiple aircraft

### Option B: GPU (CUDA) Port

For 10-100x throughput, port core simulation to GPU:
- Parallel simulation of 1000+ aircraft simultaneously
- GPU-resident GP bytecode interpreter
- Batch RPC protocol (one call per generation, not per individual)

### Batch Protocol Design (Future)

```cpp
struct BatchEvalRequest {
  std::vector<std::vector<char>> gpPrograms;  // N GPs serialized
  std::vector<std::vector<Path>> pathList;    // Shared paths
  std::vector<ScenarioMetadata> scenarios;
};

struct BatchEvalResponse {
  std::vector<double> fitnessValues;          // N fitness values
  std::vector<CrashReason> crashReasons;
};
```

**Benefits**:
- Amortizes RPC overhead across thousands of GPs
- GPU can batch-compile bytecode
- Only return detailed traces for elite candidates

---

## Scenario Caching (Future Optimization)

### Current: Redundant Path Transmission

Every EvalData message includes full path geometry (~67 KB per eval).
At 9000 evals/gen = 600 MB of redundant path data per generation.

### Proposed: Scenario Registry

Cache paths at generation start, reference by ID during evaluation:
```cpp
struct EvalDataCompact {
  std::vector<char> gp;      // ~1-2 KB only
  uint64_t gpHash;
  uint32_t scenarioId;       // Reference to cached paths
};
```

**Savings**: 67 KB → 2 KB per eval = 97% reduction in RPC payload.

---

## Metrics Summary

| Metric | Current (CPU) | Target (GPU) |
|--------|---------------|--------------|
| Sims/sec | 1,200 | 10,000-60,000 |
| Worker CPU util | 80-99% | N/A |
| RPC calls/gen | 9,000 | 1-10 |
| Deterministic | Yes | Yes |

---

## Key Source Locations

| Component | File | Notes |
|-----------|------|-------|
| Headless mode fix | SimStateHandler.cpp:67-96 | dDeltaT bypass |
| evalTask hot path | autoc.cc | GP evaluation |
| ThreadPool | threadpool.h | Task dispatch |
| RPC protocol | minisim.h | EvalData/EvalResults |
| crrcsim worker | inputdev_autoc.cpp | Simulation loop |
| PERFORMANCE_BUILD | CMakeLists.txt | Compiler flags |

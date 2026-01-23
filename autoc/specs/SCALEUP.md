# SCALEUP: Scaling AutoC for GPU-Accelerated Simulation

## Executive Summary

Current architecture processes ~1200 sims/sec with 8 CPU workers (PERFORMANCE_BUILD with `-ffast-math`).
Determinism issues have been resolved. To scale for GPU-accelerated simulation (potentially 10-100x
throughput), we need to address:

1. **Task Granularity**: Each task = 1 GP individual × 6 paths (too fine-grained for GPU batching)
2. **Serialization Overhead**: Boost binary archives for every RPC call
3. **Physics Engine**: CPU-bound simulation could be ported to SIMD/CUDA

---

## UPDATE LOG (Jan 2026)

**Target Revision:** 50,000-100,000+ sims/sec for scaled-up scenario variations (craft parameters, wind, etc.)

**Key Architecture Insight:** At GPU scale, the bottleneck shifts from simulation to coordination:
- GPU can simulate 10,000+ aircraft in parallel (~1ms per timestep batch)
- But autoc dispatch/serialize/deserialize may become Amdahl's law limit
- Solution: Move more computation to worker side, reduce round-trips

**Proposed Architecture:** Single ultra-parallel worker per GPU (not multi-threaded CPU workers)
- autoc sends batch of 1000+ GP programs per generation
- Worker simulates entire batch in parallel on GPU
- Returns only fitness values (not full state traces)

**Profiling Required First:**
1. Where does time go in autoc? (task dispatch, fitness calc, evolution, serialization)
2. Where does time go in crrcsim? (FDM physics, wind model, aero coefficients)
3. What fraction is parallelizable? (Amdahl's law)

---

## PROFILING RESULTS (Jan 2026)

Flame graph profiling completed with `perf record -g --call-graph dwarf`. Key findings:

### Overall Time Distribution

| Component | % Time | Notes |
|-----------|--------|-------|
| **crrcsim (total)** | ~85% | Dominates runtime |
| **autoc (total)** | ~2% | Dispatcher overhead very low |
| **autoc ThreadPool/serialize** | ~1.9% | Boost serialization visible but small |

### crrcsim Breakdown (within the 85%)

| Function | % Total | Category | GPU Candidate? |
|----------|---------|----------|----------------|
| **CalculateWindGust (Dryden)** | 15.6% | Wind/turbulence | Yes - embarrassingly parallel |
| **gear()/WheelSystem** | 13.7% | Ground contact | **Disable for gliders** |
| **GP evaluation (MyGene::evaluate)** | ~14% | Controller | Yes |
| **Thermal::update** | 7.0% | Wind/turbulence | Yes |
| **CalculateWindGrad** | 4.9% | Wind/turbulence | Yes |
| **aero()** | 4.9% | Aerodynamics | Yes |
| **ls_aux()** | 5.6% | EOM helpers | Yes |
| **engine()** | 5.0% | Propulsion | Yes |
| **ls_step()** | 1.5% | EOM integration | Yes |
| **ls_accel()** | 1.0% | EOM acceleration | Yes |
| Matrix33 ops (scattered) | ~10% | Linear algebra | SIMD target |

### Comparison with Pre-Profiling Estimates

| Component | PROFILING.md Expected | Actual | Delta |
|-----------|----------------------|--------|-------|
| aero() | 30-50% | **4.9%** | ❌ 6-10x lower |
| EOM (ls_step + aux + accel) | 20-30% | **8.1%** | ❌ 2-3x lower |
| Wind/turbulence | 10-20% | **27.5%** | ❌ 1.5x higher |
| GP evaluation | 10-20% | **14%** | ✓ Within range |
| Gear/wheels | Not predicted | **13.7%** | ❌ Unexpected! |
| Boost serialization | 10-30% | **~1.5%** | ❌ 10x lower |

### Key Insights

1. **Wind/Dryden turbulence is the primary bottleneck** (27.5% combined)
   - `calculate_gust()` at 14% includes expensive `Matrix33::operator*` (8% alone)
   - Matrix multiplication dominates: prime SIMD/GPU target
   - `RandGauss::Get()` calls at 2.3% (Dryden turbulence RNG)

2. **Gear calculations unexpectedly high** (13.7%)
   - `Wheel::update()` doing lots of matrix math
   - **Aircraft-specific**: gliders don't need wheel physics in flight
   - **Action**: Add config knob to disable gear physics for wheeless craft

3. **Boost serialization NOT a performance bottleneck** (~1.5%)
   - REPLACE_BOOST priority lowered from performance perspective
   - Still valuable for cross-platform portability

4. **Amdahl's Law is very favorable**
   - Serial fraction: ~2% (autoc overhead)
   - Parallel fraction: ~98% (crrcsim simulation)
   - Theoretical max speedup: 50x before serial overhead dominates
   - **Implication**: GPU scaling should achieve near-linear speedup

5. **aero() is surprisingly cheap** (4.9%)
   - Pre-profiling estimate was 30-50%
   - Simple polynomial coefficient lookups, well-optimized by compiler
   - Not the primary GPU port target

### Optimization Priority Matrix

| Target | % Time | Effort | Impact | Priority |
|--------|--------|--------|--------|----------|
| Disable gear for gliders | 13.7% | Low | High | **1** |
| SIMD Matrix33 ops | ~10% | Medium | High | **2** |
| Wind model caching | 27.5% | Medium | High | **3** |
| GPU port of wind+physics | 50%+ | High | Very High | **4** |

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

---

## Physics Engine Analysis (Jan 2026)

### crrcsim FDM Architecture

The flight dynamics model (FDM) is based on NASA LARCsim with the following structure:

```
fdm_larcsim::update()
├── env->CalculateWind()       # Wind field lookup
├── env->CalculateWindGrad()   # Wind gradient for stability derivatives
├── for n in 0..multiloop:     # Sub-stepping loop (typically 1-10 iterations)
│   ├── ls_step(dt)            # Equations of motion integration (Adams-Bashforth)
│   │   ├── Linear velocity integration
│   │   ├── Angular velocity integration
│   │   ├── Quaternion integration + normalization
│   │   └── Position update
│   ├── CalculateWindGust()    # Dryden turbulence model
│   ├── ls_aux()               # Aerodynamic state (alpha, beta, airspeed)
│   ├── aero()                 # Force/moment coefficients → forces
│   ├── engine()               # Propulsion model
│   ├── gear()                 # Ground contact
│   └── ls_accel()             # Sum forces → accelerations
```

### Computational Breakdown (MEASURED Jan 2026)

| Component | Measured % | GPU-Friendly? | Notes |
|-----------|------------|---------------|-------|
| **CalculateWindGust()** | 15.6% | Yes | Dryden turbulence, Matrix33 ops dominate |
| **gear()** | 13.7% | **Disable** | Not needed for gliders/aircraft without wheels |
| **GP eval** | 14% | Yes | Bytecode interpreter, recursive calls |
| **Thermal::update()** | 7.0% | Yes | Updraft model |
| **ls_aux()** | 5.6% | Yes | Aerodynamic state computation |
| **engine()** | 5.0% | Yes | Battery/propeller model |
| **CalculateWindGrad()** | 4.9% | Yes | Wind gradient for stability |
| **aero()** | 4.9% | Yes | Force/moment coefficients |
| **ls_step()** | 1.5% | Yes | Integration step |
| **ls_accel()** | 1.0% | Yes | Acceleration computation |

**Key finding:** Wind/turbulence (27.5%) > Gear (13.7%) > GP eval (14%) > aero (4.9%)

**Surprise:** aero() is much cheaper than expected (4.9% vs 30-50% estimated)

### Parallelism Opportunities

1. **Aircraft-level parallelism** (primary target)
   - Each aircraft is independent during simulation
   - Perfect for GPU: 10,000 aircraft × 1000 timesteps = 10M independent work items
   - Memory: ~1KB state per aircraft = 10MB for 10K aircraft (fits in GPU cache)

2. **Timestep-level parallelism** (limited)
   - Sequential dependency: state(t+1) = f(state(t))
   - Can't parallelize within single aircraft trajectory
   - But CAN pipeline: while aircraft[i] does step t+1, aircraft[i+1] does step t

3. **Path-level parallelism** (within single eval)
   - 6 paths per eval can run in parallel
   - Already exploited by current multi-worker architecture

### Amdahl's Law Analysis (MEASURED Jan 2026)

Profiling results show much better parallelism than pre-estimated:

```
MEASURED time breakdown:
- autoc ThreadPool + boost serialize/deserialize: ~2% (serial)
- crrcsim simulation (all physics + GP eval):     ~85% (parallel)
- Unaccounted (scheduling, I/O, etc):             ~13%

Parallel fraction: ~98%
Serial fraction: ~2%

Amdahl's limit with infinite parallelism: 1/0.02 = 50x speedup
```

**Implication:** Serial overhead is NOT the bottleneck. GPU scaling should achieve
near-linear speedup up to ~50x before autoc becomes the limit.

At 1200 sims/sec baseline: theoretical max = 60,000 sims/sec (vs 5,000 pre-estimated)

**This is MUCH better than expected.** The focus should be on:
1. Quick wins (disable gear for gliders: +14% instantly)
2. SIMD vectorization of Matrix33 ops
3. GPU port of wind model (primary target)

### Reducing Serial Overhead (Lower Priority Post-Profiling)

**Note:** Profiling shows serial overhead is only ~2%, not 24% as originally estimated.
These optimizations are still valuable for GPU scale but are no longer critical path.

1. **Batch protocol** (see REPLACE_BOOST.md)
   - Send all 1000 GP programs in one RPC call
   - Return only fitness values (not full state traces)
   - **Impact reduced**: Original estimate was 10x improvement, actual ~1.2x

2. **Worker-side fitness computation**
   - Worker computes fitness, not autoc
   - Still valuable for reducing bandwidth
   - **Impact reduced**: Fitness calc is small fraction of autoc time

3. **Async evolution pipelining**
   - Start gen N+1 simulation while gen N fitness is being processed
   - **Impact reduced**: Overlap is already good

**Revised priority:** Focus on simulation speedup (wind model, gear disable) over protocol optimization.

---

## GPU/TPU Scaling Strategy

### Phase 1: Profile Current System

**See [PROFILING.md](PROFILING.md) for detailed flame graph profiling guide.**

Quick start:
```bash
# Build with profiling flags (already configured in CMakeLists.txt)
cd ~/GP && rm -rf build && mkdir build && cd build
cmake -DENABLE_PROFILING=ON ../autoc && make -j$(nproc)

cd ~/crsim/crrcsim-0.9.13 && rm -rf build && mkdir build && cd build
cmake -DENABLE_PROFILING=ON .. && make -j$(nproc)

# Profile with single worker for clean results
cd ~/GP/build
perf record -g -F 99 --call-graph dwarf -o profile.data \
    ./autoc  # use config with evalThreads=1

# Generate flame graph
perf script -i profile.data | \
    ~/FlameGraph/stackcollapse-perf.pl | \
    ~/FlameGraph/flamegraph.pl > flame.svg
```

**Expected findings:** (VALIDATED - see PROFILING RESULTS section above)
- autoc: serialization (~1.5%), threadpool dispatch (~0.5%) - **MUCH LOWER than expected**
- crrcsim: wind/turbulence (27.5%), gear (13.7%), GP eval (14%), aero (4.9%) - **wind dominates, aero is small**

### Phase 2: Batch Protocol Implementation

See REPLACE_BOOST.md for wire protocol details. Key changes:

1. **EvalBatch replaces EvalData**
   - Contains 64-256 GP programs per batch
   - Single path set (cached by worker)
   - Single scenario list (shared across batch)

2. **EvalBatchResults replaces EvalResults**
   - Only fitness values (no state traces)
   - Crash reasons per eval
   - Optional: Elite indices for detailed re-eval

3. **Worker-side fitness computation**
   - Move `computePathFitness()` to worker
   - Requires sharing fitness weights via config (not per-message)

### Phase 3: CUDA Port of Physics Engine

**Target: 10,000-50,000 sims/sec on single GPU**

#### Kernel Structure

```cuda
// Each thread simulates one aircraft for one timestep
__global__ void fdm_step_kernel(
    AircraftState* states,      // [N_AIRCRAFT]
    const GPBytecode* programs, // [N_AIRCRAFT]
    const Path* paths,          // [N_PATHS] (shared)
    const WindConfig* wind,     // Single config
    float dt,
    int step_index
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N_AIRCRAFT) return;

    AircraftState& state = states[i];

    // 1. Wind calculation (can use shared memory for wind tables)
    float3 wind_vel = calculate_wind(state.pos, wind);

    // 2. GP evaluation (stack-based bytecode interpreter)
    GPContext ctx = {state, paths, step_index};
    float pitch, roll, throttle;
    eval_gp(programs[i], ctx, &pitch, &roll, &throttle);

    // 3. Physics integration
    state = integrate_step(state, pitch, roll, throttle, wind_vel, dt);
}

// Host function: simulate N aircraft for M timesteps
void simulate_batch(
    AircraftState* d_states,    // Device memory
    GPBytecode* d_programs,
    Path* d_paths,
    int n_aircraft,
    int n_steps,
    float dt
) {
    dim3 blocks((n_aircraft + 255) / 256);
    dim3 threads(256);

    for (int step = 0; step < n_steps; step++) {
        fdm_step_kernel<<<blocks, threads>>>(
            d_states, d_programs, d_paths, d_wind, dt, step
        );
        cudaDeviceSynchronize();  // Barrier between steps
    }
}
```

#### Memory Layout (Structure of Arrays for coalescing)

```cuda
// Instead of:
struct AircraftState { float pos[3]; float vel[3]; float quat[4]; ... };
AircraftState states[N];  // Bad: strided access

// Use:
struct AircraftStateSOA {
    float* pos_x;  // [N]
    float* pos_y;  // [N]
    float* pos_z;  // [N]
    float* vel_x;  // [N]
    // ... etc
};  // Good: coalesced access
```

#### Expected Performance

| Platform | Aircraft/Batch | Timesteps | Time | Sims/sec |
|----------|----------------|-----------|------|----------|
| RTX 3080 | 10,000 | 1,000 | ~100ms | 100,000 |
| A100 | 50,000 | 1,000 | ~100ms | 500,000 |
| M1 Max (Metal) | 5,000 | 1,000 | ~200ms | 25,000 |

### Phase 4: Apple Silicon / Metal Path

For M-series chips, use Metal Compute instead of CUDA:

```swift
// Metal compute shader (similar structure to CUDA)
kernel void fdm_step(
    device AircraftState* states [[buffer(0)]],
    device const GPBytecode* programs [[buffer(1)]],
    device const Path* paths [[buffer(2)]],
    constant WindConfig& wind [[buffer(3)]],
    constant float& dt [[buffer(4)]],
    uint i [[thread_position_in_grid]]
) {
    // Same logic as CUDA kernel
}
```

**Advantages of Metal:**
- Unified memory: no explicit CPU↔GPU copies
- Better for smaller batches (lower launch overhead)
- Integrated GPU shares L3 cache with CPU

**Disadvantages:**
- Smaller market (Apple only)
- Less mature tooling than CUDA

### Phase 5: Architecture for 100K+ sims/sec

At extreme scale, the architecture inverts:

**Current (CPU-centric):**
```
autoc (coordinator)
├── Worker 1 (crrcsim) → 150 sims/sec
├── Worker 2 (crrcsim) → 150 sims/sec
├── ...
└── Worker 8 (crrcsim) → 150 sims/sec
Total: ~1200 sims/sec
```

**Future (GPU-centric):**
```
autoc (coordinator)
├── GPU Worker 1 → 50,000 sims/sec
│   └── CUDA: 10,000 aircraft × 5 scenarios
└── GPU Worker 2 (optional) → 50,000 sims/sec
Total: ~100,000 sims/sec
```

**Key Differences:**
- Fewer, larger workers (1-4 GPUs instead of 8-64 CPUs)
- Batch size = generation size (1000+ evals per batch)
- autoc becomes thin coordinator, not compute node
- Evolution still on CPU (parallelizes poorly)

---

## Profiling Commands Reference

**See [PROFILING.md](PROFILING.md) for complete profiling guide with flame graph generation.**

### Quick Reference

```bash
# Build with profiling (both projects have ENABLE_PROFILING option)
cmake -DENABLE_PROFILING=ON ../autoc  # or .. for crrcsim

# Record with perf
perf record -g -F 99 --call-graph dwarf -o profile.data ./autoc

# Generate flame graph
perf script -i profile.data | stackcollapse-perf.pl | flamegraph.pl > flame.svg
```

### Memory bandwidth estimation
```bash
# Estimate memory bandwidth requirement
# 10K aircraft × 1KB state × 1000 steps × 2 (read+write) = 20 GB/step
# At 10K steps/sec = 200 GB/s (within A100 bandwidth: 1.5 TB/s)
```

---

## Risk Assessment

| Risk | Impact | Mitigation |
|------|--------|------------|
| Amdahl's law limit | Can't reach 100K | Profile first, reduce serial fraction |
| GPU memory limits | Batch size capped | Use SOA layout, compress state |
| Numerical divergence | Results differ GPU/CPU | Use deterministic math, compare carefully |
| Dryden RNG on GPU | Non-deterministic | Pre-generate wind sequences on CPU |
| Development time | Months of work | Phase incrementally, measure ROI |

---

## Recommended Next Steps (UPDATED Jan 2026)

**Profiling complete.** Priorities revised based on actual measurements.

1. **Immediate: Disable Gear for Gliders** (1-2 days)
   - Add config knob to crrcsim
   - **+14% performance with minimal risk**
   - See Implementation Proposal Phase 0

2. **Week 1: SIMD Matrix33** (high impact, medium effort)
   - Vectorize hot Matrix33 operations
   - **+8-10% performance**
   - See Implementation Proposal Phase 1

3. **Week 2-3: Wind Model Optimization**
   - Dryden turbulence caching
   - Thermal model simplification
   - **+10-15% performance**
   - See Implementation Proposal Phase 2

4. **Week 4-5: Platform Plugin Architecture**
   - Abstraction for future GPU backend
   - No performance change, but enables Phase 5+

5. **Month 2+: GPU Prototype**
   - Port wind model to CUDA (primary target, not aero)
   - GP interpreter on GPU
   - Target: 10,000+ sims/sec

5. **Month 4+: Optimization**
   - Memory layout optimization (SOA)
   - Multi-GPU support
   - Target: 50,000-100,000 sims/sec

---

## IMPLEMENTATION PROPOSAL: Incremental Working-to-Working (Jan 2026)

Based on profiling results, here is a phased approach where each step produces a working, testable system.

### Phase 0: Quick Win - Disable Gear Physics (1-2 days)

**Goal:** +14% performance by disabling unnecessary wheel physics for aircraft without landing gear.

**Changes:**
1. Add config option `gear.enabled` to crrcsim aircraft XML (default: true)
2. Add runtime flag in `CRRC_AirplaneSim_Larcsim` to skip `gear()` call
3. Set `gear.enabled=false` for glider/flying wing models

**Files:**
- `crrcsim/src/mod_fdm/fdm_larcsim.cpp` - Add conditional around gear() call
- `crrcsim/src/mod_fdm/fdm_larcsim.h` - Add `bool gearEnabled_` member
- Aircraft XML files - Add `<gear enabled="false"/>`

**Verification:**
- Run same GP with gear enabled/disabled
- Verify fitness unchanged (aircraft shouldn't be touching ground during flight)
- Measure performance improvement (~14% expected)

---

### Phase 1: SIMD Matrix33 Operations (1 week)

**Goal:** Vectorize the hot Matrix33 operations that dominate wind calculations (~8-10%).

**Approach:** Platform-specific SIMD backend with compile-time selection.

**Files to create:**
```
crrcsim/src/include/
  crrc_math_simd.h       # Unified interface
  crrc_math_scalar.h     # Fallback (current implementation)
  crrc_math_sse.h        # x86-64 SSE4.1/AVX
  crrc_math_neon.h       # ARM64 NEON (future)
```

**Key operations to vectorize:**
1. `Matrix33::operator*(Matrix33)` - 8% of total time
2. `Matrix33::operator*(Vector3)` - 2% of total time
3. `Matrix33::multrans()` - 3% of total time

**Implementation:**
```cpp
// crrc_math_simd.h - unified interface
#if defined(__AVX__)
  #include "crrc_math_avx.h"
#elif defined(__SSE4_1__)
  #include "crrc_math_sse.h"
#elif defined(__ARM_NEON)
  #include "crrc_math_neon.h"
#else
  #include "crrc_math_scalar.h"
#endif

// crrc_math_sse.h - example
inline Matrix33 Matrix33::operator*(const Matrix33& rhs) const {
    // Use _mm_* intrinsics for 4-wide SIMD
    // Pack 3 rows into 4-element vectors (pad with 0)
    __m128 row0 = _mm_loadu_ps(&m[0][0]);  // loads m[0][0..2] + garbage
    // ... vectorized multiply-add ...
}
```

**Verification:**
- Unit tests comparing scalar vs SIMD results (bit-exact or within epsilon)
- Benchmark Matrix33 operations before/after
- Full simulation regression test

**Expected improvement:** 8-10% overall (Matrix33 ops are ~10% of runtime)

---

### Phase 2: Wind Model Optimization (1-2 weeks)

**Goal:** Reduce wind/turbulence overhead from 27.5% to ~15%.

**Sub-phases:**

#### 2a: Dryden Turbulence Caching
- Pre-compute turbulence tables at simulation start
- Interpolate during timesteps instead of full Dryden calculation
- Trade memory for compute

#### 2b: Thermal Model Simplification
- `Thermal::update()` at 7% - profile for specific hotspots
- Consider reduced-fidelity mode for GP training (full fidelity for elite re-eval)

#### 2c: Wind Gradient Batching
- `CalculateWindGrad()` computes gradient via finite differences
- Cache gradient when aircraft position hasn't moved significantly

**Files:**
- `crrcsim/src/mod_windfield/windfield.cpp`
- `crrcsim/src/mod_windfield/thermalmodel.cpp`

**Verification:**
- Compare wind vectors before/after optimization (should be identical or within tolerance)
- Verify GP training still produces good controllers
- Measure performance improvement

---

### Phase 3: Platform Plugin Architecture (2 weeks)

**Goal:** Abstraction layer for physics backends (CPU scalar, CPU SIMD, GPU future).

**Design:**
```cpp
// physics_backend.h
class PhysicsBackend {
public:
    virtual ~PhysicsBackend() = default;

    // Simulate N aircraft for one timestep
    virtual void stepBatch(
        AircraftState* states,
        const GPBytecode* programs,
        const Path* paths,
        const WindConfig& wind,
        float dt,
        int n_aircraft
    ) = 0;

    // Factory
    static std::unique_ptr<PhysicsBackend> create(const std::string& type);
};

// Implementations:
class ScalarPhysicsBackend : public PhysicsBackend { ... };  // Current code
class SimdPhysicsBackend : public PhysicsBackend { ... };    // Phase 1 SIMD
class CudaPhysicsBackend : public PhysicsBackend { ... };    // Future GPU
```

**Benefits:**
- Clean separation of concerns
- Easy A/B testing of implementations
- GPU backend can be added incrementally
- Runtime selection via config

**Files:**
```
crrcsim/src/mod_fdm/
  physics_backend.h
  physics_backend_scalar.cpp
  physics_backend_simd.cpp
  physics_backend_cuda.cpp  # Future
```

---

### Phase 4: GPU Port Preparation (1 month)

**Goal:** Refactor crrcsim physics for GPU-friendly data layout.

**Sub-phases:**

#### 4a: Structure of Arrays (SOA) Conversion
- Convert `AircraftState` from AOS to SOA layout
- Enables coalesced memory access on GPU
- Keep scalar wrapper for compatibility

#### 4b: Batch Simulation API
- Single entry point: `simulateBatch(n_aircraft, n_steps)`
- Returns fitness values, not full state traces
- Worker-side fitness computation

#### 4c: GPU Kernel Prototypes
- Port wind calculation to CUDA kernel
- Port ls_step() integration to CUDA
- Benchmark with 1000+ aircraft

**This phase sets up for full GPU port but each sub-phase is independently testable.**

---

### Verification at Each Phase

| Phase | Test | Pass Criteria |
|-------|------|---------------|
| 0 | Same GP, gear on vs off | Fitness identical (within 0.1%) |
| 1 | Matrix33 unit tests | Bit-exact or <1e-6 epsilon |
| 1 | Full simulation | Fitness identical |
| 2 | Wind vector comparison | <1% difference |
| 2 | GP training | Converges to similar fitness |
| 3 | Backend A/B test | All backends produce same results |
| 4 | SOA vs AOS comparison | Bit-exact fitness |

---

### Timeline Summary

| Phase | Duration | Expected Speedup | Cumulative |
|-------|----------|------------------|------------|
| 0: Disable gear | 1-2 days | +14% | 1.14x |
| 1: SIMD Matrix33 | 1 week | +8-10% | 1.25x |
| 2: Wind optimization | 1-2 weeks | +10-15% | 1.40x |
| 3: Plugin architecture | 2 weeks | - (prep) | 1.40x |
| 4: GPU preparation | 1 month | - (prep) | 1.40x |
| 5: GPU port (future) | 2+ months | 10-50x | 15-70x |

**Conservative target after Phase 2:** 1200 × 1.4 = **1680 sims/sec**
**Target after GPU port:** 10,000-60,000 sims/sec

---

## ARCHITECTURAL CONSIDERATION: Scenario Caching & EvalTask Redesign

### Current Problem

Today's EvalTask sends paths with every evaluation:
```
EvalData = GP bytecode + pathList + scenarioList
```

With 36 scenarios × N paths × ~36 bytes/pathpoint × ~1000 points = **1-5 MB per eval**

Future scale: 1000s of scenarios (craft variants × wind × etc.) = **10-50 MB per eval**

This is wasteful because paths/scenarios are constant within a generation (or even across generations).

### Proposed Architecture: Scenario Registry

**Separate the "what to simulate" from "how to simulate":**

```
┌─────────────────────────────────────────────────────────────────────┐
│                        autoc (coordinator)                          │
├─────────────────────────────────────────────────────────────────────┤
│  ScenarioRegistry                                                   │
│  ├── pathSets: map<uint64_t hash, vector<vector<Path>>>            │
│  ├── scenarioSets: map<uint64_t hash, vector<ScenarioMetadata>>    │
│  └── craftVariants: map<uint64_t hash, CraftParams>  (future)      │
├─────────────────────────────────────────────────────────────────────┤
│  EvalTask (NEW - lightweight)                                       │
│  ├── gpBytecode: vector<char>                                      │
│  ├── gpHash: uint64_t                                              │
│  ├── pathSetHash: uint64_t      ← Reference, not data              │
│  ├── scenarioSetHash: uint64_t  ← Reference, not data              │
│  └── craftVariantHash: uint64_t ← Future                           │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     Worker (crrcsim or GPU)                         │
├─────────────────────────────────────────────────────────────────────┤
│  LocalScenarioCache                                                 │
│  ├── pathSets: LRU<hash, vector<vector<Path>>>                     │
│  ├── scenarioSets: LRU<hash, vector<ScenarioMetadata>>             │
│  └── craftVariants: LRU<hash, CraftParams>                         │
├─────────────────────────────────────────────────────────────────────┤
│  On EvalTask receive:                                               │
│  1. Check cache for pathSetHash                                    │
│  2. If miss: request from autoc (one-time per generation)          │
│  3. Simulate using cached data                                     │
└─────────────────────────────────────────────────────────────────────┘
```

### Wire Protocol Changes

**New message types:**

```cpp
// Lightweight eval request (< 2KB typical)
struct EvalTaskCompact {
    uint64_t gpHash;
    vector<char> gpBytecode;      // ~200-2000 bytes
    uint64_t pathSetHash;         // 8 bytes (reference)
    uint64_t scenarioSetHash;     // 8 bytes (reference)
    uint64_t craftVariantHash;    // 8 bytes (future)
    bool isEliteReeval;
};

// Cache request (worker → autoc, rare)
struct CacheMissRequest {
    enum Type { PathSet, ScenarioSet, CraftVariant };
    Type type;
    uint64_t hash;
};

// Cache response (autoc → worker, large but rare)
struct CacheResponse {
    uint64_t hash;
    // One of:
    vector<vector<Path>> pathSet;           // ~1-5 MB
    vector<ScenarioMetadata> scenarioSet;   // ~1-10 KB
    CraftParams craftVariant;               // ~1 KB (future)
};
```

**Savings:**
- Current: 1-5 MB per eval × 1000 evals/gen = 1-5 GB/gen
- Proposed: 2 KB per eval × 1000 evals/gen + 5 MB cache = **~7 MB/gen** (99% reduction)

### Sharding for Large Scenario Sets

When scenarios grow to 1000s:

**Option A: Scenario Groups (Recommended)**
```
ScenarioRegistry {
    // Group by similarity for cache locality
    windVariants: [calm, light, moderate, gusty]  // 4 groups
    craftVariants: [baseline, heavy, light, draggy]  // 4 groups
    pathVariants: [figure8, racetrack, landing, aerobatic]  // 4 groups

    // Total: 4 × 4 × 4 = 64 combinations, but each is a small set
    // Worker caches only the groups it's currently evaluating
}
```

**Option B: Hierarchical Hashing**
```
ScenarioSetHash = hash(windGroupHash, craftGroupHash, pathGroupHash)

Worker caches at group level, combines on-demand:
- windGroup[hash1] + craftGroup[hash2] + pathGroup[hash3]
```

**Option C: Streaming Scenarios**
```
// For very large sets, stream scenarios to workers
EvalBatch {
    gpPrograms: [...]
    scenarioStream: AsyncIterator<Scenario>  // Lazy loading
}
```

### The "GPU-in-Process" Architecture

At extreme scale, the IPC overhead becomes significant. Consider embedding the GPU simulation directly in autoc:

```
┌─────────────────────────────────────────────────────────────────────┐
│                     autoc (single process)                          │
├─────────────────────────────────────────────────────────────────────┤
│  GP Evolution (CPU)                                                │
│  ├── Population management                                         │
│  ├── Selection, crossover, mutation                                │
│  └── Fitness assignment                                            │
├─────────────────────────────────────────────────────────────────────┤
│  GPU Simulation Module (in-process)                                │
│  ├── CUDA context (persistent)                                     │
│  ├── Device memory pools (pre-allocated)                           │
│  │   ├── pathSets: 100 MB reserved                                │
│  │   ├── scenarioSets: 10 MB reserved                             │
│  │   └── aircraftStates: 1 GB reserved (10K aircraft × 1KB × batch)│
│  ├── simulateBatch(gpPrograms[], pathSetHash, scenarioSetHash)     │
│  └── Returns: fitness[] (on CPU)                                   │
└─────────────────────────────────────────────────────────────────────┘
```

**Advantages:**
- Zero IPC overhead (no serialize/deserialize)
- Zero network latency
- Shared memory for paths/scenarios (no copy needed)
- Simpler deployment (single binary)
- Perfect for single-GPU workstations

**Disadvantages:**
- GPU memory shared with autoc process
- Can't distribute across machines easily
- Harder to debug (GPU + CPU in same process)
- Requires CUDA at compile time

**Hybrid Approach:**
```cpp
// Runtime selection based on config
std::unique_ptr<SimulationBackend> backend;
if (config.gpuInProcess && cuda_available()) {
    backend = std::make_unique<CudaInProcessBackend>();
} else if (config.workers > 0) {
    backend = std::make_unique<DistributedWorkerBackend>(config.workers);
} else {
    backend = std::make_unique<SingleThreadCpuBackend>();  // Debugging
}

// All backends implement same interface
vector<float> fitness = backend->evaluateBatch(gpPrograms, pathSetHash, scenarioSetHash);
```

### Implementation Phases for Caching

**Phase A: Add hashing to current protocol (1 week)**
1. Compute hash of pathList and scenarioList in autoc
2. Send hash alongside data (no behavior change)
3. Log cache hit potential (how often hashes repeat)

**Phase B: Worker-side caching (1 week)**
1. Worker maintains LRU cache of pathSets by hash
2. autoc sends `includesPathSet: bool` flag
3. On cache hit, worker skips pathSet deserialization
4. Verify identical results

**Phase C: Lazy loading protocol (2 weeks)**
1. EvalTask sends only hashes (compact)
2. Worker requests missing data on cache miss
3. autoc responds with CacheResponse
4. Benchmark reduction in bandwidth

**Phase D: GPU memory pools (with GPU port)**
1. Pre-allocate device memory for cached data
2. Upload pathSets to GPU once per generation
3. EvalBatch references GPU-resident data by hash

### Cache Sizing Estimates

| Scenario Scale | PathSets | ScenarioSets | Total Cache | Notes |
|----------------|----------|--------------|-------------|-------|
| Current (36) | 5 MB | 1 KB | ~5 MB | Fits in L3 |
| Medium (256) | 40 MB | 10 KB | ~40 MB | Fits in RAM |
| Large (1024) | 150 MB | 40 KB | ~150 MB | Needs sharding |
| Extreme (10K) | 1.5 GB | 400 KB | ~1.5 GB | GPU memory |

### Impact on Implementation Phases

Update Phase 3 (Platform Plugin Architecture) to include:

```cpp
class SimulationBackend {
public:
    virtual ~SimulationBackend() = default;

    // Cache management
    virtual void uploadPathSet(uint64_t hash, const vector<vector<Path>>& paths) = 0;
    virtual void uploadScenarioSet(uint64_t hash, const vector<ScenarioMetadata>& scenarios) = 0;
    virtual bool hasPathSet(uint64_t hash) const = 0;
    virtual bool hasScenarioSet(uint64_t hash) const = 0;

    // Batch evaluation (uses cached data by hash)
    virtual vector<EvalResult> evaluateBatch(
        const vector<GPBytecode>& programs,
        uint64_t pathSetHash,
        uint64_t scenarioSetHash
    ) = 0;
};
```

This interface supports:
- CPU workers with local cache
- GPU workers with device memory cache
- In-process GPU with shared memory

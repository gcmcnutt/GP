# [TRANSFERRED to specs/011-gpu-native] SCALEUP: Scaling AutoC for GPU-Accelerated Simulation

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

**Profiling complete. Phase 0 and Phase 1 complete.** Priorities revised based on actual measurements.

1. ~~**Immediate: Disable Gear for Gliders**~~ **DONE** (Jan 2026)
   - Added `<wheels enabled="0">` to aircraft XML
   - **+14% achieved**

2. ~~**Immediate: Arena Thermal System**~~ **DONE** (Jan 2026)
   - Replaced global 12.8km thermal grid with bounded arena system
   - **+13% achieved** (20x more efficient thermal simulation)

3. ~~**Next: SIMD Matrix33**~~ **DONE** (Jan 2026)
   - Migrated Matrix33/Vector3 to Eigen library for SIMD vectorization
   - **+15.4% achieved** (exceeded 8-10% estimate)
   - SSE2 vectorization confirmed via perf profile

4. **Next: Wind Model Optimization**
   - Dryden turbulence caching (still ~25% of runtime post-Eigen)
   - **+10-15% performance expected**
   - See Implementation Proposal Phase 2

5. **Later: Platform Plugin Architecture**
   - Abstraction for future GPU backend
   - No performance change, but enables GPU port

6. **Future: GPU Prototype**
   - Port wind model to CUDA (primary target, not aero)
   - GP interpreter on GPU
   - Target: 10,000+ sims/sec

7. **Future: GPU Optimization**
   - Memory layout optimization (SOA)
   - Multi-GPU support
   - Target: 50,000-100,000 sims/sec

---

## IMPLEMENTATION PROPOSAL: Incremental Working-to-Working (Jan 2026)

Based on profiling results, here is a phased approach where each step produces a working, testable system.

### Phase 0: Quick Wins - Gear Disable + Arena Thermals (COMPLETED Jan 2026)

**Goal:** Eliminate unnecessary physics calculations for GP training scenarios.

#### Phase 0a: Disable Gear Physics (COMPLETED)

**Result:** Gear physics disabled for gliders via `<wheels enabled="0">` in aircraft XML.

**Impact:** 13.7% → 0% (gear no longer appears in profile)

#### Phase 0b: Arena-Bounded Thermal System (COMPLETED)

**Problem:** Global thermal grid simulated ~393 thermals across 12.8km grid, but training arena is only ±75m.
The old system spent 7% + 6.5% = 13.5% of runtime on `update_thermals()` and `Thermal::update()`
processing thermals that had zero effect on aircraft in the training arena.

**Solution:** New `arena_thermal` module that:
- Spawns 0-5 thermals within bounded arena (±150m spawn, ±75m active)
- Uses simplified bell-curve profile (vs complex v3 shell model)
- Fixed array of max 8 thermals (no linked list traversal)
- Bypasses global thermal grid entirely when enabled

**Files created:**
- `crrcsim/src/mod_windfield/arena_thermal.h` - Header with `ArenaThermalConfig`, `ArenaThermal`, `ArenaThermalField`
- `crrcsim/src/mod_windfield/arena_thermal.cpp` - Implementation
- Modified `windfield.cpp` to integrate arena thermals

**Config (in autoc_config.xml location section):**
```xml
<arena_thermals enabled="1">
  <bounds x_min="-150" x_max="150" y_min="-150" y_max="150" />
  <count_min>0</count_min>
  <count_max>5</count_max>
  <strength_mean>2.0</strength_mean>
  <strength_sigma>0.5</strength_sigma>
  <radius_mean>20</radius_mean>
  <radius_sigma>5</radius_sigma>
  <lifetime_mean>180</lifetime_mean>
  <lifetime_sigma>60</lifetime_sigma>
  <height_m>300</height_m>
</arena_thermals>
```

**Performance comparison (perf.out vs perf2.out):**

| Component | Before | After | Change |
|-----------|--------|-------|--------|
| `gear()` / `WheelSystem::update()` | 13.66% | 0% | **-13.66%** |
| `update_thermals()` | 7.02% | 0% | **-7.02%** |
| `Thermal::update()` | 6.56% | 0% | **-6.56%** |
| `calculate_wind()` | ~5.7% | 2.83% | **-2.9%** |
| `ArenaThermalField::calculateThermalWind()` | N/A | 0.65% | +0.65% |
| **Total reduction** | | | **~30%** |

**Arena thermal efficiency:** 0.65% vs 13.6% for global grid = **20x improvement**

---

### Phase 1: Eigen Integration for Matrix/Vector SIMD (1 week)

**Goal:** Replace hand-rolled `Matrix33`/`Vector3` with Eigen for automatic SIMD vectorization.

**Why Eigen:**
- Already a dependency in autoc (`gp_types.h`) and xiao-gp (`ArduinoEigen`)
- Header-only, zero runtime overhead
- Auto-detects SSE/AVX (x86) and NEON (ARM) at compile time
- Expression templates eliminate temporaries and fuse operations
- Battle-tested, well-maintained library

**Current state:**
- crrcsim uses custom `CRRCMath::Matrix33` and `CRRCMath::Vector3` (circa 2005)
- Simple scalar loops, no SIMD
- `Matrix33::operator*` is 11.75% of runtime (in `calculate_gust`)
- Additional 3% in `aero()` for total ~15%

**Approach: Option 1 - Wrapper with Eigen backend**

Keep existing API stable, replace internal storage with Eigen:

```cpp
// matrix33.h - before
class Matrix33 {
  double v[3][3];  // Row-major storage
public:
  Matrix33 operator*(const Matrix33& b) const;  // Scalar loops
  // ...
};

// matrix33.h - after
#include <Eigen/Dense>

class Matrix33 {
  // Row-major to match existing v[row][col] access pattern
  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> m_;
public:
  // Compatibility accessor for direct v[i][j] access (~100 call sites)
  double& operator()(int row, int col) { return m_(row, col); }
  double operator()(int row, int col) const { return m_(row, col); }

  // Hot path - now SIMD vectorized
  Matrix33 operator*(const Matrix33& b) const {
    Matrix33 result;
    result.m_ = m_ * b.m_;
    return result;
  }

  Vector3 operator*(const Vector3& b) const;  // Similar treatment
  // ...
};
```

**Migration steps:**

1. **Add Eigen dependency to crrcsim CMakeLists.txt**
   ```cmake
   find_package(Eigen3 3.3 REQUIRED NO_MODULE)
   target_link_libraries(crrcsim Eigen3::Eigen)
   ```

2. **Update Vector3 class** (`src/mod_math/vector3.h`)
   - Replace `double r[3]` with `Eigen::Vector3d`
   - Keep `r[i]` accessor for ~509 call sites (via `operator[]` or inline accessor)

3. **Update Matrix33 class** (`src/mod_math/matrix33.h`, `matrix33.cpp`)
   - Replace `double v[3][3]` with `Eigen::Matrix<double, 3, 3, Eigen::RowMajor>`
   - Update ~100 direct `v[i][j]` accesses to `operator()(i,j)`
   - Delete scalar loop implementations (Eigen inlines everything)

4. **Handle Quaternion** (`src/mod_math/quaternion.cpp`)
   - ~44 lines access `mat.v[i][j]` for quat↔matrix conversion
   - Update to `mat(i,j)` syntax

5. **Update other files with direct access**
   - `windfield.cpp` - gradient calculations (~15 sites)
   - `fdm_larcsim.cpp` - atmospheric rotation (~3 sites)
   - `fdm_002.cpp` - similar (~12 sites)
   - `eom01.cpp` - LocalToBody matrix (~20 sites)

**Key files to modify:**
```
crrcsim/CMakeLists.txt                    # Add Eigen dependency
crrcsim/src/mod_math/vector3.h            # Eigen::Vector3d backend
crrcsim/src/mod_math/matrix33.h           # Eigen::Matrix3d backend
crrcsim/src/mod_math/matrix33.cpp         # Delete (all inlined)
crrcsim/src/mod_math/quaternion.cpp       # v[i][j] → (i,j) syntax
crrcsim/src/mod_fdm/eom01/eom01.cpp       # v[i][j] → (i,j) syntax
crrcsim/src/mod_windfield/windfield.cpp   # v[i][j] → (i,j) syntax
crrcsim/src/mod_fdm/fdm_002/fdm_002.cpp   # v[i][j] → (i,j) syntax
crrcsim/src/mod_fdm/fdm_larcsim/fdm_larcsim.cpp  # v[i][j] → (i,j)
```

**Storage layout consideration:**
- Eigen default is column-major, crrcsim uses row-major `v[row][col]`
- Use `Eigen::RowMajor` template parameter to match existing semantics
- Eigen still vectorizes row-major matrices

**Verification:**
1. **Unit tests** - Compare old vs new Matrix33/Vector3 operations (bit-exact)
2. **Simulation regression** - Run GP training, verify fitness convergence unchanged
3. **Perf profile** - Measure reduction in `calculate_gust` and `aero` overhead

**Expected improvement:** 8-10% overall (Matrix33 ops currently ~15% of runtime)

---

#### Phase 1 Results: COMPLETED (Jan 2026)

**Approach taken:** Full Eigen migration (Option 2 - typedef to Eigen types)

```cpp
// vector3.h
using Vector3 = Eigen::Vector3d;

// matrix33.h
using Matrix33 = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
```

**Benchmark methodology:**
- Controlled A/B test with identical seeds (GPSeed=42, WindSeedBase=12121, RandomPathSeedB=13337)
- PERFORMANCE_BUILD flags (`-O3 -march=native -flto -ffast-math`)
- 21 generations, 1000 population, 9 wind scenarios = 189,180 total simulations
- Intel 4-core CPU with hyperthreading

**Results:**

| Configuration | Threads | Total Time | Sims/sec | vs Baseline |
|---------------|---------|------------|----------|-------------|
| Without Eigen | 8 | 788.3 sec | 240 | baseline |
| With Eigen | 8 | 683.2 sec | 277 | **+15.4%** |
| With Eigen | 4 | 805.7 sec | 235 | -2% vs baseline |

**Key findings:**
1. **Eigen SIMD delivers +15.4% speedup** - exceeds the 8-10% estimate
2. **Thread scaling:** 8 threads (2× cores) still optimal; hyperthreading provides ~18% boost over 4 threads
3. **Numerical equivalence:** Fitness values identical to 6+ decimal places (minor FP differences in later digits)
4. **Rule of thumb:** Use `threads = cores × 2` on SMT systems, `threads = cores` otherwise

**SIMD evidence from perf profile:**
- `double __vector(2)` operations visible (SSE2 packed doubles)
- Eigen internal ops: `pmul`, `predux`, matrix products
- Function signatures show `Eigen::Matrix<double, 3, 3, 1, 3, 3>` (RowMajor)

**Files modified:**
- `src/mod_math/vector3.h` - typedef to `Eigen::Vector3d`
- `src/mod_math/matrix33.h` - typedef to `Eigen::Matrix<double,3,3,RowMajor>` + helper functions
- `src/mod_math/CMakeLists.txt` - removed vector3.cpp/matrix33.cpp, added Eigen dependency
- ~200 call sites updated from `.r[i]`/`.v[i][j]` to `(i)`/`(i,j)` accessor syntax

---

#### Future mod_math Eigen Cleanup (Low Priority)

Analysis of remaining mod_math files for potential Eigen migration:

| File | Eigen Replacement | Effort | Benefit | Recommendation |
|------|------------------|--------|---------|----------------|
| **matrix44.h** | `Eigen::Matrix4d` | Low | Medium (not in hot path) | Could migrate |
| **quaternion.h/cpp** | Partial - `Eigen::Quaterniond` for storage | Medium | Low | Keep custom integration logic |
| **intgr.h** | None (ODE integration, not linear algebra) | N/A | N/A | Keep as-is |
| **pt1.h** | None (signal processing filter) | N/A | N/A | Keep as-is |
| **linearreg.h** | None (statistics) | N/A | N/A | Keep as-is |
| **ratelim.h** | None (value clamping) | N/A | N/A | Keep as-is |

**Notes on quaternion.h:**
- Eigen has `Eigen::Quaterniond` with `toRotationMatrix()`, normalization, slerp
- However, crrcsim's Quaternion classes include custom integration state (`IntegrationsverfahrenB<double> e0,e1,e2,e3`)
- Multiple implementations (Quaternion_001, _002, _003) with different integration schemes
- Would require keeping the integration wrapper around Eigen's quaternion storage

**Verdict:** Matrix44 is the cleanest candidate for future cleanup. Quaternion partial migration possible but not high value. Other files are not linear algebra and should stay as-is.

---

**Future Option 2:** ~~Full Eigen migration (remove wrapper classes entirely)~~ **DONE**
- `using Matrix33 = Eigen::Matrix3d;`
- `using Vector3 = Eigen::Vector3d;`
- ~~More invasive, but cleaner long-term~~
- ~~Can migrate incrementally after Option 1 is validated~~

---

### Phase 2: Dryden Turbulence Optimization (1-2 weeks)

**Goal:** Reduce Dryden turbulence overhead (currently ~25% of runtime after Phase 0/1).

**Current profile (perf2.out):**
- `calculate_gust()`: 25.74% total
  - `Matrix33::operator*`: 11.75% → addressed by Phase 1 Eigen
  - `RandGauss`: 5.36% → RNG calls per timestep
  - `pow`: 2.25% → power calculations

**Note:** Thermal optimization is DONE (Phase 0b achieved 20x improvement via arena thermals).

**Sub-phases:**

#### 2a: Dryden State Caching
The Dryden model computes turbulence as a filtered noise process. Some state can be cached:
- Filter coefficients depend only on airspeed and altitude bands
- Pre-compute coefficient tables for discrete altitude/airspeed bins
- Interpolate coefficients instead of recomputing

#### 2b: RNG Optimization
`RandGauss` at 5.36% is significant:
- Consider faster Gaussian approximation (Ziggurat method)
- Or pre-generate Gaussian sequences per scenario (determinism benefit too)
- Trade memory for compute

#### 2c: Power Function Optimization
`pow()` at 2.25%:
- Identify which powers are constant (e.g., `pow(x, 0.5)` → `sqrt(x)`)
- Use fast approximations for non-critical paths
- Consider lookup tables for common exponents

**Files:**
- `crrcsim/src/mod_windfield/windfield.cpp` - Dryden implementation
- `crrcsim/src/mod_misc/crrc_rand.cpp` - RNG functions

**Verification:**
- Compare wind vectors before/after optimization (should be identical or within tolerance)
- Verify GP training still produces good controllers
- Measure performance improvement

**Expected improvement:** 10-15% overall (targeting 25% → ~15% for Dryden)

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
| 0a | Same GP, gear on vs off | Fitness identical (within 0.1%) |
| 0b | Arena thermals vs global grid | Wind vectors equivalent, 20x faster |
| 1 | Eigen vs old Matrix33/Vector3 | Bit-exact or <1e-10 epsilon |
| 1 | Full simulation | Fitness identical to pre-Eigen |
| 2 | Dryden output comparison | <1% difference |
| 2 | GP training | Converges to similar fitness |
| 3 | Backend A/B test | All backends produce same results |
| 4 | SOA vs AOS comparison | Bit-exact fitness |

---

### Timeline Summary

| Phase | Duration | Expected Speedup | Actual Speedup | Cumulative | Status |
|-------|----------|------------------|----------------|------------|--------|
| 0a: Disable gear | 1-2 days | +14% | +14% | 1.14x | **DONE** |
| 0b: Arena thermals | 2-3 days | +13% | +13% | 1.30x | **DONE** |
| 1: Eigen integration | 1 week | +8-10% | **+15.4%** | 1.50x | **DONE** |
| 2: Dryden optimization | 1-2 weeks | +10-15% | - | ~1.65x | Pending |
| 3: Plugin architecture | 2 weeks | - (prep) | - | ~1.65x | Pending |
| 4: GPU preparation | 1 month | - (prep) | - | ~1.65x | Pending |
| 5: GPU port (future) | 2+ months | 10-50x | - | 15-75x | Future |

**Achieved after Phase 0:** ~30% reduction in crrcsim execution time
**Achieved after Phase 1:** ~50% total improvement (Eigen exceeded expectations)
**Conservative target after Phase 2:** 1200 × 1.65 = **~2000 sims/sec**
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

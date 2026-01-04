# Attention to Detail: Fitness Calculation NaN Bug Resolution

## Summary of Bug Found and Fixed (2026-01-02)

### The Problem
GP evolution showed negative fitness values (down to -8.08746 by gen 41), which violated the design invariant that fitness >= 0 with 0 being perfect. Additionally, tiny NaN warnings appeared occasionally due to floating-point precision at dot product boundaries.

### Root Cause Analysis

**Primary Issue: Throttle Energy Calculation Bug**
- **Line**: autoc.cc:773-776 (and duplicate in determinism test section)
- **Bug**: Throttle command range is [-1, 1] but fitness calculation treated it as [0, 1]
  ```cpp
  // WRONG (old code):
  gp_scalar throttlePercent = current_throttle * 100.0f;
  // If throttle = -0.5: throttlePercent = -50.0
  // throttle_energy_sum += pow(-50.0, 1.0) = -50.0  // NEGATIVE!
  ```
- **Evidence**:
  - `aircraft_state.h:252` clamps throttle to [-1.0, 1.0]
  - `inputdev_autoc.cpp:764` maps to crrcsim: `inputs->throttle = throttleCommand / 2.0 + 0.5`
  - This maps [-1,1] → [0,1] where -1 = motor off (zero power), +1 = full power
- **Impact**: GP individuals using negative throttle (attempting low-power gliding) got negative fitness contributions, causing total fitness to go negative

**Fix Applied**:
```cpp
// CORRECT (new code):
gp_scalar throttleNormalized = (current_throttle + 1.0f) * 0.5f;  // Maps [-1,1] -> [0,1]
gp_scalar throttlePercent = throttleNormalized * 100.0f;          // Range: 0 to 100
throttle_energy_sum += pow(throttlePercent, THROTTLE_EFFICIENCY_WEIGHT);
```

Now:
- Throttle = -1.0 → energy = 0.0 (perfect efficiency - motor off)
- Throttle = 0.0 → energy = 50.0 (half power)
- Throttle = +1.0 → energy = 100.0 (full power)

**Secondary Issue: Dot Product Floating-Point Precision**
- **Line**: autoc.cc:761-763
- **Bug**: Dot product of two normalized vectors should be in [-1, 1], but floating-point rounding can produce 1.00000001
  ```cpp
  direction_alignment = aircraft_movement.normalized().dot(frame.tangent);
  // Could be 1.00000001 due to floating-point error
  movementDirectionError = (1.0 - 1.00000001) * 50.0 = -0.000005 * 50.0 = -0.00025
  // pow(-0.00025, 1.3) = NaN
  ```
- **Fix**: Clamp to valid range
  ```cpp
  direction_alignment = std::clamp(direction_alignment, -1.0f, 1.0f);
  ```
- **Justification**: This is NOT masking a bug - it's handling legitimate floating-point precision limits at mathematical boundaries

### Files Changed

1. **autoc.cc:774-777** (throttle energy calculation - main evaluation)
2. **autoc.cc:1243-1246** (throttle energy calculation - determinism test section)
3. **autoc.cc:763** (dot product clamp - main evaluation)
4. **autoc.cc:1234** (dot product clamp - determinism test section)
5. **autoc.cc:905** (normalization_factor threshold increased from 0.0f to 0.1f to prevent divide-by-near-zero)
6. **autoc.cc:1363** (normalization_factor threshold - determinism test section)

### Testing Results

**Before Fix**:
- Gen 0-6: Fitness ~37 (stuck)
- Gen 7: Breakthrough to 7.62
- Gen 11: First negative fitness -0.18
- Gen 21: -5.47 (new best)
- Gen 41: -8.09 (continued "improvement" into negative territory)
- Occasional errors: `ERROR: movementDirectionError=-5.96046e-06 direction_alignment=1`

**After Fix** (minisim, multi-threaded, non-demetic):
- ✅ Monotonically decreasing fitness (as expected for "lower is better")
- ✅ Fitness never goes negative
- ✅ No NaN errors
- ✅ Multi-threaded determinism working perfectly

## Next Steps: Path to Production-Ready GP Evolution

### Phase 1: CRRCSIM Perfect Determinism (CURRENT PRIORITY)
**Goal**: Achieve same deterministic behavior in crrcsim as we now have in minisim

**Status**:
- Minisim: ✅ Perfect determinism (multi-threaded, non-demetic)
- CRRCsim: ⚠️ Still showing variance after throttle fix (2026-01-02)

**Current Test Configuration** (autoc.ini):
```ini
PopulationSize                  = 1000
NumberOfGenerations             = 20
DemeticGrouping                 = 0        # Non-demetic
WindScenarios                   = 2        # Testing with 2 wind scenarios
WindSeedBase                    = 4337
WindSeedStride                  = 23
EvalThreads                     = 8        # Multi-threaded
MinisimProgram                  = ./crrcsim.sh  # Using crrcsim
```

**Observed Behavior** (from data.stc):
- Gen 0: Best = 189.673843
- Gen 1: Best = 1.9e+02 (appears unchanged, but scientific notation)
- Gen 2: Best = 1.9e+02 (appears unchanged)
- Gen 3: Best = 1.9e+02 (appears unchanged)
- **Observation**: Best fitness appears stable across first few generations, but need to verify if there are small fluctuations hidden by formatting

**Current Variance Root Causes** (from previous investigation):
- Fixed: RandGauss Box-Muller state persistence ✅
- Fixed: eta1-4 Dryden turbulence RNG state ✅
- Fixed: Thermal generation RNG state ✅
- Fixed: Throttle energy calculation (negative fitness bug) ✅
- **Remaining**: Unknown source of variance

**Next Debugging Steps**:
1. **WAITING**: User is testing single-threaded to isolate multi-threading issues
2. **Enable determinism test**: Change `if (false && gen == 2)` to `if (gen == 2)` at autoc.cc:1543
3. **Consider disabling thermals**: User suggested disabling thermal engine in crrcsim to narrow down variance sources
   - Thermals add random updraft positions/strengths even with fixed seed
   - May contribute to variance if not being reset properly
4. **Measure actual variance** with determinism test enabled
5. **Hunt down remaining sources**:
   - Check for uninitialized variables in FDM (EOM01)
   - Look for floating-point ordering differences across threads
   - Investigate SDL_GetTicks() or other time-dependent calls
   - Profile variance patterns: is it consistent across certain paths? Wind scenarios?
   - Consider if wind scenario assignment is truly deterministic across threads
6. **Target**: Variance < 0.0001 units (needed for 4th-order energy optimization to work)

**Debug Facilities**:
- **Determinism Test** (autoc.cc:1543): Clones best individual to entire population at gen 2, re-evaluates, reports variance
  - Currently disabled: `if (false && gen == 2)`
  - To enable: `if (gen == 2)` or `if (true && gen == 2)`
  - Reports: min, max, avg, variance, and shows fitness groupings
- **Single-threaded test**: Set `EvalThreads = 1` to eliminate thread ordering issues
- **Disable thermals**: May need to modify crrcsim scenery config or wind initialization

**Why This Matters**:
- Current 0.005 variance masks energy cost signal (4th-order metric)
- Evolution cannot optimize throttle efficiency until variance is eliminated
- Elite re-evaluation exposes variance, causing non-monotonic fitness

### Phase 2: Demetic Mode Elite Preservation
**Goal**: Fix fitness tracking in demetic mode where elite's aggregated fitness is lost

**Current Status**:
- Non-demetic: ✅ Working (with small 0.001 regressions due to 0.005 crrcsim variance)
- Demetic: ❌ Broken (fitness jumps wildly: 135→191→192→191→174→180→173→177)

**Root Cause**:
- Elite has aggregated fitness across all scenarios
- Regular population members have single-scenario fitness
- When elite is re-evaluated on single scenario, aggregated fitness is lost
- Attempted solution exists (autoc.cc:623-661, 665-744) but partially working

**Actions Required**:
1. **Verify current behavior**: Add logging to confirm elite is being re-evaluated
2. **Option A**: Skip elite in evaluation entirely
   - Mark elite with flag (hasAggregatedFitness)
   - Skip scenario assignment and evalTask() for elite
   - Preserve fitness across generations
3. **Option B**: Restore aggregated fitness after evaluation
   - After GPPopulation::evaluate() returns
   - Scan for individual with hasAggregatedFitness flag
   - Restore bestOfPopulation pointer
4. **Test**: Best fitness should be strictly monotonically improving (no jumps)

**Why Defer Until After Phase 1**:
- If we achieve perfect crrcsim determinism (variance=0), elite re-evaluation won't cause regressions
- Easier to debug demetic issues without variance noise
- Can verify solution works by checking monotonic fitness

### Phase 3: Path Length Normalization
**Goal**: Ensure all paths contribute equally to fitness regardless of duration

**Current Problem**:
- 30-second path has 2.5x more samples than 12-second paths
- Evolution optimizes for long path performance, ignores short paths
- Total error sum biases toward longer paths

**Proposed Solution**:
- Normalize by path duration or expected sample count
- Use per-second error metrics instead of total sum
- Options:
  - `normalized_error = total_error / path_duration_seconds`
  - `normalized_error = total_error / expected_sample_count`
  - Weight by `1/duration` to equalize importance

**Why This Matters**:
- Real-world flights vary in duration
- Controller should perform well on all path types equally
- Current bias toward long paths may produce controllers that fail on short maneuvers

**Actions Required**:
1. Add path duration or sample count to Path struct
2. Modify fitness calculation to normalize by duration
3. Verify all 6 paths contribute roughly equally to fitness variance
4. Test that evolution improves on both short and long paths

**Why Defer Until After Phase 2**:
- Requires changing fitness calculation (risky)
- Want deterministic baseline first
- Demetic mode must work correctly before path normalization
- Easy to test: compare fitness contribution variance across paths

## Lessons Learned: Attention to Detail

### 1. Comment Accuracy Matters
**Issue**: `autoc.cc:774` said "Throttle is 0-1 range" but was actually -1 to 1
- This misleading comment delayed diagnosis
- Always verify ranges by checking actual clamping code
- Cross-reference with sim interface (inputdev_autoc.cpp:764)

### 2. Domain Knowledge is Critical
**Issue**: Initially tried to use `abs()` for throttle, treating negative as "reverse thrust"
- User corrected: negative throttle means motor OFF (zero energy), not reverse
- Physical reality: -1 = idle, 0 = half power, +1 = full power
- Always validate assumptions against physical constraints

### 3. Don't Hide Bugs with Guards
**Issue**: Initially added `if (localFitness < 0.0f) localFitness = 0.0f;` guard
- User correctly identified: "a guard isn't right -- because this means a better solution is negative"
- Better approach: Trace root cause and fix the source
- Guards should only handle legitimate edge cases, not mask bugs

### 4. Floating-Point Precision vs. Real Bugs
**Issue**: `direction_alignment = 1.00000001` causing tiny negative errors
- This IS a legitimate use of clamping (mathematical boundary, not a bug)
- vs. `direction_alignment = 2.0` which would indicate serious normalization bug
- Know when to clamp (precision) vs. when to error (logic bug)

### 5. Verify Semantics Across System Boundaries
**Issue**: GP outputs throttle [-1,1], crrcsim expects [0,1], fitness calculated wrong range
- Always check interface contracts at system boundaries
- Cross-reference: GP output → sim input → fitness calculation
- One inconsistent interpretation breaks the entire chain

### 6. Test Incremental Changes
**Progression**:
1. Found NaN errors → investigated movement direction
2. Found negative fitness → traced to throttle calculation
3. Fixed throttle → verified against crrcsim interface
4. Confirmed minisim determinism → ready for crrcsim determinism phase

Each fix was tested before moving to next issue.

## Historical Context

### Previous Session Work (2025-12-31)
- Achieved 99.98% determinism in crrcsim (variance 20-50 → 0.005 units)
- Fixed RandGauss Box-Muller state persistence
- Reset eta1-4 Dryden turbulence RNGs
- Reset thermal generation RNGs
- Attempted SimStateHandler refactoring → REVERTED (increased variance to 120 units)
- Total changes: 36 lines across 3 files (crrc_rand.h, crrc_rand.cpp, windfield.cpp)

### This Session Work (2026-01-02)
- Fixed throttle energy calculation (negative fitness bug)
- Added dot product clamping (NaN prevention)
- Improved normalization_factor threshold (divide-by-near-zero prevention)
- Achieved perfect determinism in minisim (multi-threaded, non-demetic)
- Identified clear path forward for crrcsim determinism → demetic mode → path normalization

## Current Configuration for Testing

From `autoc.ini`:
```ini
PopulationSize                  = 1000
NumberOfGenerations             = 200
DemeticGrouping                 = 0        # Non-demetic mode for testing
WindScenarios                   = 1         # Single scenario for faster iteration
SimNumPathsPerGeneration        = 6
EvalThreads                     = 8         # Multi-threaded
PathGeneratorMethod             = aeroStandard
MinisimProgram                  = [commented out]  # Using minisim for determinism testing
```

## Debug Facilities Available

1. **Determinism Test** (autoc.cc:1531):
   - Change `if (false && gen == 2)` to `if (gen == 2)`
   - Clones best individual to entire population at generation 2
   - Measures fitness variance across 8 worker threads
   - Reports: mean, stddev, min, max fitness values

2. **RNG Tracing** (environment variable):
   ```bash
   export AUTOC_RNG_TRACE=1
   ```
   - Logs every rand() call with index, value, context
   - Useful for debugging PRNG state issues

3. **NaN Detection** (autoc.cc:920):
   - Increments `nanDetector` counter when NaN detected
   - Currently passive (no logging) - can add stderr output if needed

## Success Criteria

### Phase 1 Complete When:
- [ ] CRRCsim determinism test shows variance < 0.0001 units (same as minisim)
- [ ] 50 identical individuals produce identical fitness across all 8 threads
- [ ] Best fitness is strictly monotonically improving (no regressions)

### Phase 2 Complete When:
- [ ] Demetic mode shows monotonically improving fitness (no jumps)
- [ ] Elite's aggregated fitness is preserved across generations
- [ ] Fitness improves at similar rate to non-demetic mode

### Phase 3 Complete When:
- [ ] All 6 paths contribute roughly equally to fitness variance
- [ ] Short paths (12s) and long paths (30s) both show improvement
- [ ] Per-second error metrics are balanced across path types

## Final Note

This document captures the journey from "negative fitness WTF?" to "ah, throttle range was wrong all along" to "now let's achieve perfect determinism in crrcsim just like we have in minisim."

The key insight: **Fix the obvious bugs first (negative fitness), then tackle the subtle issues (0.005 variance), then optimize (path normalization).**

Each phase builds on the previous. Trying to fix demetic mode or path normalization while fitness could go negative would have been a nightmare to debug.

2026-01-03 13:02:52: <info> 5 251.692688 272.395206 303.057617     14 109.880000 53     7 7.320000 10
2026-01-03 13:03:46: <info> 5 262.849518 270.419659 296.679965     51 176.380000 52     10 7.660000 13
2026-01-03 13:04:16: <info> 5 260.415955 270.932515 310.773496     1491 175.040000 2     17 6.980000 2

2026-01-03 14:40:55: <info> 5 230.479385 273.530349 316.228165     305 196.400000 92     11 8.780000 14
2026-01-03 14:40:06: <info> 5 230.479385 273.530349 316.228165     305 196.400000 92     11 8.780000 14
2026-01-03 14:39:17: <info> 5 230.479385 273.530349 316.228165     305 196.400000 92     11 8.780000 14

## 2026-01-03 Afternoon: Multi-Process Non-Determinism Investigation

### Progress Made
1. ✅ **Disabled Wind/Turbulence** (autoc_config.xml:82): `velocity="0" turbulence="0"`
   - Eliminates Dryden turbulence RNG calls (`eta1-4.Get()`)
   - Early return in `calculate_gust()` at windfield.cpp:809 prevents RNG usage

2. ✅ **Fixed Initial RNG Seed** (crrc_main.cpp:500): Changed `srand(time(0))` → `srand(42)`
   - Ensures all processes start with same RNG state
   - Previous: Each process had different `time(0)` → different `rand()` sequence

3. ✅ **Valgrind Analysis**: No uninitialized variables in simulation code
   - 7,044 errors found, ALL in WSL graphics libraries (libdxcore.so, libd3d12core.so, etc.)
   - Zero uninitialized variable warnings in crrcsim FDM or simulation code
   - Created valgrind-graphics.supp to suppress graphics noise

### Current Status: Multi-Process Variance Still Present

**Single-Threaded (EvalThreads=1)**: ✅ PERFECT determinism
```
2026-01-03 14:40:55: <info> 5 230.479385 273.530349 316.228165
2026-01-03 14:40:06: <info> 5 230.479385 273.530349 316.228165
2026-01-03 14:39:17: <info> 5 230.479385 273.530349 316.228165
```
All stats identical down to the last digit across multiple runs.

**Multi-Process (EvalThreads=8)**: ❌ Still non-deterministic

Run 1:
- Variance: 1.5e-05 (0.0000064%)
- Groups: [3 individuals @ 230.479370] vs [47 @ 230.479385]
- Group correlation with worker IDs: 0,3,4 vs 1,2,5,6,7...

Run 2:
- Variance: 3.1e-02 (0.014%)
- Groups: [1 @ 228.340942], [4 @ 228.348877], [1 @ 228.348892], [44 @ 228.372330]
- Different grouping pattern than Run 1

**Key Observation**: Variance magnitude and grouping CHANGES between runs, and variance is often >3% in practice with larger samples.

### Why This Matters
- Minisim achieves perfect determinism with 8 worker processes
- CRRCsim shows variance with 8 worker processes
- Both use identical transport: autoc multi-threaded task dispatch to worker processes
- The 4th-order energy optimization signal is very small
- Current variance (up to 3%+) completely masks the energy cost signal
- Evolution cannot optimize throttle efficiency until variance is eliminated

### Root Cause: CRRCsim-Specific Non-Determinism

**Architecture**: Autoc is multi-threaded, dispatching tasks to worker processes via socket transport. Both minisim and crrcsim use identical process launch and communication. The difference is **internal to the simulators**.

**Not caused by**:
- ✅ RNG seeding (fixed with `srand(42)`)
- ✅ Wind/turbulence RNGs (disabled via velocity=0, turbulence=0)
- ✅ Uninitialized variables (valgrind found none in simulation code)
- ✅ Multi-process architecture (minisim uses same architecture, works perfectly)

**Likely candidates** (CRRCsim-specific):
1. **Per-process initialization differences**: CRRCsim initializes something differently per process
2. **Graphics/SDL initialization side-effects**: Even with graphics disabled, init may affect state
3. **Static/global FDM variables**: Not initialized consistently across processes
4. **Residual state from previous evaluations**: State not properly reset between evals in crrcsim

**Evidence**:
- Variance groups correlate with worker indices (0,3,4 vs 1,2,5,6,7...)
- Single crrcsim worker (EvalThreads=1, sequential evals) = perfect determinism
- Multiple crrcsim workers (EvalThreads=8, parallel evals) = variance
- Multiple minisim workers (EvalThreads=8, parallel evals) = perfect determinism
- **Conclusion**: The bug is in crrcsim, not in the multi-process architecture

### TODO: Systematic Elimination Process

**Phase 1A: Multi-Process Determinism (CURRENT BLOCKER)**
- [ ] **Next Step**: Compare DTEST logs between workers to find exact divergence point
- [ ] Check for static/global variables in FDM that aren't initialized consistently
- [ ] Add determinism logging to FDM state variables (position, velocity, quaternion)
- [ ] Check if process launch order or timing affects initial state
- [ ] Investigate floating-point environment (rounding modes, denormals, SSE flags)
- [ ] Compare minisim vs crrcsim: what's different in multi-threading vs multi-process?

**Phase 1B: Re-enable Wind/Turbulence** (after 1A complete)
- [ ] With multi-process determinism working, re-enable wind gradually
- [ ] Verify Dryden turbulence RNG resets are working correctly
- [ ] Test with `velocity > 0, turbulence = 0` first
- [ ] Then test with both enabled

**Phase 1C: Re-enable Thermals** (after 1B complete)
- [ ] Set thermal `density > 0` gradually
- [ ] Verify thermal generation RNGs are deterministic
- [ ] Check thermal position/strength consistency across processes

**Phase 2: Demetic Mode Elite Preservation** (after Phase 1 complete)
- [ ] Fix fitness tracking in demetic mode
- [ ] Implement elite skip or fitness restoration
- [ ] Test that best fitness is monotonically improving

**Phase 3: Path Length Normalization** (after Phase 2 complete)
- [ ] Normalize fitness by path duration
- [ ] Verify all 6 paths contribute equally

### Files Modified (2026-01-03)
1. `autoc_config.xml:82` - Disabled wind/turbulence
2. `crrc_main.cpp:500` - Fixed `srand(42)` instead of `srand(time(0))`
3. `crrcsim.sh:31-32` - Added valgrind suppressions
4. `valgrind-graphics.supp` - Created to suppress WSL graphics errors

### Critical Insight
The fact that we get **different variance amounts AND groupings** between runs suggests this is NOT just FPU differences or compiler optimizations. Those would produce consistent groupings. Something varies per process launch that affects simulation results.

## 2026-01-03 Late Afternoon: DTEST Analysis - Divergence During Simulation

### Test Run (16:18-16:19, 8 workers)
- Variance: 0.031 (3.1%)
- 3 fitness groups: [1 @ 228.340942], [5 @ 228.348877], [44 @ 228.372330]
- Base PID: 3552359

### DTEST_RESET Analysis: Perfect Initialization Across All Workers

All 8 worker processes show **IDENTICAL** initial state:
```
Position:   [0.000000, 0.000000, -25.001221]  ✅ All 8 workers
Quaternion: [0.000000, 0.000000,   0.000000,   1.000000]  ✅ All 8 workers
Velocity:   [-7.015587, 0.000000,  -0.000000]  ✅ All 8 workers
PRNG state: [4337, 4337]  ✅ All 8 workers
```

**Verified via**: `/home/gmcnutt/GP/autoc/analyze_dtest.sh`
- Extracted position, quaternion, velocity from all 8 worker DTEST_RESET logs
- Used `sort | uniq -c` to confirm all values identical
- 8 unique workers, 1 unique value for each state variable

### Key Finding: **Initialization is Deterministic, Simulation is Not**

**What This Eliminates:**
- ❌ ~~Per-process RNG state differences~~ (all start with [4337, 4337])
- ❌ ~~Per-process initialization bugs~~ (position/velocity/quaternion identical)
- ❌ ~~Graphics/SDL init side-effects on state~~ (would affect initial values)
- ❌ ~~Process launch order affecting initial state~~ (all identical)

**What This Proves:**
- ✅ The bug is in the **FDM simulation loop**, not initialization
- ✅ All workers start from identical state but diverge during integration
- ✅ The divergence accumulates over the ~6-second flight simulation

### Why Minisim Works But CRRCsim Doesn't

Both use identical:
- Multi-process architecture (8 worker processes)
- Transport protocol (socket-based task dispatch)
- Initial state (same position/velocity/quaternion)
- RNG seeding (same windSeed)

**The difference must be in the FDM integration loop itself:**
- Minisim: Custom simple FDM, likely uses basic Euler integration
- CRRCsim: Complex FDM (EOM01) with gear model, aerodynamics, etc.

Hypothesis: CRRCsim's FDM has some per-process state that minisim doesn't have.

## 2026-01-03 16:30-16:31: CRITICAL BREAKTHROUGH - DTEST_FIRST Analysis

### Fixed DTEST_FIRST Logging

**Problem**: DTEST_FIRST logging code existed but never triggered because it checked for `pathIndex == 0 && simTimeMsec == 0`, but by the time `update()` is first called, `simTimeMsec` has already advanced past 0.

**Solution** ([inputdev_autoc.cpp:469-470](src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp#L469-L470), [inputdev_autoc.h:110](src/mod_inputdev/inputdev_autoc/inputdev_autoc.h#L110)):
1. Added member variable `bool needDtestFirstLog` to class
2. Set flag to `true` in reset function (after DTEST_RESET log)
3. Check flag in update() and log first timestep, then clear flag
4. Also updated [analyze_dtest.sh](analyze_dtest.sh) to analyze DTEST_FIRST entries

### Test Run (16:30-16:31, 8 workers)
- Variance: 0.023 (2.3%)
- 2 fitness groups: [7 @ 228.348877], [43 @ 228.372330]
- Base PID: 3555579

### **CRITICAL FINDING: Workers Diverge Within First 117ms!**

**DTEST_RESET at t=0** (all 8 workers **IDENTICAL**):
```
Position:   [0.000000,  0.000000, -25.001221]
Quaternion: [0.000000,  0.000000,   0.000000,  1.000000]
Velocity:   [-7.015587, 0.000000,  -0.000000]
PRNG:       [4337, 4337]
```

**DTEST_FIRST at t=117ms** (workers **ALREADY DIVERGED**):

Worker 3555594: `pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000223,0.999954]` ✅
Worker 3555603: `pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000223,0.999954]` ✅
Worker 3555611: `pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000223,0.999954]` ✅
Worker 3555619: `pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000223,0.999954]` ✅
Worker 3555627: `pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000224,0.999954]` ⚠️ quat[2] differs
Worker 3555634: `pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000223,0.999954]` ✅
Worker 3555643: `pos[-0.811070,-0.000000,-24.952831] quat[0.000000,0.009610,0.000101,0.999954]` ❌ **BOTH pos AND quat differ!**
Worker 3555651: `pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000223,0.999954]` ✅

### Analysis

**Divergence magnitude at t=117ms:**
- Worker 3555627: Quaternion component differs by **0.000001** (1 in the last decimal place)
- Worker 3555643: Position differs by **~0.0002m**, quaternion by **~0.0001**

**Pattern:**
- 6 workers match perfectly (majority group)
- 1 worker (3555627) has minor quaternion difference
- 1 worker (3555643) has both position AND quaternion differences

**This proves:**
1. ✅ Divergence happens **WITHIN THE FIRST 117ms** of FDM integration
2. ✅ The non-determinism is **NOT accumulation error** - it starts immediately
3. ✅ Something in the FDM integration loop produces different results **per process**
4. ✅ Worker 3555643 shows the largest divergence - might correlate with which fitness group it's in

### Next Steps
- Add more frequent DTEST logging (every 10-20ms) to pinpoint exact frame of divergence
- Investigate EOM01 FDM integration loop for per-process state
- Check for floating-point environment differences (denormals, rounding modes)
- Look for static/global variables in FDM that aren't reset properly
- Compare worker PIDs to fitness groups to see if same workers consistently diverge

## 2026-01-03 16:45: EUREKA - First-Run Initialization Bug Identified!

### Discovery Process

Noticed that workers 3555594 and 3555643 had **175 log lines** while others had **160 lines**, suggesting they did more work. Checked `DTEST_FIRST path=0` entries across all workers:

**Critical Pattern Found:**
- Each worker evaluates path=0 multiple times (6 or 7 times for 50-individual population)
- Worker's **FIRST** evaluation of path=0: Shows UNIQUE divergent values
- Worker's **SUBSEQUENT** evaluations of path=0: All converge to **IDENTICAL** stable value `pos[-0.811070,-0.000000,-24.952831] quat[0.000000,0.009610,0.000101,0.999954]`

### Evidence

**Worker 3555594 (7 evaluations of path=0):**
```
Eval 1: pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000223,0.999954]  ← UNIQUE
Eval 2-7: pos[-0.811070,-0.000000,-24.952831] quat[0.000000,0.009610,0.000101,0.999954]  ← STABLE
```

**Worker 3555627 (6 evaluations of path=0):**
```
Eval 1: pos[-0.810833,-0.000002,-24.952824] quat[0.000002,0.009614,0.000224,0.999954]  ← UNIQUE (quat[2] differs!)
Eval 2-6: pos[-0.811070,-0.000000,-24.952831] quat[0.000000,0.009610,0.000101,0.999954]  ← STABLE
```

**Worker 3555643 (7 evaluations of path=0):**
```
Eval 1-7: pos[-0.811070,-0.000000,-24.952831] quat[0.000000,0.009610,0.000101,0.999954]  ← STABLE from start!
```

### Root Cause: **Uninitialized State on First Simulation Run**

**Why Single-Threaded Works:**
- One crrcsim worker process
- First evaluation (1/50): Produces "first-run" result with uninitialized state
- Remaining evaluations (2-50/50): Use stable initialized state
- **Determinism test happens after 5 generations**, so all 50 clones use stable state
- Result: Perfect determinism

**Why Multi-Process Fails:**
- 8 crrcsim worker processes, each freshly launched
- Each worker's first evaluation uses uninitialized state
- Population of 50 distributed across 8 workers: ~6-7 evals per worker
- **During determinism test**: Each worker is doing early evaluations, some are still on "first run"
- Result: Mix of "first-run unstable" and "warmed-up stable" results → variance

### Worker Evaluation Counts (path=0 during determinism test)
```
Worker 3555594: 7 evaluations (did more work)
Worker 3555603: 6 evaluations
Worker 3555611: 6 evaluations
Worker 3555619: 6 evaluations
Worker 3555627: 6 evaluations
Worker 3555634: 6 evaluations
Worker 3555643: 7 evaluations (did more work)
Worker 3555651: 6 evaluations
```

Workers 3555594 and 3555643 handled more evaluations, explaining their larger log files (175 vs 160 lines).

### What This Means

1. ✅ **NOT** floating-point environment differences
2. ✅ **NOT** race conditions or threading bugs
3. ✅ **NOT** random number generator issues
4. ✅ **IS** uninitialized state in crrcsim that gets initialized after first simulation run
5. ✅ Problem is in crrcsim initialization, NOT in autoc or the multi-process architecture

### Next Investigation

Find what state in crrcsim/EOM01 FDM is:
- Uninitialized at process startup
- Gets initialized during first `Simulation->reset()` or first FDM integration
- Remains stable for subsequent evaluations

Candidates:
- Static variables in FDM that are initialized on first use
- Lazy initialization in EOM01 or aerodynamics calculations
- State that persists across `Simulation->reset()` calls
- Graphics/OpenGL state that affects FDM calculations somehow

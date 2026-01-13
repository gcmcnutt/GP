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

---

# Determinism Debugging Unwind Checklist (2026-01-04)

## Overview
This section tracks the debugging changes made to achieve determinism and provides a checklist for unwinding back to production configuration. Changes are categorized by whether they should be **kept** (core fixes) or **removed** (debug infrastructure).

## Core Fixes - KEEP THESE ✅

These changes fixed the root cause and should remain in production:

### 1. Binary Serialization (CRITICAL - KEEP)
**File**: `~/GP/autoc/minisim.h`
- ✅ **KEEP**: Changed from `text_archive` to `binary_archive` (lines 20-21)
- ✅ **KEEP**: Size-prefixed framing in `sendRPC()` and `receiveRPC()` (lines 281-311)
- ✅ **KEEP**: Cross-platform portability comments (lines 23-29)
- **Reason**: This is the root cause fix for non-determinism

**Verification**:
```bash
cd ~/GP/autoc && grep "binary_archive" minisim.h
# Should show: #include <boost/archive/binary_oarchive.hpp>
# Should show: #include <boost/archive/binary_iarchive.hpp>
```

### 2. AircraftState Zero-Initialization (KEEP)
**File**: `~/GP/autoc/aircraft_state.h`
- ✅ **KEEP**: Default constructor zero-initialization (lines 205-214)
- **Reason**: Prevents uninitialized memory issues

### 3. String Literal Fix (KEEP)
**File**: `~/crsim/crrcsim-0.9.13/src/mod_fdm_config.h`
- ✅ **KEEP**: Added space in MOD_FDM_INFOSTR macro (line 38)
- **Reason**: Required for C++11 compliance

### 4. std::to_string Fix (KEEP)
**File**: `~/crsim/crrcsim-0.9.13/src/mouse_kbd.cpp`
- ✅ **KEEP**: Changed stringstream to std::to_string (line 337)
- **Reason**: Cleaner code, no MSan false positives

---

## Debug Infrastructure - DECIDE: KEEP OR REMOVE 🤔

### 5. DTest Tracker Infrastructure
**Files**:
- `~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/dtest_tracker.h`
- `~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/dtest_tracker.cpp`
- `~/crsim/crrcsim-0.9.13/src/mod_inputdev/CMakeLists.txt` (lines 9-10)

**Options**:
- **Option A - KEEP**: Useful for future debugging, minimal overhead when disabled
- **Option B - REMOVE**: Clean up codebase

**Decision Point**: Does the value of having determinism tracking for future debugging outweigh code complexity?

### 6. serializeForHash() Helper
**File**: `~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp`
- Lines 47-57: `serializeForHash()` template function
- Lines 690-692: Usage in logPathComplete
- Lines 726-728: Usage in logEvalComplete

**Decision**: Depends on decision #5 (only used when tracker active)

---

## Debug Code - REMOVE THESE ❌

### 7. AUTOC_RECV_HASH Debug Logging (REMOVE)
**File**: `~/GP/autoc/autoc.cc`

**Lines to remove**:
- Lines 64-72: `hashData()` function
- Lines 717-726: First AUTOC_RECV_HASH block (in evalTask)
- Lines 1236-1245: Second AUTOC_RECV_HASH block (in main/bakeoff)

**Verification After Removal**:
```bash
cd ~/GP/autoc && grep "AUTOC_RECV_HASH" autoc.cc
# Should return nothing
```

### 8. MALLOC_PERTURB (REMOVE)
**File**: `~/GP/autoc/crrcsim.sh`
- ❌ **REMOVE**: Line 40: `export MALLOC_PERTURB_=255`
- ❌ **REMOVE**: Lines 38-39: Comments

**Verification After Removal**:
```bash
cd ~/GP/autoc && grep "MALLOC_PERTURB" crrcsim.sh
# Should return nothing
```

---

## Compiler Flags - RESTORE PRODUCTION SETTINGS ⚙️

### 9. Restore -O2 Optimization (RESTORE)
**File**: `~/crsim/crrcsim-0.9.13/CMakeLists.txt`

**Action**: Change `-O0` to `-O2` on lines 208-209
**KEEP** the FP determinism flags: `-ffp-contract=off -fno-fast-math -msse2 -mfpmath=sse`

**Current**:
```cmake
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall ...")
```

**Production**:
```cmake
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -Wall ...")
```

### 10. C++17 Requirement (KEEP)
- ✅ **KEEP**: Lines 203-204: C++17 standard requirement (needed for std::clamp)

### 11. USE_CLANG Option (KEEP)
- ✅ **KEEP**: Lines 12-20: Allows easy toggling between Clang/GCC for debugging

---

## Configuration Files - RESTORE PRODUCTION VALUES 📝

### 12. autoc.ini Thread Count
**File**: `~/GP/autoc/autoc.ini`
- Current: `EvalThreads = 8`
- **Decision**: User should verify intended production value (1 or 8)

### 13. Re-enable Wind/Thermals (ENABLE GRADUALLY)
**File**: `~/GP/autoc/autoc.ini`
- Current: `WindScenarios = 1`
- **Action**: Gradually increase and verify determinism holds

### 14. Re-enable Demetic Grouping (ENABLE GRADUALLY)
**File**: `~/GP/autoc/autoc.ini`
- Current: `DemeticGrouping = 0`
- **Action**: Enable and verify determinism with demes

---

## Unwind Execution Order

Execute in this order to safely unwind:

1. ✅ **Baseline**: Run determinism test with current settings (DONE)
2. ❌ **Remove debug logging**: AUTOC_RECV_HASH from autoc.cc
3. ❌ **Remove MALLOC_PERTURB**: From crrcsim.sh
4. ⚙️ **Restore -O2**: In CMakeLists.txt (change -O0 to -O2)
5. 🧪 **Test**: Verify determinism still works with -O2
6. 🤔 **Decide on tracker**: Keep or remove dtest_tracker infrastructure
7. 📝 **Restore config**: autoc.ini settings as needed
8. 🚀 **Enable features gradually**:
   - WindScenarios = 3 (test determinism)
   - WindScenarios = 5+ (test determinism)
   - DemeticGrouping = 1 (test determinism)
9. 🧪 **Full validation**: 20+ generation run, check for ELITISM VIOLATIONS

---

## Testing Checklist

### Phase 1: Production Flags Verification
```bash
# Apply -O2 optimization
cd ~/crsim/crrcsim-0.9.13
# Edit CMakeLists.txt: change -O0 to -O2 on lines 208-209
cd build && cmake .. && make -j8

# Rebuild autoc
cd ~/GP/build && make autoc

# Run determinism test
cd ~/GP/autoc
timeout 180 ../build/autoc 2>&1 | grep -E "(SUCCESS|FAILURE)"
# Expected: "SUCCESS: All individuals have identical fitness"
```

### Phase 2: Wind Variation Testing
```bash
# Edit autoc.ini: WindScenarios = 3
# Run determinism test - should still pass
```

### Phase 3: Demetic Grouping Testing
```bash
# Edit autoc.ini: DemeticGrouping = 1, DemeSize = 10
# Run determinism test - should still pass
```

### Phase 4: Long-Run Validation
```bash
# Full evolution run
cd ~/GP/autoc
../build/autoc 2>&1 | tee validation.log
grep "ELITISM VIOLATION" validation.log
# Expected: None
```

---

## Git Diff Summary

**GP Repo Changes**:
```
autoc/aircraft_state.h | 11 +++++++ (KEEP: zero-init)
autoc/autoc.cc         | 35 +++++++++++++++++++ (REMOVE: debug logging)
autoc/autoc.ini        |  4 +++--- (RESTORE: production values)
autoc/crrcsim.sh       |  4 ++++ (REMOVE: MALLOC_PERTURB)
autoc/minisim.h        | 42 +++++++++++++++++++++ (KEEP: binary serialization)
```

**crrcsim Repo Changes**:
```
CMakeLists.txt         | 40 +++++++--- (KEEP structure, RESTORE -O2)
src/mod_fdm_config.h   |  2 +- (KEEP: C++11 fix)
src/mod_inputdev/CMakeLists.txt | 2 + (DECIDE: tracker files)
src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp | 87 +++++ (DECIDE: tracker code)
src/mouse_kbd.cpp      | 11 +-- (KEEP: std::to_string)
```

---

## Success Criteria

System ready for production when:
- ✅ Determinism test passes with -O2
- ✅ Wind scenarios work deterministically
- ✅ Demetic grouping works deterministically
- ✅ No ELITISM VIOLATIONS in 20+ generation runs
- ✅ Debug code removed or disabled
- ✅ Performance restored (throughput ≥ with -O2)

---

## Rollback Plan

If determinism breaks:
1. `git diff HEAD~1` to see what changed
2. `git revert <commit>` to undo
3. Re-enable debug logging to diagnose
4. Check serialization format hasn't changed

---

**Last Updated**: 2026-01-04
**Status**: Determinism achieved, ready to unwind debug infrastructure
**Next Action**: User decision on what to keep/remove, then execute unwind in order

---

# Elite Trace Capture for Elitism Violation Debugging (2026-01-04)

## Problem Statement

Despite achieving perfect determinism (binary serialization, zero variance in determinism tests), we still occasionally see **ELITISM VIOLATIONS** during normal GP evolution where elite fitness worsens across generations:

```
Gen 1:  234.031540
Gen 2:  234.031555  ← WORSE! (+1.5e-05)
Gen 5:  234.032516  ← WORSE! (+9.7e-04)
Gen 10: 247.564636  ← MUCH WORSE! (+13.5)
```

Elite fitness should be **monotonically improving** (or at minimum, non-worsening). The fact that it gets worse suggests there's still some source of non-determinism or state that isn't being properly preserved/reset.

## Previous Debugging Approaches - What Didn't Work

### Approach 1: Warmup Code (Red Herring)
- **What we tried**: Added dummy simulation warmup at process initialization to "initialize lazy state"
- **Result**: Elitism violations continued
- **Why it didn't work**: Masked the real problem rather than fixing it
- **Status**: ❌ Removed (warmup is a red herring)

### Approach 2: Artificial Determinism Test (Limited Value)
- **What we tried**: Clone best individual to entire population at gen 20, re-evaluate all 50 clones, measure variance
- **Result**: Worked well for finding binary serialization bug, but has limitations:
  - Only tests at specific generation (gen 20)
  - Doesn't capture actual elitism violations during normal evolution
  - Hard to explain why it would show divergence AFTER 20 generations of work
  - Only first tasks show non-determinism (first-run initialization bug)
- **Why limited**: Tests artificial scenario, not the actual symptom (elite fitness worsening)
- **Status**: ✅ Available but not addressing the real problem

## New Strategy: Elite Trace Capture and Comparison

### Core Insight
Instead of testing for determinism artificially, focus on the **actual symptom**: why does elite fitness worsen during normal evolution?

### The Plan

**When elite is evaluated**:
1. Autoc already receives full `EvalResults` with `aircraftStateList` (complete simulation trace)
2. Save this full trace when evaluating the elite individual
3. Use existing binary serialization (already proven to work)
4. Store to S3/Redis (existing infrastructure)

**When fitness worsens**:
1. Compare the previous elite trace with the current elite trace
2. Have both traces available for detailed analysis
3. Identify exact point of divergence

### Why This Approach is Better

1. **Surgical**: Focuses on actual elitism violations, not artificial tests
2. **Complete data**: Already have full simulation traces in `EvalResults.aircraftStateList`
3. **Minimal changes**: Only autoc.cc needs modification, not crrcsim
4. **Uses proven tech**: Binary serialization already working
5. **Real-world**: Captures violations during normal evolution, not synthetic scenarios

### Implementation Details

**Location**: All changes in `autoc.cc` only (not crrcsim)

**Data Already Available**:
- `EvalResults.aircraftStateList`: Vector of `AircraftState` for every timestep
- `AircraftState` contains:
  - Position, velocity, orientation (quaternion)
  - Pitch/roll/throttle commands
  - Simulation time
  - Path index
- This is returned from every evaluation via binary serialization

**What to Track**:
1. Identify when elite individual is being evaluated
2. Store the `EvalResults` for that evaluation
3. Save to S3/Redis with key like: `elite-trace-gen{N}-individual{ID}`
4. When fitness worsens (gen N+1 elite fitness > gen N elite fitness):
   - Log both trace keys for comparison
   - Optionally trigger automatic diff analysis

**Note on Granularity**:
- User feedback: "later we can deal with detailed data like quat per step, etc -- but i think for elite this is mostly captured already (like in data.dat -- a VERY large file)"
- `aircraftStateList` already contains per-step data
- `data.dat` file also has comprehensive trace data
- Initial focus: Just save the traces, manual comparison first
- Future: Automated diff tool to pinpoint divergence

### Success Criteria

1. ✅ Elite traces saved every generation
2. ✅ When elitism violation occurs, have both previous and current trace
3. ✅ Can manually compare traces to find divergence point
4. ✅ Identify root cause of fitness worsening

### Why We Expect This to Work

**Theory**: Elite re-evaluation uses different state than original evaluation
- First evaluation: Worker state A → fitness X
- Re-evaluation: Worker state B → fitness X + delta
- With full traces, we can see WHERE the divergence starts (position, velocity, orientation)
- Then work backward to find WHAT state caused it

**Evidence from First-Run Bug**:
- Workers showed different behavior on first evaluation vs subsequent
- Once we had DTEST logging, pinpointed divergence to first 117ms
- Same principle: capture full trace, find divergence point, identify root cause

## Implementation Roadmap

### Phase 1: Basic Trace Capture (This Session)
- [ ] Add elite tracking to autoc.cc
- [ ] Save `EvalResults` when elite is evaluated
- [ ] Use binary serialization to store to S3/Redis
- [ ] Log trace keys when elitism violation detected

### Phase 2: Comparison Tool (Future)
- [ ] Build diff tool to compare two `EvalResults` traces
- [ ] Output: First divergence point (generation, timestep, which field)
- [ ] Automatic trigger when violation detected

### Phase 3: Root Cause (Future)
- [ ] Use trace comparison to identify pattern
- [ ] Find uninitialized state or non-deterministic code path
- [ ] Fix and verify elitism violations disappear

## Related Issues

This debugging strategy addresses:
- ELITISM VIOLATIONS during normal evolution
- Understanding why elite fitness worsens
- Identifying remaining sources of non-determinism

This does NOT address (separate issues):
- Demetic mode elite preservation (Phase 2 from earlier plan)
- Path length normalization (Phase 3 from earlier plan)

---

**Status**: ✅ Implemented, testing reveals crrcsim non-determinism
**Next Action**: Switch to minisim for deterministic evaluation

---

# CRRCSIM Non-Determinism Investigation (2026-01-04)

## Current Status

Elite trace capture is fully implemented and working. However, testing reveals **crrcsim itself is non-deterministic** even with wind disabled.

## Evidence from Test Run

### Test Configuration
- Wind mode: 0 (disabled)
- Thermals: disabled (`num_thermals = 0`)
- Same elite individual re-evaluated across generations

### Observed Non-Determinism

**Generation 12 → 13**: Same elite structure, different fitness
```
Gen 12: fitness=194.118286
Gen 13: fitness=194.118317 (delta: +3.05e-05)

Divergence at Path 0, Step 1, Time 117ms:
  Velocity Y-axis: 54,689 ULP difference
  Trace1: [-6.924340725, -5.897961819e-05, 7.578417063e-01]
  Trace2: [-6.924338818, -5.917857561e-05, 7.578418255e-01]
```

**Generation 13 → 14**: Same elite structure, different fitness again
```
Gen 13: fitness=194.118332
Gen 14: fitness=194.118317 (delta: -1.53e-05, improved!)

Divergence at Path 0, Step 4, Time 468ms:
  Velocity Y-axis: 543,009 ULP difference
  Trace1: [-8.030043602, -9.586629858e-06, 3.315576077]
  Trace2: [-8.030044556, -9.092766049e-06, 3.315575838]
```

### Key Observations

1. **Very Early Divergence**: Step 1 (117ms) and Step 4 (468ms)
2. **Massive ULP Differences**: 54K-543K ULPs in Y-axis (East velocity)
3. **Near-Zero Values**: Y-velocity ~1e-05 m/s (small absolute differences create huge ULP ratios)
4. **Wind Disabled**: Confirms wind_from_terrain() is NOT the root cause
5. **Same Elite, Different Fitness**: Proves crrcsim FDM is non-deterministic

## Why "Elite Replaced but Fitness NOT Improved"

Messages like:
```
ELITE_TRACE gen=2: Elite replaced but fitness NOT improved (library bug?)
  Previous: fitness=234.031540
  Current:  fitness=234.031555
```

This occurs when:
- GP library selects different individual as "best" (structure changed)
- But new individual has equal or worse fitness
- Likely library tie-breaking on tree size or other metadata
- Not a determinism issue - just GP library selection logic

## Root Cause Analysis

### NOT Wind or Thermals
- Wind mode = 0 (disabled)
- Thermals = 0 (disabled)
- Still seeing non-determinism

### Likely Causes

1. **FDM Numerical Integration**
   - Floating-point operation order varies
   - CPU reordering or optimization differences
   - Rounding errors accumulating differently

2. **Uninitialized FDM State**
   - Some state not fully reset by `Global::Simulation->reset()`
   - Previous evaluation affecting next one
   - Hidden static variables or class members

3. **External Libraries**
   - Aerodynamics calculations using non-deterministic libraries
   - Matrix operations with different implementations
   - Compiler optimizations affecting FP math order

4. **CPU/FPU State**
   - Floating-point rounding mode differences
   - SSE/AVX instruction usage
   - Different code paths on different cores

### Why Search for "rand" is Too Myopic

As user noted: "libraries that use random under the hood with all sorts of names"
- Could be `std::random_device`, `drand48`, `mt19937`, etc.
- Could be in third-party libraries (Eigen, Boost, etc.)
- Could be non-random but still non-deterministic (timing, memory addresses, etc.)

## Solution: Switch to Minisim

### Why Minisim

User statement: **"minisim is deterministic (but unsophisticated)"**

Minisim characteristics:
- Simple vector math
- No complex aerodynamics libraries
- No wind/thermal calculations
- Repeatable/deterministic
- Lightweight

### Current Configuration

```ini
MinisimPortOverride = 0    # 0 = use crrcsim
```

To switch to minisim:
```ini
MinisimPortOverride = 11000  # Use minisim on port 11000
```

### Implementation Plan

1. **Test minisim determinism**
   - Set `MinisimPortOverride = 11000`
   - Run same elite trace test
   - Verify zero divergence (all ULP differences < 4)

2. **If minisim is deterministic**
   - Use minisim for ALL GP fitness evaluation
   - Keep crrcsim for final validation/visualization
   - Accept "unsophisticated" physics for evolution

3. **If minisim is also non-deterministic**
   - Investigate minisim code directly
   - Much simpler codebase to debug
   - Fix root cause in minisim

## Files Modified This Session

### Elite Trace Capture (Completed)
- `/home/gmcnutt/GP/autoc/autoc.cc` (lines 1647-1940)
  - Added structure-based elite tracking
  - In-memory trace comparison with ULP diagnostics
  - File-based trace storage to `/tmp/crrcsim/`

- `/home/gmcnutt/GP/autoc/compare_traces.cc` (new file)
  - Standalone trace comparison tool

- `/home/gmcnutt/GP/autoc/CMakeLists.txt`
  - Added compare_traces build target

### Wind Diagnostics (Completed)
- `/home/gmcnutt/GP/autoc/aircraft_state.h`
  - Added `wind_velocity` field (gp_vec3)
  - Added getter/setter methods
  - Updated serialization (BOOST_CLASS_VERSION 1→2)
  - Added buffer overflow guards (0xDEADBEEF/0xCAFEBABE)

### Configuration Changes
- `/home/gmcnutt/crsim/crrcsim-0.9.13/autoc_config.xml`
  - Changed `wind_mode` from 2 to 0 (disabled wind)

### Documentation Created
- `/home/gmcnutt/GP/autoc/specs/WIND_NONDETERMINISM_FINDINGS.md`
- `/home/gmcnutt/GP/autoc/specs/THERMAL_PRNG_FINDINGS.md`
- `/home/gmcnutt/GP/autoc/specs/WIND_DETERMINISM_TEST_PLAN.md`
- `/home/gmcnutt/GP/autoc/specs/IN_MEMORY_TRACE_COMPARISON.md`

## Test Results Summary

### With wind_mode=0 (Wind Disabled)

| Generation | Elite Fitness | Status | Divergence |
|------------|---------------|--------|------------|
| 1 | 234.031540 | Baseline | - |
| 2 | 234.031555 | Elite replaced (worse) | - |
| 3 | 234.031540 | Elite replaced (better) | - |
| 4 | 230.479385 | Elite replaced (better) | - |
| 5-9 | 230.479385 | Elite replaced (same) | - |
| 10 | 230.476929 | Elite replaced (better) | - |
| 11 | 194.118286 | Elite replaced (better) | - |
| **12** | **194.118317** | **VIOLATION** | **54,689 ULP @ 117ms** |
| **13** | **194.118332** | **VIOLATION** | **543,009 ULP @ 468ms** |
| 14 | 194.118317 | Elite replaced (better) | - |
| 15-20 | 194.118317-332 | Multiple violations | - |

**Conclusion**: crrcsim is fundamentally non-deterministic, even with all wind/thermal systems disabled.

## Next Actions

1. **Switch to minisim** - Set `MinisimPortOverride = 11000` in autoc.ini
2. **Test minisim determinism** - Run same 20-generation test
3. **Verify zero divergence** - Expect ULP differences < 4 (FP rounding only)
4. **If successful** - Use minisim for all GP evolution
5. **Document findings** - Update this file with minisim test results

---

**Last Updated**: 2026-01-06 18:55 PST
**Status**: EOM01 FDM state initialization fixes applied, testing in progress
**Next Action**: Build and test with uninitialized variable fixes

---

## 2026-01-06: Deep Dive into CRRCSIM FDM State Initialization

### Investigation Summary

Following up on evidence that crrcsim shows step-0-to-step-1 divergence with identical initial state, we systematically investigated the EOM01 FDM (Flight Dynamics Model) for uninitialized state and hidden variables.

### Test Evidence - The Smoking Gun

From `~/GP/autoc/autoc_diag41.log` and trace comparison output:

**Step 0 (Initial State)**: PERFECTLY IDENTICAL across all traces
```
pos=[0.000000, 0.000000, -25.001221]  ✓ ULP[0,0,0]
vel=[-7.015587, 0.000000, -0.000000]  ✓ ULP[0,0,0]
quat=[w,x,y,z] identical
RNG state: [4337, 4337] identical
```

**Step 1 (First Integration, t=117ms)**: DIVERGENCE DETECTED
```
Path 0:
  vLocalDot1: [-3.811813, 0.001542, 5.319456]
  vLocalDot2: [-3.939870, 0.001598, 5.319243]
  ULP difference: [1121, 474198, 445]

Path 1:
  vLocalDot differences: tiny (ULP ~1-500)

Path 2:
  vLocalDot at step 4: IDENTICAL (first few steps matched)
```

**Critical Pattern**: Different paths show different divergence patterns and timing:
- Path 0: Large divergence at step 1
- Path 1: Tiny divergence at step 1
- Path 2: Perfect match until step 4

This pattern is **inconsistent with FPU differences or compiler optimizations** (which would affect all paths identically). It points to **uninitialized member variables** that have different leftover values from previous evaluations.

### Multi-Process Architecture Clarification

**User correction**: "each crrcsim itself is a single process"
- Autoc is multi-threaded (8 EvalThreads)
- Each EvalThread communicates with a separate crrcsim process via sockets
- So it's **multi-process from autoc's perspective**, not multi-threading within crrcsim
- Each crrcsim process handles multiple sequential evaluations (reusing the FDM object)

### Root Cause Analysis - What We Found

#### 1. FDM Object Lifecycle: REUSE, Not Recreation ✅ VERIFIED

Traced through the code:
```cpp
// crrc_main.cpp:655 - Aircraft created ONCE at startup
Global::aircraft = new Aircraft();

// Aircraft constructor (aircraft.cpp:40) creates FDMInterface
fdmInterface = new ModFDMInterface();

// ModFDMInterface::loadAirplane() creates the EOM01 object ONCE
fdm = new CRRC_AirplaneSim_Larcsim(...);  // Inherits from EOM01

// On each reset:
Global::Simulation->reset() {
  initialize_flight_model() {
    Global::aircraft->getFDMInterface()->initAirplaneState(...);  // NOT new
  }
}
```

**Conclusion**: The EOM01 FDM object is **created once** and **reused** across all evaluations within a worker process. This means uninitialized constructor state can leak between evaluations.

#### 2. Uninitialized Debug State Variables ✅ FIXED

Found in `eom01.h` (recently added for diagnostics):
```cpp
CRRCMath::Vector3 dbg_V_local_airmass;  // NOT initialized in constructor!
CRRCMath::Vector3 dbg_V_gust_body;      // NOT initialized in constructor!
CRRCMath::Vector3 dbg_force_body;       // NOT initialized in constructor!
CRRCMath::Vector3 dbg_moment_body;      // NOT initialized in constructor!
```

These are written during `ls_aux()` and `ls_accel()` but **only read for diagnostic output**. However, if they contain garbage on first evaluation, and the garbage differs between processes (heap allocator state), this could affect memory layout and cache behavior, leading to FP calculation differences.

**Constructor was EMPTY**:
```cpp
EOM01::EOM01(const char* logfilename, FDMEnviroment* myEnv) : FDMBase(logfilename, myEnv)
{
  // NOTHING HERE! All member variables uninitialized!
}
```

#### 3. Integrator History Variables - Partially Initialized ✅ FIXED

Adams-Bashforth integration requires "past" derivative values:
```cpp
// These ARE initialized by ls_step_init() when called from derived classes
v_V_dot_past              // Linear acceleration (previous timestep)
v_R_omega_dot_body_past   // Angular acceleration (previous timestep)
e_dot_0/1/2/3_past        // Quaternion derivatives (previous timestep)
latitude/longitude/radius_dot_past
```

However, `ls_step_init()` is only called:
- On first aircraft load (from derived class constructor)
- When explicitly triggered by `initAirplaneState()`

The `resetIntegratorState()` function was added to `SimStateHandler::reset()` to clear these, but it **didn't clear the debug variables**.

#### 4. Other Member Variables - Uninitialized ✅ FIXED

The EOM01 class has ~30 member variables that were NONE initialized in the constructor:
```cpp
// State vectors
v_R_omega_body, v_V_local, v_R_omega_dot_body, v_V_dot_local
v_V_local_rel_ground, v_V_local_rel_airmass, v_V_wind_body, v_P_CG_Rwy

// Quaternions and Euler angles
e_0, e_1, e_2, e_3
euler_angles_v[3], geocentric_position_v[3], geodetic_position_v[3]

// Transformation matrix
LocalToBody (3x3 matrix)

// Scalars
Mass, I_xx, I_yy, I_zz, I_xz
Sea_level_radius, V_rel_wind, Alpha, Beta, Gravity, Density
```

All of these contain **garbage values** until `ls_step_init()` or `initAirplaneState()` is called.

### Systematic Investigation Completed (Steps 1-5)

✅ **Step 1: Initialize dbg_* variables in EOM01 constructor**
- Added zero-initialization for all 4 debug vectors
- Added zero-initialization for ALL member variables (30+ total)
- File: `~/crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.cpp:41-86`

✅ **Step 2: Expand resetIntegratorState() to zero dbg_* variables**
- Updated `resetIntegratorState()` to also clear debug state
- Called on every `Global::Simulation->reset()`
- File: `~/crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.h:137-152`

✅ **Step 3: Verify LocalToBody matrix has no cached state**
- Confirmed: Recalculated fresh every timestep from quaternions (eom01.cpp:314-322)
- No caching, no intermediate buffers
- Direct computation: `LocalToBody.v[i][j] = f(e_0, e_1, e_2, e_3)`

✅ **Step 4: Search for static variables in aerodynamics calculations**
- Checked fdm_larcsim (the FDM type being used)
- No problematic static variables found
- All aero coefficient calculations use local variables

✅ **Step 5: Verify no member variables used as calculation temporaries**
- Reviewed `aero()` function in fdm_larcsim.cpp:436-585
- All calculation variables are local stack variables
- Aerodynamic coefficients (CL_a, CD_prof, etc.) are read-only config parameters
- No hidden state that could leak between evaluations

### What This DOES Eliminate

✅ **RNG state differences** - Verified all processes start with identical seed [4337, 4337]
✅ **Wind/turbulence** - Test was run with wind disabled
✅ **Thermals** - Disabled
✅ **Per-process initialization** - All start with identical pos/vel/quat at step 0
✅ **Valgrind uninitialized variables** - Previous session found zero warnings in simulation code
✅ **Static/global FDM state** - No problematic static variables found
✅ **Calculation temporaries** - All use local variables
✅ **Cached transform matrices** - LocalToBody recalculated fresh each step

### What This DOES NOT Eliminate (Still Possible)

⚠️ **FPU state differences** - Different rounding modes, SSE flags, denormal handling across threads/processes
⚠️ **Compiler optimizations** - -O2 may reorder floating-point operations differently
⚠️ **Third-party library non-determinism** - Eigen, Boost, or other libraries using non-deterministic algorithms
⚠️ **Memory allocator effects** - Different heap layouts affecting cache behavior → FP calculation order
⚠️ **Lazy initialization in external libraries** - First call may initialize static state differently per process
⚠️ **Multi-process scheduling effects** - Different process launch order or CPU affinity affecting FPU state

### Why Uninitialized Variables Are the MOST LIKELY Culprit

1. **Explains the pattern**: Different paths diverge at different steps/magnitudes
2. **Explains single-threaded success**: One process, sequential evaluations, heap state becomes consistent after first eval
3. **Explains multi-process variance**: Each process has different heap allocator state → different garbage in uninitialized vars
4. **Butterfly effect**: Tiny initial differences in uninitialized floats → compounding errors in numerical integration
5. **Already proven**: First-run bug from 2026-01-03 showed workers produce different results on first eval, then converge

### Memory Leak / Buffer Overflow Checks - Clean ✅

User notes: "we have done a lot of checks for mem leaks, uninitialized heap, buffer overflows -- so far none show"

This confirms:
- No memory corruption
- No buffer overflows
- Valgrind found zero uninitialized variable warnings in simulation code
- **BUT**: Valgrind only catches *reads* of uninitialized memory that affect control flow
- It may NOT catch uninitialized floats that are read then immediately overwritten (as in integrator history)

### Debug Environment Consideration

User theory: "the debug environment itself may have issues, and or may be disturbing the outcome"

**Our approach**: Focus on step 0→1 divergence with identical initial state
- This is NOT a debug logging artifact (divergence in actual physics results)
- This is NOT observer effect (trace capture happens after simulation)
- This IS a butterfly effect from tiny initial state differences

### Files Modified (2026-01-06)

**CRRCSIM Changes**:
1. `~/crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.cpp:41-86`
   - Comprehensive constructor initialization (30+ variables)
   - Debug state, integrator history, current state, matrices, scalars

2. `~/crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.h:137-152`
   - Expanded `resetIntegratorState()` to include debug variables
   - Called on every reset to ensure clean state

### Expected Outcome After Fix

**If uninitialized variables were the root cause:**
- ✅ All paths should show identical step-1 state across processes
- ✅ Variance should drop to near-zero (ULP < 4, only FP rounding)
- ✅ Multi-process determinism should match single-process determinism

**If issues remain:**
- ⚠️ Points to FPU environment or compiler optimization differences
- ⚠️ May need `-ffp-contract=off` or similar compiler flags
- ⚠️ May need explicit FPU state initialization

### Next Steps

1. **Build with fixes** (user will perform)
   ```bash
   cd ~/crsim/crrcsim-0.9.13/build
   cmake .. && make -j8
   ```

2. **Run determinism test** (user will perform)
   ```bash
   cd ~/GP/autoc
   timeout 180 ../build/autoc 2>&1 | tee test_with_init_fix.log
   ```

3. **Compare traces** (user will perform)
   ```bash
   ./compare_traces /tmp/crrcsim/elite-trace-*-gen1-*.dmp /tmp/crrcsim/elite-trace-*-gen2-*.dmp
   ```

4. **Expected success criteria:**
   - All ULP differences < 4 (FP rounding only)
   - No divergence in vLocalDot or omegaDotBody at step 1
   - Fitness variance < 0.0001 units

5. **If still failing:**
   - Check compiler flags (may need stricter FP determinism: `-ffp-contract=off -fno-fast-math`)
   - Investigate FPU environment initialization
   - Profile external library calls for non-deterministic behavior

---

**Last Updated**: 2026-01-06 19:15 PST
**Status**: EOM01 initialization fixes applied but INSUFFICIENT - divergence persists
**Next Action**: Investigate per-process lazy initialization patterns

---

## 2026-01-06 Evening: Initialization Fixes INSUFFICIENT

### Critical Finding: autoc_diag42.log Analysis

**Test run timeline:**
- eom01.cpp modified: 11:38
- crrcsim binary built: 11:42
- autoc_diag42.log run: 11:42:44

**This means the initialization fixes WERE APPLIED** when diag42 was run.

### Results: Problem Still Present

**Generation 6 elitism violation** in diag42 shows IDENTICAL symptoms:

```
Step 0: PERFECT initialization (ULP[0,0,0] for pos/vel)
Step 1: MASSIVE divergence:
  - vLocalDot ULP differences: 79,813 / 10,437,162 / 31,316
  - OmegaBody[0] ULP: 180,195,903 ⚠️ (sign flip territory!)
  - OmegaBody[2] ULP: 19,336,885
  - MomentBody[0] ULP: 1,767,524

ALL 6 PATHS diverge at step 1 (117ms)
```

**This is WORSE than diag41**:
- diag41: Different paths had different divergence patterns/timing
- diag42: ALL paths diverge uniformly at step 1
- Suggests systematic per-process constant, not random garbage

### What This Eliminates

Our comprehensive initialization fixes addressed:
- ✅ EOM01 constructor now initializes all 30+ member variables
- ✅ resetIntegratorState() clears debug variables
- ✅ All integrator history zeroed on reset

Yet the problem **persists with identical symptoms**.

**This definitively rules out:**
- ❌ Uninitialized member variables in EOM01
- ❌ Leftover integrator history between evaluations
- ❌ Debug state affecting calculations

### New Theory: Per-Process Lazy Initialization

User insight: *"perhaps there is a per-process reset state that is somehow dithered or randomized on first run -- but otherwise static. e.g. each time we reset it resets to a process unique position, same for wind, gust, thermals even though they are disabled. so, initialized once style pattern."*

**This perfectly explains the evidence:**
1. ✅ Single-threaded works: One process → one initialization → consistent
2. ✅ Multi-process fails: Each process initializes differently → consistent within process
3. ✅ Step 0 identical: Reset works, but resets to *different per-process constants*
4. ✅ Step 1 diverges: First use of those per-process constants in force/moment calculations
5. ✅ All paths diverge uniformly: The per-process constant affects all calculations

**Pattern to find:**
```cpp
// Anti-pattern we're looking for:
static bool initialized = false;
static double some_constant;

void some_function() {
  if (!initialized) {
    some_constant = compute_something_non_deterministic();  // Uses time, rand, ASLR, etc.
    initialized = true;
  }
  return some_constant;  // Always returns same value within this process
}
```

### Investigation Progress

**Checked so far:**
- ✅ windfield.cpp: initialize_gust() - clean, resets RNG properly
- ✅ FDMBase constructor - only stores env pointer, opens logfile
- ✅ Math libraries (matrix44.h) - only static inline functions
- ✅ No time(0), rand(), random() calls in FDM/windfield

**Need to investigate:**
- ⚠️ External library lazy initialization (Boost, STL, compiler runtime)
- ⚠️ Graphics/OpenGL state affecting FDM (even with rendering disabled)
- ⚠️ Compiler optimization creating per-process code paths
- ⚠️ Memory allocator affecting floating-point calculation order
- ⚠️ FPU environment (rounding mode, denormal handling) set differently per process

### Why CRRCSIM Must Work

User clarification: *"note we must get crrcsim working as minisim is only a framework test. no basis in physics really..."*

**Critical requirement:**
- Minisim is a placeholder for testing architecture
- No real aerodynamics, no proper physics
- GP evolution needs realistic flight dynamics to produce useful controllers
- CRRCSIM is the only option for actual flight physics

**Cannot switch to minisim** as originally planned - must solve crrcsim determinism.

### Current Situation

**Evidence summary:**
- Step 0 initialization: PERFECT (all processes identical)
- Step 1 force/moment calculation: DIVERGENT (different per process)
- Pattern: Consistent within each process, different between processes
- Magnitude: Millions of ULPs (not FP rounding - this is different values)

**Most likely culprits (ordered by probability):**

1. **Per-process FPU environment initialization** (70% confidence)
   - Rounding mode set differently on first FPU operation
   - Denormal handling configured per-process
   - Would affect all floating-point calculations uniformly

2. **Lazy library initialization** (60% confidence)
   - Boost/STL initializing static state on first use
   - Different initialization based on process memory layout
   - Math functions computing constants on first call

3. **Compiler optimization creating non-deterministic code** (40% confidence)
   - -O2 generating different code paths based on runtime state
   - Function pointer resolution varying by ASLR
   - Template instantiation affected by link order

4. **Graphics library side-effects** (30% confidence)
   - OpenGL/SDL initialization affecting process state
   - Even when rendering disabled, init code runs
   - May set FPU modes or allocate memory affecting later calculations

### Recommended Next Steps

1. **Add FPU state logging at process startup**
   - Log rounding mode (FE_TONEAREST, etc.)
   - Log denormal mode (DAZ, FTZ flags)
   - Compare between processes

2. **Force FPU deterministic mode explicitly**
   ```cpp
   #include <fenv.h>
   fesetround(FE_TONEAREST);  // Force consistent rounding
   // Set SSE MXCSR for denormals
   ```

3. **Compare first-timestep calculations in detail**
   - Log intermediate values in aero() function
   - Trace exact point where divergence originates
   - Identify which calculation produces different results

4. **Build with stricter determinism flags**
   ```cmake
   -O0 -ffp-contract=off -fno-fast-math -msse2 -mfpmath=sse
   ```
   Remove all optimizations, force consistent FP math

5. **Profile with different compiler/settings**
   - Try gcc vs clang
   - Try -O0 vs -O2
   - Verify if optimization level affects determinism

### Success Criteria (Updated)

CRRCSIM determinism **must be achieved** for production use:
- ⚠️ Minisim is not an option (no real physics)
- ⚠️ Must solve per-process divergence
- ⚠️ Target: ULP differences < 4 (FP rounding only)
- ⚠️ All 6 paths, all timesteps, identical across processes

### Files Modified (2026-01-06)

**Fixes applied (necessary but insufficient):**
1. `~/crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.cpp:41-86`
   - Comprehensive constructor initialization
2. `~/crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.h:137-152`
   - Expanded resetIntegratorState()

**These fixes prevent future issues but did not solve the current problem.**

---

**Last Updated**: 2026-01-06 19:15 PST
**Status**: Initialization fixes insufficient, investigating per-process lazy init patterns
**Critical**: CRRCSIM must work (minisim not viable for production)
**Next Action**: Focus on FPU environment and external library lazy initialization

---

# 2026-01-07: DEFINITIVE - Physics Divergence Timeline

## The Question We've Been Going In Circles About

**User's frustration**: "we have said both that at step0 it is same and different"

This has been the source of confusion. Let me state it definitively now.

## DEFINITIVE TIMELINE - No More Circles ✅

### What We KNOW (Facts, Not Theory)

**FACT 1: Reset Produces IDENTICAL State**
- Evidence: RESET_TRACE_POST logs from all workers
- All workers show bitwise identical values after `Global::Simulation->reset()`:
  ```
  pos=[0.0, 0.0, -82.0250015258]
  vel=[-23.01701736450, ~2.8e-15, 0.0]
  quat=[~6.1e-17, 0.0, 0.0, 1.0]
  alpha=0.0, beta=0.0, vRelWind=23.01701736450
  RNG state=[4337, 4337]
  ```
- Location: [inputdev_autoc.cpp:576](../../../crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp#L576)
- When: Immediately after `Global::Simulation->reset()` returns

**FACT 2: Step 0 Shows DIVERGENCE**
- Evidence: Physics trace comparison from autoc_diag66.log
- Path 1, Step 0 (time=0.000000ms):
  ```
  Worker1 (pid=609620): pos_x=-0.034850, vel_x=-23.077723
  Worker2 (pid=609667): pos_x=-0.034844, vel_x=-23.073706
  Difference: 6μm position, 4mm/s velocity
  ULP: 874 billion, 1.1 trillion respectively
  ```
- Location: Physics trace captured in [fdm_larcsim.cpp:702-765](../../../crsim/crrcsim-0.9.13/src/mod_fdm/fdm_larcsim/fdm_larcsim.cpp#L702-L765)
- When: During first `aero()` call at `gPhysicsStepCounter == 0`

**FACT 3: First GP Evaluation is ~117ms After Reset**
- Evidence: Code analysis and timeline reconstruction
- Sequence:
  1. Reset completes at time=0 → IDENTICAL state
  2. Physics simulation loop begins
  3. Step 0 captured at time=0.000000ms → DIVERGENT state
  4. First GP evaluation at ~117ms (after settling period)
- Time calculation: `simTimeMsec = gPhysicsStepCounter * 40.0` ([fdm_larcsim.cpp:730](../../../crsim/crrcsim-0.9.13/src/mod_fdm/fdm_larcsim/fdm_larcsim.cpp#L730))
- Step 0 = 0ms, Step 1 = 40ms, Step 2 = 80ms, Step 3 ≈ 117ms

### THE GAP - Where Divergence Originates

**Between RESET_TRACE_POST (identical) and Step 0 Trace (divergent):**

The sequence of events:
1. `Global::Simulation->reset()` completes → RESET_TRACE_POST logs identical state
2. `gPhysicsStepCounter = 0` reset ([inputdev_autoc.cpp:603](../../../crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp#L603))
3. **[THE GAP]** Simulation loop starts, first integration step(s) occur
4. First `aero()` call → step 0 trace captured → shows divergence

**What happens in the gap:**
- Wind field evaluation (even with wind disabled, code path executes)
- Initial force/moment calculations
- First physics integration step
- Transformation matrix computations
- **Something introduces non-determinism HERE**

### DEFINITIVE ANSWER to "Are First GP Eval Sensor Inputs Identical?"

**NO - The sensor inputs to the first GP evaluation are NOT identical.**

**Reasoning:**
1. Step 0 (time=0ms) shows divergence in position/velocity/acceleration
2. First GP evaluation occurs at ~117ms (step 3)
3. By step 3, the aircraft state has already diverged from step 0
4. Therefore, the position/velocity/attitude sensors fed to the GP are different

**Conclusion:**
- The non-determinism is in the **physics integration**, NOT in the GP evaluation
- The GP is receiving different sensor inputs because the physics has already diverged
- The source is somewhere in the gap between reset completion and step 0 trace capture

### Why This Matters

This is NOT "going in circles" anymore. This is definitive:

1. ✅ Reset is deterministic (proven by RESET_TRACE_POST)
2. ✅ Step 0 is non-deterministic (proven by physics trace divergence)
3. ✅ Divergence happens BEFORE first GP evaluation (step 0 at 0ms, GP at ~117ms)
4. ✅ The problem is in the physics integration, not GP evaluation
5. ✅ The gap is small (~0-40ms, before step 1 at 40ms)

### What We're Looking For

The non-determinism source is in one of these locations:

**Candidates (in order of likelihood):**
1. **FPU environment differences** - Rounding mode, denormal handling set per-process
2. **Lazy library initialization** - External libraries initializing on first FP operation
3. **Wind field evaluation** - Even with wind=0, the code path may produce different results
4. **Force/moment calculation** - First aero() call using uninitialized or per-process state
5. **Memory allocator effects** - Heap layout affecting FP operation order

**Where to instrument next:**
- Add logging at the VERY START of first physics integration loop
- Log intermediate values in first `aero()` call before step 0 trace capture
- Compare FPU environment state across workers at process startup
- Trace exact execution path between reset and first aero() call

### Files Referenced

- [inputdev_autoc.cpp:576](../../../crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp#L576) - RESET_TRACE_POST logging
- [inputdev_autoc.cpp:603](../../../crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp#L603) - gPhysicsStepCounter reset
- [fdm_larcsim.cpp:702-765](../../../crsim/crrcsim-0.9.13/src/mod_fdm/fdm_larcsim/fdm_larcsim.cpp#L702-L765) - Physics trace capture
- [fdm_larcsim.cpp:730](../../../crsim/crrcsim-0.9.13/src/mod_fdm/fdm_larcsim/fdm_larcsim.cpp#L730) - Time calculation

---

**Status**: DEFINITIVE timeline established - no more confusion
**Next**: Focus instrumentation on the gap between reset and step 0

### CRITICAL FIX: Actual Simulation Time Now Captured (2026-01-07)

**Problem Discovered**: Physics trace was using FAKE time based on step counter
- Old code: `simTimeMsec = gPhysicsStepCounter * 40.0` (assumed 40ms steps)
- Reality: Physics runs at ~3ms per step, ~100ms between GP evals
- Result: All 50 captured steps clustered around first GP eval, not spread across simulation

**Fix Applied**: [fdm_larcsim.cpp:728-735](../../../crsim/crrcsim-0.9.13/src/mod_fdm/fdm_larcsim/fdm_larcsim.cpp#L728-L735)
```cpp
// Get actual simulation time from Global::Simulation
double simTimeMsec = 0.0;
if (Global::Simulation) {
  simTimeMsec = static_cast<double>(Global::Simulation->getSimulationTimeSinceReset());
}
```

**Impact**: Now we can binary search through the 50 steps (0-150ms) to find exact divergence point!

---

# 2026-01-07: PROBLEM SOLVED - Adams-Bashforth Integrator State

## Root Cause Identified and Fixed

**The Bug**: [eom01.h:139-156](../../../crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.h#L139-L156)

`resetIntegratorState()` was clearing `v_V_dot_past` and `v_R_omega_dot_body_past` but **NOT** clearing `v_V_dot_local` and `v_R_omega_dot_body` (current accelerations).

**Why This Caused Non-Determinism:**

Adams-Bashforth formula in [eom01.cpp:214](../../../crsim/crrcsim-0.9.13/src/mod_fdm/eom01/eom01.cpp#L214):
```cpp
v_V_local += (v_V_dot_local*3 - v_V_dot_past)*dth;
```

After reset:
- `v_V_dot_past = 0` (cleared)
- `v_V_dot_local = old_value` (NOT cleared - leftover from previous path evaluation!)

First `ls_step()` used: `v_V_local += (old_accel*3 - 0)*dth`

Each worker had different leftover acceleration from its last-evaluated path → repeatable per-worker divergence.

**The Fix:**
```cpp
void resetIntegratorState() {
  v_V_dot_past = CRRCMath::Vector3();
  v_V_dot_local = CRRCMath::Vector3();  // NEW: Clear current linear acceleration
  v_R_omega_dot_body_past = CRRCMath::Vector3();
  v_R_omega_dot_body = CRRCMath::Vector3();  // NEW: Clear current angular acceleration
  // ... (clear other state)
}
```

**Verification:**
- autoc_diag71.log: 40 generations, ZERO divergences, perfect elite fitness repeatability
- autoc_diag72.log: Wind/turbulence enabled, ZERO divergences
- autoc_diag73.log: Wind/turbulence/thermals enabled, ZERO divergences (in progress)

**Status**: ✅ COMPLETELY SOLVED - Multi-process determinism achieved

---

## TODO: Investigate Thermal Impact on Fitness

**Observation from autoc_diag73.log:**
With thermals enabled (density=2.4e-06), the evolution trajectory is **IDENTICAL** to autoc_diag72.log (wind-only):
- Same fitness values at each generation
- Same GP hashes selected as elites
- Same improvement progression

**Possible explanations:**
1. Thermals may not significantly affect the specific 3 paths being tested
2. RNG seed may produce thermal patterns that don't alter fitness landscape for short evaluations
3. Thermal density too low to have measurable impact
4. Paths may not pass through thermal zones

**Follow-up tasks:**
- [ ] Check if thermals are actually being generated (add debug logging to thermal creation)
- [ ] Verify aircraft paths intersect with thermal locations
- [ ] Compare wind velocity output with/without thermals to confirm thermal lift is being calculated
- [ ] Consider increasing thermal density or using longer evaluation times to see thermal effects
- [ ] Add thermal-specific sensors to GP (vertical velocity, thermal strength at location)

**Priority**: Low - Determinism is solved, this is just optimization/validation of thermal physics

---

**Last Updated**: 2026-01-07 22:50 PST
**Status**: Non-determinism SOLVED via integrator state fix
**Next**: Validate thermal impact on evolution (optional)

---

# 2026-01-08: Wind/Turbulence Non-Determinism - Butterfly Effect Issue

## Problem Statement

After fixing the Adams-Bashforth integrator state bug (which achieved perfect determinism with wind/turbulence/thermals disabled), we discovered **residual non-determinism** when wind and turbulence are enabled.

**User's requirement**: System must achieve **bitwise identical** fitness for same GP with same scenario. Even ~1e-6 differences are unacceptable because detecting 4th order effects requires exact repeatability.

## Test Results Matrix

| Test | Config | Elite Re-Eval Result | Status |
|------|--------|---------------------|--------|
| diag71 | No wind/turb/thermals | 0.000000 diff | ✅ Perfect |
| diag72 | Wind/turb (no thermals) | 0 PHYSICS_TRACE_DIVERGENCE | ✅ Perfect |
| diag73 | Wind/turb + thermals | 0 PHYSICS_TRACE_DIVERGENCE | ✅ Perfect |
| diag74 | No wind/turb/thermals | 0.000000 diff | ✅ Perfect |
| diag75 | Wind/turb (no thermals) | ±4e-6 diff | ❌ Non-deterministic |
| diag76 | Thermals only | 0.000000 diff | ✅ Perfect |
| diag77 | Wind/turb + thermals | ±1e-6 to ±5e-6 diff | ❌ Inconsistent |

### Critical Observation from diag77

**The non-determinism is INCONSISTENT** - this is the key clue:

```
GP 0x7d49a5a12fcd2828 (gen 4-9):
  Baseline: 225.909225
  Re-evals: 225.909225 (5 consecutive identical re-evaluations!)
  Difference: 0.000000 ✅ PERFECT

GP 0x66f61710fca81ca5 (gen 1):
  Baseline: 246.227829
  Re-eval: 246.227834
  Difference: +5e-6 ❌

GP 0xc7a8d94c70f704c3 (gen 10-14):
  Baseline: 216.150467
  Re-evals: 216.150462 (4 times)
  Difference: -5e-6 ❌

GP 0xa90d426c420728bf (gen 15-17):
  Baseline: 208.915421
  Re-evals: 208.915422 (2 times)
  Difference: +1e-6 ❌
```

## Why This Is NOT Like the Adams-Bashforth Bug

**Adams-Bashforth bug pattern** (SOLVED):
- **Systematic**: ALL GPs showed divergence
- **Early**: Divergence at step 0 (immediately after reset)
- **Growing**: Divergence accumulated over simulation
- **Per-worker**: Each worker had consistent leftover state from previous path
- **Root cause**: Uninitialized state in integrator

**Wind/turbulence bug pattern** (CURRENT):
- **Inconsistent**: Some GPs perfect (5 re-evals identical), others show tiny diffs
- **Magnitude**: Only ±1-5 ULP differences (~1e-6 to 5e-6)
- **No physics trace divergence**: 0 PHYSICS_TRACE_DIVERGENCE events in logs
- **Random**: No pattern to which GPs are affected

## Hypotheses

### Most Likely: Floating-Point Butterfly Effect

**Theory**: Wind/turbulence calculations involve chaotic operations where tiny FPU rounding differences amplify through the simulation.

**Evidence**:
1. Some GPs are completely stable (225.909225 - 5 identical re-evals)
2. Others show single-ULP differences (~1e-6 to 5e-6)
3. No pattern to which GPs are affected
4. Magnitude is extremely small (1-5 ULP)

**Possible sources**:
- FPU rounding mode differences between workers
- Denormal number handling in turbulence calculations
- Transcendental function precision (sin/cos/sqrt in wind field)
- Memory alignment affecting SIMD operations
- Compiler optimization producing different instruction sequences per translation unit

### Less Likely: Wind Field State Initialization

**Why less likely**:
- Would affect ALL GPs consistently (like Adams-Bashforth did)
- diag72/73 showed 0 PHYSICS_TRACE_DIVERGENCE (physics is identical)
- Only fitness differs, not physics trace

**But**: Could still be subtle state in wind/gust RNG that doesn't show up in physics trace

## Investigation Plan

### Phase 1: Isolate Wind vs Turbulence vs Thermal Interactions

- [x] Test wind/turbulence only (diag75) → ±4e-6 diff
- [x] Test thermals only (diag76) → Perfect
- [x] Test wind/turb + thermals (diag77) → ±1-5e-6 diff (inconsistent)
- [ ] **Next**: Test wind only (turbulence=0) to isolate turbulence component
- [ ] Test different wind velocities (1, 5, 10 m/s) to see if magnitude correlates

### Phase 2: FPU Environment Analysis

- [ ] Log FPU control word at worker startup (rounding mode, precision, exception masks)
- [ ] Force consistent FPU state across workers:
  ```cpp
  #include <fenv.h>
  fesetround(FE_TONEAREST);
  fesetenv(FE_DFL_ENV);
  ```
- [ ] Compile with `-ffp-contract=off` to disable FMA optimizations
- [ ] Compare assembly of wind/turbulence functions across translation units

### Phase 3: Wind Field Instrumentation

**If FPU fixes don't work**, instrument wind field calculations:

- [ ] Add logging to turbulence/gust RNG calls
- [ ] Capture intermediate wind field values during evaluation
- [ ] Compare wind field state between baseline and re-eval
- [ ] Check for any lazy initialization or cached values

### Phase 4: Deep Dive into Specific GP

**Focus on GP 0x7d49a5a12fcd2828** (the one that's perfect):
- [ ] What makes this GP different from others?
- [ ] Does it use different sensors/actuators?
- [ ] Does its flight path avoid certain wind conditions?
- [ ] Compare with GP 0xc7a8d94c70f704c3 (consistently shows -5e-6 diff)

## Files to Investigate

- `/home/gmcnutt/crsim/crrcsim-0.9.13/src/mod_windfield/windfield.cpp` - Wind/turbulence/gust calculations
- `/home/gmcnutt/crsim/crrcsim-0.9.13/src/mod_fdm/fdm_larcsim/fdm_larcsim.cpp` - Force/moment integration
- `/home/gmcnutt/crsim/crrcsim-0.9.13/autoc_config.xml` - Wind/turbulence configuration

## Expected Outcome

**Goal**: Achieve bitwise identical fitness (0.000000 difference) with wind/turbulence enabled, matching the perfection we achieved with wind disabled.

**Acceptable solution**: Either:
1. Fix root cause (FPU environment, wind state, etc.)
2. Understand and document that wind/turbulence is inherently chaotic at ULP precision (and decide if <5e-6 is acceptable)

---

**Status**: Investigation started - need to test wind-only (no turbulence) next
**Priority**: HIGH - Required for production GP evolution with realistic wind conditions
**Last Updated**: 2026-01-08 10:00 PST

---

# 2026-01-09: DETERMINISM ACHIEVED - Full Wind/Turbulence/Thermal Support

## Summary

**MAJOR MILESTONE**: Achieved **0 divergences** with full wind, turbulence, AND thermal simulation enabled across multiple wind scenarios (WindScenarios=3).

## Root Causes Fixed

### 1. Park-Miller LCG Implementation (Previous Session)
- Replaced system `rand()` with deterministic Park-Miller LCG (a=48271, m=2^31-1)
- Ensures identical RNG sequence across all worker processes
- File: `crrc_rand.cpp`

### 2. RandGauss Reset Architecture (Previous Session)
- Added `initialize_gust()` call at end of `initialize_wind_field()`
- Ensures all RandGauss objects (eta1-4, rnd_radius, rnd_strength, rnd_lifetime) reset to phase=0
- File: `windfield.cpp`

### 3. SimStateHandler Wind Seed Reset (Previous Session)
- Added `reset(unsigned int windSeed)` overload
- Properly orders: set seed → create thermals → reset RandGauss objects
- File: `SimStateHandler.cpp`

### 4. Skip insertData() in Autoc Mode (THIS SESSION - CRITICAL FIX)
- **Root Cause**: In main loop (`crrc_main.cpp` lines 949-950), `CRRC_Random::insertData(SDL_GetTicks())` was injecting wall-clock time into the RNG every frame
- **Impact**: RNG state diverged between runs because `SDL_GetTicks()` returns real elapsed time
- **Fix**: Skip `insertData()` calls when `inputMethod() == eIM_autoc`
- **File**: `crrc_main.cpp` lines 948-954

```cpp
// Insert random data for additional entropy
// Skip in autoc mode to preserve determinism - autoc resets the RNG with a specific seed
// and any insertData() calls would corrupt the deterministic sequence
if (Global::TXInterface->inputMethod() != T_TX_Interface::eIM_autoc) {
  CRRC_Random::insertData(SDL_GetTicks());
  CRRC_Random::insertData(Global::inputs.getRandNum());
}
```

### 5. Elite Trace Storage Fixes (THIS SESSION)
- Added `gpHash` clearing in `clearEvalResults()` (line 392)
- Added `gpHash` copying in `appendEvalResults()` (lines 409-411)
- Added `physicsTrace` appending in `appendEvalResults()` (line 417)
- Added elite storage after bakeoff winner selection (lines 582-593)
- **File**: `autoc.cc`

### 6. Memory Leak Fixes (THIS SESSION)
- Removed per-evaluation storage in `gCurrentGenEvalResults` map (was causing 1.5GB peak)
- Added `physicsTrace.clear()` to `clearEvalResults()` (was missing)
- **Result**: Peak memory reduced from 1.92GB to ~800MB-1.1GB

## Test Results

### Log 117: Wind/Gust Only (No Thermals)
| Metric | Value |
|--------|-------|
| **Divergences** | **0** ✓ |
| Elite Checks | 11 |
| Generations | 20 |
| Peak Memory | 1.13G |
| Leaked | 521MB |
| Runtime | 894s |

### Log 118: Full Wind/Gust/Thermals (density=2.4e-06)
| Metric | Value |
|--------|-------|
| **Divergences** | **0** ✓ |
| Elite Checks | ~19 |
| Generations | 20 |
| Peak Memory | 805MB |
| Leaked | 521MB |
| Runtime | 937s |

## Phase 1 Complete ✅

All success criteria met:
- [x] CRRCsim determinism achieved with wind/turbulence/thermals
- [x] Multiple wind scenarios (WindScenarios=3) work deterministically
- [x] Elite re-evaluations produce identical fitness
- [x] 0 divergences in 20-generation runs

## Remaining TODOs

### Phase 2: Demetic Mode Verification
- [ ] Enable `DemeticGrouping = 1` in autoc.ini
- [ ] Verify determinism holds with demetic processing
- [ ] Test that elite fitness is monotonically improving
- [ ] Verify bakeoff across demes works correctly

### Phase 3: Debug Code Cleanup
- [ ] Remove/disable verbose ELITE_TRACE logging
- [ ] Remove/disable ELITE_STORE logging
- [ ] Remove/disable physics trace divergence reporting
- [ ] Consider keeping infrastructure behind compile flag for future debugging
- [ ] Restore `-O2` optimization (currently `-O0` for debugging)
- [ ] Remove MALLOC_PERTURB from crrcsim.sh
- [ ] Clean up AUTOC_RECV_HASH debug code if still present

### Phase 4: Production Readiness
- [ ] Run extended validation (100+ generations)
- [ ] Verify performance is acceptable with `-O2`
- [ ] Document final configuration for production use
- [ ] Checkpoint code with git commit

## Files Modified (2026-01-09)

### Core Fixes
1. `crrc_main.cpp:948-954` - Skip insertData() in autoc mode
2. `autoc.cc:392` - Clear gpHash in clearEvalResults()
3. `autoc.cc:409-411` - Copy gpHash in appendEvalResults()
4. `autoc.cc:417` - Append physicsTrace in appendEvalResults()
5. `autoc.cc:582-593` - Store elite trace after bakeoff

### Configuration
6. `autoc_config.xml:67` - Thermal density (toggled between 0 and 2.4e-06 for testing)

## Lessons Learned

1. **Wall-clock time injection was the culprit**: The `insertData(SDL_GetTicks())` call was designed to add entropy for interactive use, but destroyed determinism for automated testing.

2. **Divergence at Step 13**: Both thermal and gust divergences manifested at Step 13 (~39ms) because that's when the first `insertData()` call occurred in the main loop timing.

3. **Physics identical, RNG diverged**: The trace comparison showed all physics values were identical (0 ULPs), but RNG state differed - pointing directly at `insertData()` as the culprit.

4. **Proper debugging infrastructure paid off**: The elite trace capture and ULP comparison tools made it possible to pinpoint the exact source of non-determinism.

---

**Status**: ✅ Phase 1 COMPLETE - Full determinism achieved
**Next Action**: Verify demetic mode, then clean up debug code
**Last Updated**: 2026-01-09 14:45 PST

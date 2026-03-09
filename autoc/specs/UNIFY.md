# Unifying GP and Bytecode Evaluation Code Paths

## Executive Summary

The autoc.cc codebase contains two parallel evaluation paths for computing fitness:
1. **MyGP::evalTask()** - Used during normal GP evolution (lines 691-1134)
2. **BytecodeEvaluationGP::evalTask()** - Used for bytecode verification mode (lines 1320-1654)

These paths duplicate approximately 350 lines of nearly identical fitness computation logic, creating a maintenance burden where changes must be applied twice and can easily diverge.

This document proposes a unified architecture that extracts the common evaluation logic into shared components while preserving the flexibility needed for each mode's unique requirements.

---

## UPDATE LOG (Jan 2026)

**Fitness Function v2 Changes:**
The fitness function has been significantly simplified since this document was written. Key changes:

1. **Metrics reduced from 4 to 3:**
   - Waypoint distance (meters, raw) with **asymmetric power** (ahead/behind rabbit)
   - Movement direction alignment (scaled to meter-equivalent)
   - **Energy deviation** (J/kg) replaces throttle efficiency - uses specific energy (0.5*v² + g*h)
   - Cross-track error removed (subsumed by waypoint distance)

2. **Raw units, no percentage scaling:**
   - Old: Everything normalized to 0-100% then power-weighted
   - New: Raw meters and J/kg with direct power-law weighting

3. **Simpler aggregation:**
   - Old: Time-weighted average across paths
   - New: Simple sum with per-step time weight (`STEP_TIME_WEIGHT`)

4. **Soft lexicographic crash penalty:**
   - Old: Complex power-law scaling based on fraction completed
   - New: Additive penalty `CRASH_COMPLETION_WEIGHT * (1 - fraction_completed)`

5. **Demetic mode shelved:**
   - All demetic branching removed from the codebase
   - Phase 5 (ScenarioStrategy abstraction) is no longer needed

**Impact on Unification:**
- **Simpler extraction** - fewer metrics, cleaner code
- **Reduced risk** - less branching logic to preserve
- **Phase 5 eliminated** - no deme mode to abstract

---

---

## 1. Current State Analysis

### 1.1 Duplicated Code Locations

| Component | MyGP::evalTask() | BytecodeEvaluationGP::evalTask() | Notes |
|-----------|------------------|----------------------------------|-------|
| Scenario metadata building | Lines 709-781 | Lines 1328-1389 | Nearly identical |
| RPC send/receive | Lines 790-822 | Lines 1399-1417 | Identical |
| Path loop initialization | Lines 855-898 | Lines 1430-1457 | Identical |
| Fitness metric accumulators | Lines 893-898 | Lines 1452-1457 | Identical |
| Step loop processing | Lines 900-984 | Lines 1459-1532 | Identical |
| Cross-track calculation | Lines 920-931 | Lines 1477-1485 | Identical |
| Movement direction | Lines 933-950 | Lines 1487-1498 | Identical |
| Throttle efficiency | Lines 952-959 | Lines 1501-1506 | Identical |
| Fitness normalization | Lines 1081-1101 | Lines 1615-1632 | Identical |
| Crash penalty | Lines 1107-1122 | Lines 1638-1646 | Identical |

### 1.2 Key Differences Between Paths

| Aspect | MyGP::evalTask() | BytecodeEvaluationGP::evalTask() |
|--------|------------------|----------------------------------|
| GP serialization | Uses save() to serialize GP tree | Pre-serialized bytecode buffer |
| Hash computation | hashString(lispForm) | hashByteVector(bytecodeBuffer) |
| Elite reeval detection | Checks gLastEliteGpHash | Not implemented |
| Physics trace capture | Stores to gPendingEliteEvalResults | Not implemented |
| Detailed logging (printEval) | Full logging with header | Simpler format, different header |
| bakeoffSequence handling | Uses bakeoffPathCounter | Uses bakeoffPathCounter |

### 1.3 Data Flow Through Each Path

**MyGP::evalTask() Flow:**
```
1. Serialize GP tree to buffer
2. Get scenario from scenarioForIndex()
3. Build EvalData with:
   - gp buffer
   - gpHash (from Lisp string)
   - pathList from scenario
   - scenarioList with metadata
4. Set isEliteReeval if matching gLastEliteGpHash
5. Send EvalData via RPC
6. Receive EvalResults
7. For each path in results:
   - Walk timesteps computing metrics
   - Accumulate waypoint distance, cross-track, movement direction, throttle
   - Apply normalization by path odometer
   - Apply crash penalty if crashed
   - Sum to stdFitness
8. Store physics traces if elite reeval
9. Set fitnessValid = 1
```

**BytecodeEvaluationGP::evalTask() Flow:**
```
1. Get scenario from scenarioForIndex()
2. Build EvalData with:
   - bytecodeBuffer (pre-serialized)
   - gpHash (from bytecode hash)
   - pathList from scenario
   - scenarioList with metadata
3. Send EvalData via RPC
4. Receive EvalResults
5. For each path in results:
   - Walk timesteps computing metrics (SAME as above)
   - Apply normalization (SAME as above)
   - Apply crash penalty (SAME as above)
   - Sum to stdFitness
6. Set fitnessValid = 1
```

### 1.4 Scenario Metadata Building (Duplicated Pattern)

Both paths contain nearly identical logic for building scenarioList (approximately 50 lines each).

**UPDATE (Jan 2026):** The demetic mode code has been largely removed/simplified. Scenario metadata building is now simpler with `pathVariantIndex` defaulting to -1 (unset/aggregated). The complexity around demetic vs non-demetic branching in metadata construction is no longer a primary concern.

---

## 2. Proposed Architecture

### 2.1 Core Abstraction: FitnessComputer

Extract the fitness computation loop into a standalone function that operates on evaluation results:

```cpp
// fitness_computer.h
#ifndef FITNESS_COMPUTER_H
#define FITNESS_COMPUTER_H

#include "minisim.h"

// Fitness computation result for a single path
struct PathFitnessResult {
    gp_fitness fitness;           // Computed fitness for this path
    int simulationSteps;          // Number of steps processed
    bool hadNaN;                   // Whether NaN was detected

    // Per-metric breakdown (for debugging/analysis) - UPDATED for new fitness function
    gp_fitness waypointComponent;     // Now: raw meters with asymmetric power
    gp_fitness movementComponent;     // Direction alignment scaled to meters
    gp_fitness energyComponent;       // NEW: replaces throttle (J/kg deviation)
    // NOTE: crossTrackComponent removed - now subsumed by waypoint distance
    // NOTE: throttleComponent removed - replaced by energy-based metric
};

// Compute fitness for a single path
// This is the core shared logic extracted from evalTask()
// UPDATED: Now uses raw units, asymmetric power-law weighting, and energy-based altitude metric
PathFitnessResult computePathFitness(
    const std::vector<Path>& path,
    const std::vector<AircraftState>& aircraftStates,
    CrashReason crashReason
);

// Aggregate fitness across all paths in an evaluation result
// UPDATED: Simple sum (no time-weighted average), crash penalty is per-path additive
gp_fitness computeTotalFitness(const EvalResults& results);

// Helper: compute path frame at given index (used in fitness computation)
// Already exists as computePathFrame() in anonymous namespace
PathFrame computePathFrame(const std::vector<Path>& path, int index);

#endif
```

### 2.2 Core Abstraction: EvalDataBuilder

Extract scenario metadata construction into a builder:

```cpp
// eval_data_builder.h
#ifndef EVAL_DATA_BUILDER_H
#define EVAL_DATA_BUILDER_H

#include "minisim.h"
#include "autoc.h"

class EvalDataBuilder {
public:
    EvalDataBuilder(const ScenarioDescriptor& scenario, bool isBakeoff);

    // Set the GP payload (either serialized tree or bytecode)
    EvalDataBuilder& setGpPayload(const std::vector<char>& buffer, uint64_t hash);

    // Enable elite reeval tracing
    EvalDataBuilder& setEliteReeval(bool enabled);

    // Enable deterministic logging
    EvalDataBuilder& setDeterministicLogging(bool enabled);

    // Build the EvalData struct
    EvalData build();

private:
    const ScenarioDescriptor& scenario_;
    bool isBakeoff_;
    std::vector<char> gpBuffer_;
    uint64_t gpHash_;
    bool isEliteReeval_ = false;
    bool enableLogging_ = false;

    void buildScenarioList(EvalData& data);
};
```

### 2.3 Unified Evaluation Flow

```cpp
// In both evalTask() implementations:
void evalTaskCommon(WorkerContext& context,
                    const std::vector<char>& gpBuffer,
                    uint64_t gpHash,
                    int scenarioIndex,
                    bool isBakeoff,
                    bool checkEliteReeval) {
    const ScenarioDescriptor& scenario = scenarioForIndex(scenarioIndex);

    // Build evaluation data using shared builder
    EvalDataBuilder builder(scenario, isBakeoff);
    builder.setGpPayload(gpBuffer, gpHash);

    if (checkEliteReeval && gpHash != 0 && gpHash == gLastEliteGpHash) {
        builder.setEliteReeval(true);
    }

    EvalData evalData = builder.build();
    evalData.sanitizePaths();

    // RPC exchange
    sendRPC(*context.socket, evalData);
    context.evalResults = receiveRPC<EvalResults>(*context.socket);

    // Backfill missing data (shared logic)
    backfillEvalResults(context.evalResults, evalData);

    // Compute fitness using shared function
    gp_fitness totalFitness = computeTotalFitness(context.evalResults);

    return totalFitness;
}
```

### 2.4 Detailed Logging Extraction

The printEval logging block is also duplicated but has minor formatting differences. Extract into a configurable logger:

```cpp
// eval_logger.h
class EvalLogger {
public:
    enum Format { Full, Compact };

    EvalLogger(std::ostream& output, Format format = Full);

    void printHeader(int pathVariant, int windVariant);
    void logStep(const LogStepData& data);

private:
    std::ostream& out_;
    Format format_;
    bool headerPrinted_ = false;
};
```

---

## 3. Deme vs Non-Deme Processing Unification

**UPDATE (Jan 2026):** Demetic mode has been shelved. The codebase now operates exclusively in non-demetic mode, significantly simplifying this section. The `pathVariantIndex` now defaults to -1 (unset/aggregated) and the demetic branching logic has been removed.

### 3.1 Current State (Simplified)

With demetic mode removed:
- **Scenario building**: Creates ONE scenario containing ALL path variants × ALL winds
- **Scenario assignment**: Always scenario 0
- **Bakeoff**: Only evaluates bestOfPopulation across all scenarios

### 3.2 Implications for Unification

The removal of demetic mode **simplifies the unification effort considerably**:
- No more conditional branching on `DemeticGrouping`
- Scenario metadata building is now straightforward
- The `ScenarioStrategy` abstraction is **no longer needed** for Phase 5

This reduces the refactoring scope and risk.

---

## 4. Specific Refactoring Steps

### Phase 1: Extract Fitness Computation (Low Risk)

1. **Create fitness_computer.h/cc**
   - Extract `computePathFrame()` from anonymous namespace
   - Create `computePathFitness()` for single-path fitness
   - Create `computeTotalFitness()` for aggregating path results

2. **Update MyGP::evalTask()**
   - Replace inline fitness computation with `computeTotalFitness()`
   - Keep all RPC and logging logic in place

3. **Update BytecodeEvaluationGP::evalTask()**
   - Replace inline fitness computation with same function

4. **Test**: Verify identical fitness values for same inputs

### Phase 2: Extract Scenario Metadata Building (Medium Risk)

1. **Create eval_data_builder.h/cc**
   - Implement `EvalDataBuilder` class
   - Encapsulate demetic vs non-demetic metadata logic

2. **Update both evalTask() implementations**
   - Replace manual EvalData construction with builder

3. **Test**: Verify metadata identical for same scenarios

### Phase 3: Extract Logging (Low Risk)

1. **Create eval_logger.h/cc**
   - Implement `EvalLogger` with format options
   - Unify header printing logic

2. **Update evalTask() logging blocks**
   - Replace sprintf/fout with EvalLogger calls

3. **Test**: Compare output files for identical content

### Phase 4: Unify Evaluation Core (Higher Risk)

1. **Create common evaluation helper**
   - Single function that both evalTask() implementations call
   - Parameterized by GP buffer, hash, and feature flags

2. **Reduce BytecodeEvaluationGP**
   - Minimize to just bytecode buffer handling
   - Inherit common evaluation logic

3. **Test**: Full regression testing with multiple scenarios

### Phase 5: Strategy Pattern for Deme Modes (Optional)

**UPDATE (Jan 2026):** This phase is **no longer needed**. Demetic mode has been shelved, so there is no need to abstract the deme/non-deme branching. Skip this phase.

---

## 5. Risk Assessment

### 5.1 High-Risk Areas

| Area | Risk | Mitigation |
|------|------|------------|
| Fitness computation extraction | Medium | Extensive unit tests comparing old vs new outputs |
| Floating-point precision | High | Use gp_fitness (double) consistently, compare bit-exact |
| Crash penalty calculation | Medium | Test with various crash points (0%, 50%, 99%). Note: penalty is now additive (soft lexicographic) |
| Deme/non-deme divergence | ~~High~~ Low | Demetic mode removed - no longer a concern |
| printEval formatting | Low | Diff output files before/after |

### 5.2 Testing Strategy

1. **Unit Tests**
   - `computePathFitness()` with known inputs/outputs
   - `EvalDataBuilder` produces correct metadata
   - Crash penalty edge cases

2. **Integration Tests**
   - Run identical GP through both old and new paths
   - Compare fitness values bit-for-bit
   - Compare output file contents

3. **Regression Tests**
   - Full GP evolution run with identical seeds
   - Verify generation statistics match
   - Verify best-of-generation fitness matches

### 5.3 Rollback Plan

- Keep original code in `#ifdef UNIFIED_EVAL` blocks initially
- Run both paths in parallel during transition
- Compare results before committing to unified path

---

## 6. Migration Strategy

### 6.1 Phased Rollout

| Phase | Changes | Testing | Notes |
|-------|---------|---------|-------|
| 1 | Extract fitness_computer | Unit + integration | **Simpler now** - fewer metrics (3 vs 4), no time-weighted average |
| 2 | Extract eval_data_builder | Unit + integration | **Simpler** - demetic branching removed |
| 3 | Extract eval_logger | Output comparison | Log format simplified (fewer columns) |
| 4 | Unify evalTask core | Full regression | Primary goal unchanged |
| 5 | ~~Strategy pattern~~ | ~~Mode-specific tests~~ | **Skipped** - demetic mode shelved |

### 6.2 Compatibility Considerations

- EvalData/EvalResults serialization version unchanged
- S3 output format unchanged
- data.dat/data.stc format unchanged
- Worker RPC protocol unchanged

### 6.3 Documentation Updates

After unification:
- Update this document with final architecture
- Add API documentation for new modules
- Update ARCHITECTURE.md (if exists) with new component diagram

---

## 7. Code Examples

### 7.1 Extracted computePathFitness()

**UPDATE (Jan 2026):** The fitness function has been significantly simplified. The new architecture uses:

1. **Raw units instead of percentages** - Waypoint distance in meters, energy deviation in J/kg
2. **Asymmetric power-law weighting** - Position and altitude errors use different exponents based on direction
3. **Per-step accumulation with time weighting** - `step_fitness_sum += step_fitness * STEP_TIME_WEIGHT`
4. **Energy-based altitude metric** - Replaces throttle efficiency with specific energy deviation (0.5*v² + g*h)
5. **Soft lexicographic crash penalty** - Uses `CRASH_COMPLETION_WEIGHT * (1 - fraction_completed)` additive penalty

The new fitness metrics are:
- **Waypoint distance (meters)**: Asymmetric power based on ahead/behind rabbit (`WAYPOINT_AHEAD_POWER=1.8`, `WAYPOINT_BEHIND_POWER=1.2`)
- **Movement direction alignment**: Scaled to ~meter equivalent via `DIRECTION_SCALE`
- **Energy deviation (J/kg)**: Asymmetric by altitude (`ALTITUDE_LOW_POWER=1.5`, `ALTITUDE_HIGH_POWER=1.0`)

```cpp
PathFitnessResult computePathFitness(
    const std::vector<Path>& path,
    const std::vector<AircraftState>& aircraftStates,
    CrashReason crashReason) {

    PathFitnessResult result{};

    if (path.empty() || aircraftStates.empty()) {
        result.fitness = 0.0;
        return result;
    }

    gp_fitness step_fitness_sum = 0.0;
    int stepIndex = 0;

    while (++stepIndex < static_cast<int>(aircraftStates.size())) {
        const auto& stepAircraftState = aircraftStates.at(stepIndex);
        int pathIndex = std::clamp(stepAircraftState.getThisPathIndex(), 0,
                                   static_cast<int>(path.size()) - 1);
        const Path& currentPathPoint = path.at(pathIndex);
        PathFrame frame = computePathFrame(path, pathIndex);

        gp_vec3 aircraftPosition = stepAircraftState.getPosition();

        // Metric 1: Waypoint distance (meters, asymmetric power)
        gp_vec3 craftOffset = aircraftPosition - currentPathPoint.start;
        gp_scalar waypointDistance = craftOffset.norm();
        gp_scalar alongTrack = craftOffset.dot(frame.tangent);
        gp_scalar waypoint_power = (alongTrack > 0)
            ? WAYPOINT_AHEAD_POWER : WAYPOINT_BEHIND_POWER;

        // Metric 2: Movement direction alignment
        gp_scalar movementDirectionError = 0.0f;
        if (stepIndex > 1) {
            gp_vec3 movement = stepAircraftState.getPosition() -
                               aircraftStates.at(stepIndex-1).getPosition();
            if (movement.norm() > 0.1f) {
                gp_scalar alignment = std::clamp(movement.normalized().dot(frame.tangent),
                                                 -1.0f, 1.0f);
                movementDirectionError = (1.0f - alignment) * DIRECTION_SCALE;
            }
        }

        // Metric 3: Energy deviation (J/kg = m^2/s^2)
        gp_scalar craft_speed = stepAircraftState.getVelocity().norm();
        gp_scalar craft_altitude = -aircraftPosition.z();  // NED: -z = altitude
        gp_scalar craft_energy = 0.5f * craft_speed * craft_speed + 9.81f * craft_altitude;

        gp_scalar rabbit_speed = SIM_RABBIT_VELOCITY;
        gp_scalar rabbit_altitude = -currentPathPoint.start.z();
        gp_scalar rabbit_energy = 0.5f * rabbit_speed * rabbit_speed + 9.81f * rabbit_altitude;

        gp_scalar energy_deviation = std::abs(craft_energy - rabbit_energy);
        gp_scalar altitude_power = (craft_altitude < rabbit_altitude)
            ? ALTITUDE_LOW_POWER : ALTITUDE_HIGH_POWER;

        // Per-step fitness with time weight
        gp_fitness step_fitness = pow(waypointDistance, waypoint_power) +
                                  pow(movementDirectionError, MOVEMENT_DIRECTION_WEIGHT) +
                                  pow(energy_deviation, altitude_power);
        step_fitness_sum += step_fitness * STEP_TIME_WEIGHT;

        result.simulationSteps++;
    }

    result.fitness = step_fitness_sum;

    // Crash penalty: soft lexicographic
    if (crashReason != CrashReason::None) {
        gp_fitness total_path_distance = path.back().distanceFromStart;
        gp_fitness fraction_completed = std::max(
            path.at(aircraftStates.back().getThisPathIndex()).distanceFromStart / total_path_distance,
            0.001);
        result.fitness = CRASH_COMPLETION_WEIGHT * (1.0 - fraction_completed) + result.fitness;
    }

    result.hadNaN = std::isnan(result.fitness);
    return result;
}
```

### 7.2 Extracted EvalDataBuilder

```cpp
EvalData EvalDataBuilder::build() {
    EvalData data;
    data.gp = gpBuffer_;
    data.gpHash = gpHash_;
    data.isEliteReeval = isEliteReeval_;
    data.pathList = scenario_.pathList;

    uint64_t scenarioSeq = globalScenarioCounter.fetch_add(1, std::memory_order_relaxed) + 1;

    buildScenarioList(data);

    if (!data.scenarioList.empty()) {
        data.scenario = data.scenarioList.front();
    } else {
        data.scenario.scenarioSequence = scenarioSeq;
    }

    return data;
}

void EvalDataBuilder::buildScenarioList(EvalData& data) {
    const GPVariables& gpCfg = ConfigManager::getGPConfig();
    uint64_t scenarioSeq = globalScenarioCounter.load(std::memory_order_relaxed);

    data.scenarioList.reserve(scenario_.pathList.size());

    for (size_t idx = 0; idx < scenario_.pathList.size(); ++idx) {
        ScenarioMetadata meta;
        meta.scenarioSequence = scenarioSeq;
        meta.enableDeterministicLogging = enableLogging_;

        if (isBakeoff_) {
            meta.bakeoffSequence = bakeoffPathCounter.fetch_add(1, std::memory_order_relaxed) + 1;
        }

        // Delegate to strategy (simplified inline version shown)
        if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0) {
            meta.pathVariantIndex = scenario_.pathVariantIndex;
            if (idx < scenario_.windScenarios.size()) {
                meta.windVariantIndex = scenario_.windScenarios[idx].windVariantIndex;
                meta.windSeed = scenario_.windScenarios[idx].windSeed;
            }
        } else {
            size_t numWinds = scenario_.windScenarios.size();
            if (numWinds > 0) {
                meta.pathVariantIndex = static_cast<int>(idx / numWinds);
                size_t windIdx = idx % numWinds;
                meta.windVariantIndex = scenario_.windScenarios[windIdx].windVariantIndex;
                meta.windSeed = scenario_.windScenarios[windIdx].windSeed;
            }
        }

        data.scenarioList.push_back(meta);
    }
}
```

---

## 8. Conclusion

The proposed unification reduces code duplication from ~350 duplicated lines to a single shared implementation. Key benefits:

1. **Maintainability**: Changes to fitness computation apply once
2. **Consistency**: No risk of paths diverging
3. **Testability**: Extracted functions are unit-testable
4. **Extensibility**: New fitness metrics require one change point

The phased approach minimizes risk by allowing verification at each step before proceeding.

### Recommended First Step

Begin with Phase 1 (fitness_computer extraction) as it:
- Has the highest duplication impact (~100 lines)
- Is lowest risk (pure computation, no side effects)
- Provides immediate testing benefit
- Does not require changes to RPC or logging

---

## Appendix A: Line Count Summary

| Component | MyGP Lines | BytecodeGP Lines | Duplicated |
|-----------|------------|------------------|------------|
| Scenario metadata | 72 | 62 | ~60 |
| RPC handling | 32 | 28 | ~25 |
| Fitness loop | 180 | 170 | ~150 |
| Normalization/crash | 45 | 40 | ~40 (simplified) |
| Logging | 70 | 80 | ~50 |
| **Total** | **~400** | **~380** | **~325** |

## Appendix B: Files to Create

```
/home/gmcnutt/GP/autoc/
  fitness_computer.h    # Core fitness computation (updated for new metrics)
  fitness_computer.cc
  eval_data_builder.h   # EvalData construction (simplified, no demetic branching)
  eval_data_builder.cc
  eval_logger.h         # Unified logging (simplified output format)
  eval_logger.cc
  # scenario_strategy.h/cc - NO LONGER NEEDED (demetic mode shelved)
```

## Appendix C: Related Specifications

- specs/LARGE_SCENARIO_STRATEGY.md - Deme mode architecture
- specs/PARETO.md - Multi-objective fitness (future extension point)
- specs/FASTMATH.md - Numerical precision considerations

## Appendix D: Bug Fixes (Mar 2026)

### D.1 Temporal History Wipe (inputdev_autoc.cpp)

**Problem:** `aircraftState` was reconstructed from scratch every eval tick (line 747),
destroying the temporal history ring buffer. All `GETDPHI_PREV(n)` for n>=1,
`GETDTHETA_PREV(n)`, `GETDPHI_RATE`, and `GETDTHETA_RATE` nodes returned 0.0 or
garbage during training — making them dead weight in the GP search space.

**Fix:** Replace aggregate initialization with in-place setter calls that preserve
the history fields. Add `clearHistory()` call on path reset.

**Files changed:**
- `crrcsim/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp` — in-place update + clearHistory on path reset

### D.2 `-ffast-math` Breaks NaN Guards (CMakeLists.txt)

**Problem:** Performance builds used `-ffast-math` which implies `-ffinite-math-only`,
causing the compiler to optimize away `std::isnan()` checks. GP trees routinely produce
NaN via `DIV(x,0)`, `SQRT(-1)`, etc., and NaN guards in `getInterpolatedTargetPosition()`
and elsewhere were silently eliminated.

**Fix (v1):** Append `-fno-finite-math-only` after `-ffast-math` — preserves NaN checks
but remaining sub-flags (`-funsafe-math-optimizations`, `-fassociative-math`,
`-freciprocal-math`) distort the fitness landscape on aarch64 via NEON FMA rounding.

**Fix (v2 — current):** Remove `-ffast-math` entirely from both CMakeLists.txt files.
Perf flags are now: `-O3 -g -march=native -mtune=native -funroll-loops -flto`.
On aarch64 (NVIDIA Grace / NEON), `-O3` alone provides sufficient speedup. The remaining
sub-flags interact with NEON FMA instructions to produce different rounding than x86 AVX,
which shifts the fitness landscape enough to break evolution.

**Throughput regression:** Removing `-ffast-math` drops throughput from ~5000 sims/sec to
~200 sims/sec on short paths (longSequential). Investigation shows this is not a codegen
issue — CPU utilization drops to ~5% per worker, indicating the bottleneck is dispatch
latency not computation speed. Longer paths (progressiveDistance) peg CPUs normally.
Root cause TBD — may be TCP round-trip overhead dominating short eval times on the
aarch64/Linux 6.17 kernel, or `-flto` pessimizing the dispatch path.

**Files changed:**
- `autoc/CMakeLists.txt`
- `crrcsim/CMakeLists.txt`

### D.3 Path Interpolator Float Precision (gp_evaluator_portable.cc, pathgen.h)

**Problem:** The path interpolation system has several 32-bit float precision issues that
may cause subtle non-determinism, especially on longer paths or under different optimizer
flags:

1. **`unsigned long` → `gp_scalar` (float) cast loses precision at high sim times.**
   At `gp_evaluator_portable.cc:354`, `aircraftState.getSimTimeMsec()` is `unsigned long`
   but cast to `gp_scalar` (float, 7 digits). After ~100s of sim time (100,000ms), float
   cannot represent millisecond differences. The binary search at line 328 uses exact `<=`
   comparison — a 1-ULP difference in the cast selects a different waypoint pair.

2. **Accumulated float odometer in path generation (`pathgen.h:67-94`).**
   Each waypoint's `simTimeMsec` is derived from `odometer / velocity * 1000`, where
   `odometer` is accumulated via `odometer += newDistance` in a loop of ~50+ iterations.
   Each iteration adds transcendental error from `(point - lastPoint).norm()`. The
   accumulated error is deterministic per-seed but fragile to optimizer reordering.

3. **Binary search has no epsilon tolerance (`gp_evaluator_portable.cc:322-333`).**
   Exact float comparison `path[mid].simTimeMsec <= goalTimeMsec` means a 1-ULP
   difference in either side selects a different bracket. Combined with (1) and (2),
   this can produce different interpolated positions on re-evaluation.

**Constraint:** `double` is too slow on target embedded hardware (Xiao-GP). Must stay
32-bit float (`gp_scalar`).

**Design considerations for fix:**
- The rabbit moves in one direction only along the path — no backtracking
- Search area is small: current position ± MAX_OFFSET_STEPS (±10 steps = ±1 second)
- Sample rate is 10Hz (100ms steps), so waypoint spacing is ~1.6m at 16 m/s
- Path is generated once per scenario and reused — accumulated error is consistent
  within a single evaluation, but may differ across optimizer flags or platforms
- A forward-scanning approach from last known index (instead of binary search from
  scratch each time) would avoid the float comparison sensitivity entirely
- Alternatively, use integer millisecond timestamps in the path and avoid float
  time comparison altogether

**Status:** Analysis only — no code changes yet. Revisit as part of this feature.

### D.4 Future TODO: Eval Loop Integration Tests

The current test suite validates individual sensor nodes but not the eval loop
orchestration (construct state → record history → evaluate → advance). This class of
bug (history wipe) would be caught by an integration test that:

1. Injects a known GP tree (e.g., `GETDPHI_PREV(1)`) into the eval loop
2. Runs 3+ ticks with known aircraft/path state
3. Asserts that temporal nodes return non-zero values after the first tick

This should be addressed as part of the unification effort (Phase 1/4), where the
eval loop can be decomposed into testable components. See also the minisim eval loop
which has the same pattern and could serve as a simpler test harness.

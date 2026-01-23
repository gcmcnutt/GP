# Pareto Fitness for Autoc

> **STATUS: PROPOSED (Jan 2026)**
>
> Multi-objective optimization to replace fragile weighted-sum fitness.
> Eliminates need for hand-tuned weights during evolution.

---

## The Problem: Knife-Edge Tuning

Current weighted-sum fitness is fragile:
- Small weight changes (e.g., throttle 1.0 -> 1.5) can destabilize evolution
- Weights tuned for one craft/path/wind combo don't transfer
- Evolution finds local optima that exploit weight imbalances (e.g., spiraling at full throttle)

**Observed**: With throttle weight 1.5, paths became very precise but craft couldn't complete them (crash/OOB). The selection pressure shift crushed previously viable strategies.

---

## Pareto Approach

### Core Idea

Instead of summing metrics with weights, keep them separate and use **Pareto dominance** for selection:

- Individual A **dominates** B if A is <= B on ALL objectives AND < on at least one
- No weights during evolution - selection purely based on "don't be dominated"
- Weights only applied at extraction time (picking bytecode for production)

### Benefits

1. **No weights during evolution** - selection pressure is purely "don't be dominated"
2. **Diversity preserved** - efficient-but-loose trackers coexist with tight-but-hungry trackers
3. **Robust to metric scaling** - dominance is ordinal, not magnitude-dependent
4. **Policy weights only at extraction** - easy to change what you ship without re-evolving

---

## Objective Hierarchy

### Priority Ordering (not weights, just understanding)

1. **Survivability** - Don't crash, don't get ahead of target
2. **Tracking** - Stay on path, move in correct direction
3. **Efficiency** - Minimize energy consumption

**Key insight**: Survivability isn't a separate Pareto objective - it's embedded in tracking/alignment. A crash truncates the run, degrading tracking and alignment averages naturally.

---

## Data Structures

### ParetoObjectives (3 objectives)

```cpp
// In autoc.h
struct ParetoObjectives {
    gp_fitness tracking;      // O1: waypoint_distance + movement_direction
    gp_fitness alignment;     // O2: angle to target (don't outrun rabbit)
    gp_fitness efficiency;    // O3: throttle energy

    int pareto_rank = 0;             // 0 = non-dominated front
    gp_fitness crowding_distance = 0.0;
};
```

### Why 3 Objectives Instead of 4?

| Original 4 Metrics | Issue | Resolution |
|-------------------|-------|------------|
| waypoint_distance | Correlated with cross_track | Combine into **tracking** |
| cross_track | Correlated with waypoint | Combine into **tracking** |
| movement_direction | Key for path following | Include in **tracking** |
| throttle_energy | Independent efficiency goal | Keep as **efficiency** |
| (new) angle_to_target | Don't get ahead of rabbit | New **alignment** objective |

Fewer objectives = smaller Pareto fronts = faster convergence.

---

## Core Algorithms

### Dominance Check

```cpp
bool dominates(const ParetoObjectives& a, const ParetoObjectives& b) {
    bool all_leq = (a.tracking <= b.tracking &&
                    a.alignment <= b.alignment &&
                    a.efficiency <= b.efficiency);
    bool any_lt = (a.tracking < b.tracking ||
                   a.alignment < b.alignment ||
                   a.efficiency < b.efficiency);
    return all_leq && any_lt;
}
```

### Non-Dominated Sorting (NSGA-II style)

```
Rank 0 (front 0): All individuals not dominated by anyone
Rank 1 (front 1): All individuals dominated only by rank 0
Rank 2 (front 2): All individuals dominated only by rank 0-1
... etc
```

### Selection

Fill new population by rank. Rank 0 first, then rank 1, etc.
Within a rank that doesn't fully fit, use **crowding distance** (prefer individuals in sparse regions of objective space).

---

## Evaluation Changes

Instead of summing to scalar fitness, store 3 objective values:

```cpp
// In evalTask(), at the end:

// O1: Tracking = waypoint distance + movement direction (both contribute to staying on path)
paretoObj.tracking = normalized_waypoint_distance + normalized_movement_direction;

// O2: Alignment = angle to target (don't get ahead of rabbit)
// angle_to_target: 0° = pointing at target (good), 180° = pointing away (bad)
paretoObj.alignment = normalized_angle_to_target;

// O3: Efficiency = throttle energy
paretoObj.efficiency = normalized_throttle_energy;

// For library compatibility, still compute stdFitness as simple sum
// (used for statistics display, not selection)
stdFitness = paretoObj.tracking + paretoObj.alignment + paretoObj.efficiency;
```

### Crash Handling

Crashes affect objectives naturally:

```cpp
// Crash truncates run - tracking/alignment averages degrade automatically
// since the run is shorter and final position is likely off-path

if (crashReason != CrashReason::None) {
    // Optional: Add penalty to alignment objective specifically
    // (OOB = failure to track the rabbit)
    gp_fitness fraction_remaining = 1.0 - (distanceTraveled / total_path_distance);
    paretoObj.alignment += OOB_ALIGNMENT_PENALTY * fraction_remaining;
}
```

**Key insight**: We don't need a separate crash penalty objective. Crash truncates the run, which:
1. Degrades tracking average (less path covered)
2. Degrades alignment average (likely off-target when crashed)
3. Optional small penalty on alignment for OOB specifically

---

## Selection Override

In `endOfEvaluation()`:

```cpp
void MyPopulation::endOfEvaluation() {
    // 1. Collect all ParetoObjectives from population
    std::vector<std::pair<int, ParetoObjectives*>> pop;
    for (int i = 0; i < containerSize(); i++) {
        pop.push_back({i, &NthGP(i)->paretoObj});
    }

    // 2. Non-dominated sorting
    assignParetoRanks(pop);

    // 3. Compute crowding distance within each rank
    computeCrowdingDistances(pop);

    // 4. Sort by (rank ASC, crowding DESC)
    std::sort(pop.begin(), pop.end(), [](auto& a, auto& b) {
        if (a.second->pareto_rank != b.second->pareto_rank)
            return a.second->pareto_rank < b.second->pareto_rank;
        return a.second->crowding_distance > b.second->crowding_distance;
    });

    // 5. Override stdFitness with rank-based value for library selection
    // Lower rank = better fitness, crowding as tiebreaker
    for (int i = 0; i < pop.size(); i++) {
        NthGP(pop[i].first)->stdFitness =
            pop[i].second->pareto_rank * 10000.0 -
            pop[i].second->crowding_distance;
    }
}
```

---

## Best Extraction (for S3/Production)

Pick from rank 0 (Pareto front) using policy weights:

```cpp
MyGP* selectBestForProduction(std::vector<MyGP*>& front) {
    // Policy: survivability first (tracking + alignment), efficiency secondary
    const float w_tracking = 2.0f;    // Stay on path (safety critical)
    const float w_alignment = 2.0f;   // Don't outrun rabbit (safety critical)
    const float w_efficiency = 0.5f;  // Energy (nice to have)

    MyGP* best = nullptr;
    float best_score = FLT_MAX;
    for (auto* gp : front) {
        float score = w_tracking * gp->paretoObj.tracking +
                      w_alignment * gp->paretoObj.alignment +
                      w_efficiency * gp->paretoObj.efficiency;
        if (score < best_score) {
            best_score = score;
            best = gp;
        }
    }
    return best;
}
```

**Key insight**: Weights only appear when you need a single answer for deployment, never during evolutionary selection pressure.

---

## Tracking "Best" Across Generations

With Pareto, there's no single "best" - you have a Pareto front. Options:

| Approach | Description | Recommendation |
|----------|-------------|----------------|
| **Hypervolume** | Single scalar capturing front quality | Complex, reference-point dependent |
| **Weighted sum for tracking** | Apply policy weights only for S3 storage | Simple, recommended |
| **Achievement scalarization** | Distance to target goals | More intuitive than weights |
| **Lexicographic** | Strict priority ordering | Simple but inflexible |

**Recommended**: Use weighted sum only at extraction time. Evolution is adaptive, but we still get a single candidate for production.

---

## Implementation Phases

### Phase 1: Data Structure
- Add `ParetoObjectives` struct to `autoc.h` (3 objectives: tracking, alignment, efficiency)
- Add angle-to-target calculation in `evalTask()` for alignment metric
- Populate objectives in `evalTask()` alongside existing fitness
- Keep current selection unchanged
- Verify values are reasonable in logs

### Phase 2: Pareto Ranking
- Add dominance check function
- Add non-dominated sorting in `endOfEvaluation()`
- Override `stdFitness` with rank-based value
- Verify selection behavior changes

### Phase 3: Crash Handling
- Add optional OOB penalty to alignment objective
- Verify crash truncation naturally degrades tracking/alignment
- Tune `OOB_ALIGNMENT_PENALTY` if needed

### Phase 4: Best Extraction
- Update S3 storage to use policy-weighted selection from front
- Add Pareto front size to logging
- Optionally log hypervolume for front quality tracking

---

## Expected Behavior Changes

| Current (weighted sum) | Pareto |
|------------------------|--------|
| Single "best" per generation | Pareto front of non-dominated solutions |
| Weight changes destabilize evolution | Weight-free selection is robust |
| Specialists get crushed | Diverse strategies coexist |
| Tuning required per craft/path | No tuning - policy weights only at extraction |

---

## Logging Additions

```
PARETO_FRONT: gen=42 size=127 hypervolume=0.847
PARETO_BEST: gen=42 tracking=18.5 alignment=12.3 efficiency=45.6 policyScore=89.4
```

---

## TODO / Future Investigations

### 4th Objective (Brevity/Duration) - FAILED EXPERIMENT (Jan 2026)

**Attempted:** Add flight duration as 4th Pareto objective (brevity = 1/duration, minimize).

**Problems encountered:**
1. With 4 objectives, Pareto front grew to ~900/1000 individuals - no selection pressure
2. Policy score weighting couldn't balance duration vs tracking - either duration dominated (picked poor trackers that survived long) or tracking dominated (picked short-lived good trackers)
3. Priority-based weights (1000/100/10/1) didn't help because objective value ranges differ wildly (tracking: 10-160, brevity: 0.08-0.6)
4. Fundamental tension: bakeoff/elite selection uses policy score which must somehow weight all objectives, defeating the "no weights during evolution" goal

**Conclusion:** Duration as a Pareto objective doesn't work well. Consider:
- Constraint-based approach (must survive N seconds to be feasible)
- Multiplicative fitness (tracking * duration) rather than additive
- Lexicographic pre-filtering before Pareto ranking
- Revisit when progressive difficulty is implemented

---

### Interleaved Step-Wise Fitness Evaluation

**Current approach:** Run individual on scenario A to completion, measure fitness. Run on scenario B to completion, measure fitness. Average the fitnesses.

**Alternative approach:** Interleave steps across scenarios. Step 1 of A, step 1 of B, step 1 of C... then step 2 of A, step 2 of B... Sum/average the per-step fitness contributions across all scenarios simultaneously.

**Potential benefits:**
- Smooths out knife-edge crashes - one early crash doesn't dominate the signal
- Reduces scenario-specific overfitting
- More stable fitness gradients
- Better for progressive difficulty - harder scenarios contribute from step 1

**Potential drawbacks:**
- Loses temporal coherence for diagnostics
- Masks scenario-specific failure modes
- Implementation complexity (lockstep, different-length scenarios)
- May reward mediocrity over specialists

**Context:** At 10Hz over 15 seconds = 150 steps × 6 scenarios = 900 fitness samples per individual. Interleaved approach has much higher resolution for distinguishing "almost survived" from "crashed immediately" across the ensemble. May pair well with progressive difficulty.

---

### Trimmed Mean Fitness Aggregation

**Current approach:** Average fitness across all N scenarios (path × wind × ...).

**Problem:** One catastrophic failure or one lucky easy scenario can skew the average. Encourages "safe mediocrity" rather than robust competence.

**Alternative: Trimmed mean (like figure skating judging)**
- Sort scenario fitnesses
- Drop worst X% and best Y%
- Average the remaining middle scenarios

**Trim policy options:**
```cpp
struct TrimPolicy {
    float trim_bottom_pct = 0.10f;  // Drop worst 10%
    float trim_top_pct = 0.05f;     // Drop best 5%
};
```

With 36 scenarios (6 paths × 6 winds): drop ~4 worst, ~2 best, average remaining 30.

**Future dimensions:** path × wind × craft × entry × ...
- Each dimension can vary independently (6 paths, 5 winds, 7 crafts, 4 entries, etc.)
- Could do stratified trimming within dimensions
- Could weight dimensions differently

**Advanced options:**
1. **Stratified trimming** - trim within each dimension, then combine
2. **Weighted dimensions** - craft variation weighted more than entry point
3. **Progressive trimming** - generous early (20%), tight late (5%)
4. **Failure-aware** - never trim crashes, only trim "meh" scenarios

---

### Fitness Normalization / Difficulty Factors

**Current:** Normalize by path length (time-weighted average). Longer paths contribute proportionally more.

**Not addressed:** Path difficulty varies. A figure-8 through turbulence is harder than a straight line in calm air.

**Options to consider:**
1. **Difficulty factors** - like skating, multiply raw score by path complexity
   - Path already computes `total_angular_change` - could use this
   - Harder paths get higher difficulty multiplier
2. **Percentile ranking within paths** - compare GP's performance to population baseline on each path
3. **Historical normalization** - track average fitness per scenario over time, normalize against it

**Tradeoff:** Difficulty factors add tunables. May be better to let trimmed mean handle outliers naturally.

---

## References

- NSGA-II: Deb et al., "A Fast and Elitist Multiobjective Genetic Algorithm"
- Crowding distance maintains diversity on Pareto front
- Tournament selection using Pareto rank + crowding as tiebreaker

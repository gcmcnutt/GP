# Making Demetic Mode Work Correctly

## The Problem

In demetic mode, the GP library's `calculateStatistics()` overwrites `bestOfPopulation`
based on single-scenario `stdFitness` values, which are not comparable across demes.

### Root Cause

1. Each deme evaluates on a **different scenario** (1 path × 6 winds = 6 sims)
2. The bakeoff evaluates deme winners on **all scenarios** (6 paths × 6 winds = 36 sims)
3. After bakeoff, we correctly set `bestOfPopulation = gGlobalEliteIndex`
4. But `GPPopulation::generate()` calls `calculateStatistics()` which scans all individuals
   and picks whoever has the lowest `stdFitness` - comparing apples to oranges

### Evidence from autoc_130.log

```
gen=6: 158.940704 (best fitness from bakeoff)
gen=7: 166.824600 (WORSE - elite replaced incorrectly)
gen=8: 169.020538 (even worse)
gen=9: 170.631287 (continuing to degrade)
```

The `worse_replacements` counter correctly detected these violations.

---

## Proposed Fix: Fitness Tuple

### Option A: Library Modification (Clean)

Add secondary fitness field to GP class:

```cpp
// In GP class (~/GP/include/gp.h)
class GP {
public:
  double stdFitness;           // Existing: deme-local fitness
  double aggregatedFitness;    // NEW: bakeoff fitness (inf until evaluated)
  bool hasAggregatedFitness;   // NEW: true after bakeoff evaluation
};
```

Modify `calculateStatistics()` in `~/GP/src/generate.cc`:

```cpp
void GPPopulation::calculateStatistics() {
  for (int n = 0; n < containerSize(); n++) {
    GP* current = NthGP(n);

    // For best: prefer aggregated fitness if available
    double currentEffective = current->hasAggregatedFitness
      ? current->aggregatedFitness
      : current->stdFitness;

    double bestEffective = best->hasAggregatedFitness
      ? best->aggregatedFitness
      : best->stdFitness;

    if (currentEffective < bestEffective ||
        (currentEffective == bestEffective &&
         current->length() < best->length())) {
      bestOfPopulation = n;
      best = current;
    }
  }
}
```

### Option B: Override in MyPopulation (Quick)

```cpp
class MyPopulation : public GPPopulation {
  void calculateStatistics() override {
    GPPopulation::calculateStatistics();  // Let library run

    // Restore global elite if it's better than library's pick
    if (gGlobalEliteIndex >= 0 && std::isfinite(gGlobalEliteFitness)) {
      MyGP* libraryBest = NthMyGP(bestOfPopulation);
      if (gGlobalEliteFitness < libraryBest->getFitness()) {
        bestOfPopulation = gGlobalEliteIndex;
      }
    }
  }
};
```

---

## Alternative: Rethink Bakeoff Strategy

The current bakeoff approach has fundamental issues:

1. **Fitness semantics are confused**: `stdFitness` means different things in different contexts
2. **Selection pressure is inconsistent**: Deme selection uses 6-sim fitness, but elitism uses 36-sim
3. **Migration is awkward**: Propagating a generalist to specialist demes may not be optimal

### Possible Redesigns

#### A. Full Evaluation Mode (Current non-deme)
- Every individual evaluated on all 36 scenarios
- Simple, correct, but 6x slower per generation
- Current performance: ~1100 sims/sec, ~5 min/gen

#### B. Two-Phase Fitness
- Phase 1: Deme evaluation (6 sims) for selection within deme
- Phase 2: All deme winners get full evaluation (36 sims) for global ranking
- Store both fitness values explicitly
- Selection uses deme fitness, elitism uses global fitness

#### C. Rolling Evaluation
- Each generation, evaluate on a rotating subset of scenarios
- Accumulate fitness over multiple generations
- More exploration, less exploitation per generation

#### D. Island Model (True Demes) — Recommended

The purest approach: let demetic migration do its natural work without forcing a global champion.

**Core idea**: Each deme evolves independently with its own local elite. Migration spreads
good genetic material naturally. No per-generation bakeoff needed.

**Key question**: What individual should we store for production (xiao flight code)?

**Options for Production Output**:

1. **Periodic Full Evaluation** (Recommended)
   - At checkpoint intervals (e.g., every 10 gens), evaluate each deme's best on all scenarios
   - Store whoever has best aggregated fitness as the production candidate
   - Decouples evolution speed from production selection
   ```cpp
   if (gen % checkpointInterval == 0) {
     for (each deme) {
       fullEval(demeBest, allScenarios);
     }
     storeToS3(bestOverallGeneralist);
   }
   ```

2. **Store All Deme Bests**
   - Store each deme's champion at every generation
   - Post-process offline: load all, run full evaluation tournament
   - Most flexible but requires offline pipeline

3. **Ensemble/Voting** (Future)
   - Don't pick one — use all deme specialists with a meta-controller
   - Each specialist votes on control output; use median/weighted average
   - Exploits specialization rather than fighting it

**Why Island Model works better**:
- No fitness semantic confusion (all comparisons are within-deme, apples-to-apples)
- Migration spreads good genes naturally without forced "generalist wins" logic
- Evolution runs at full demetic speed (~560 sims/sec)
- Production selection is a separate concern (periodic full eval)

---

## Recommendation

1. **Short term**: Use non-deme mode for production runs (it works correctly, ~240 sims/sec)
2. **Medium term**: Implement Island Model (Option D) with periodic full evaluation checkpoints
   - Remove per-gen bakeoff overhead
   - Let demes evolve independently with natural migration
   - Every N generations, run full eval tournament to identify production candidate
3. **Skip**: Option A (library fitness tuple) — adds complexity without solving the fundamental issue

The current demetic implementation tries to get benefits of both approaches but ends up
with the downsides of both. The Island Model is the cleanest path forward:
- **Evolution**: Pure demetic with migration (fast, ~560 sims/sec)
- **Production selection**: Separate concern (periodic full eval checkpoint)

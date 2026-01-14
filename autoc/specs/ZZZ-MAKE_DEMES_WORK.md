# Making Demetic Mode Work Correctly

> **STATUS: DEPRIORITIZED (Jan 2026)**
>
> After analysis, the demetic approach with specialized demes (one path per deme) may not
> be the best path forward. The core issue is that specialists don't transfer well between
> demes, so migration adds noise rather than spreading useful genes.
>
> **Alternative approaches to explore instead:**
>
> 1. **Progressive Fitness Checks** - Start with simple scenarios, add complexity over time
>    as the population improves. Early generations use fewer/simpler paths, later generations
>    add harder scenarios. This naturally focuses compute on promising individuals.
>
> 2. **Sparse Sampling** - Rather than full Cartesian product (paths × winds × craft × ...),
>    randomly sample combinations each generation. Over time, good generalists emerge
>    because they're tested on varied (but not exhaustive) scenarios.
>
> 3. **Adaptive Scenario Selection** - Focus evaluation on scenarios where the population
>    is weakest. If most individuals handle path-1 well but fail on path-3, spend more
>    compute on path-3 variants.
>
> These approaches keep the single-population model (no deme complexity) while improving
> efficiency. See future specs for implementation details.

---

## Current State (as of Jan 2026)

### What Gets Emitted Currently

**data.stc** (statistics file):
```
#Gen|       Fitness     |      Length       |   Depth
#   |  Best|Avg.|Worst  |  Best|Avg.|Worst  |  Best|Avg.|Worst
0 158.940704 234.123456 512.345678    42 156 892    8 12 24
1 155.234567 228.567890 498.765432    38 148 856    7 11 23
#GenSimStats gen=1 sims=12800 total=12800 durationSec=45.23 rate=283.02
...
```

**data.dat** (full GP dump):
```
Best of generation 0 (Fitness = 158.940704, Structural Complexity = 42)
(+ (* x y) (sin z) ...)

Best of generation 1 (Fitness = 155.234567, Structural Complexity = 38)
...
```

**Logger output** (ELITE_STORE lines):
```
ELITE_STORE: gen=1 gpHash=0xabc123->0xdef456 fitness=158.940704->155.234567 len=42->38 depth=8->7 REPLACED
ELITE_STATUS: gen=1 reevals=0 divergences=0
```

**S3 storage**: `{timestamp}/gen{N}.dmp` containing serialized `EvalResults` with:
- GP bytecode
- Fitness
- Path list and scenario list
- Aircraft state trajectories

### The Problem (Unchanged)

In demetic mode, the GP library's `calculateStatistics()` overwrites `bestOfPopulation`
based on single-scenario `stdFitness` values, which are not comparable across demes.

Each deme evaluates on a **different scenario** (1 path × N winds = N sims).
The library then compares these incomparable fitness values when picking `bestOfPopulation`.

---

## Proposed Solution: True Island Mode with Best-of-Best Tracking

### Core Design

**Principle**: Let demes evolve independently with their own local elites. No per-generation
bakeoff that forces elite propagation. Migration spreads genetic material naturally.
The "best-of-best" is tracked separately for production output but doesn't feed back into evolution.

### Two Fitness Concepts

| Concept | Description | Used For |
|---------|-------------|----------|
| **Deme Elite** | Best individual within each deme by local fitness | Selection, crossover, mutation within deme |
| **Best-of-Best** | Best individual across all demes by aggregated fitness | Production output (S3 storage), logging |

### Key Difference from Current Approach

**Current (broken)**: Bakeoff winner's aggregated fitness overwrites `stdFitness`, then
propagated to other demes. This corrupts local selection pressure.

**Proposed**: Bakeoff determines best-of-best for storage/logging ONLY. Each deme keeps
its own elite with local fitness intact. No cross-deme elite propagation.

### Outputs in Island Mode

#### data.stc Changes

Add a separate line for best-of-best after the per-generation statistics:

```
#Gen|       Fitness     |      Length       |   Depth
#   |  Best|Avg.|Worst  |  Best|Avg.|Worst  |  Best|Avg.|Worst
0 158.940704 234.123456 512.345678    42 156 892    8 12 24
#BestOfBest gen=0 fitness=162.345678 len=42 depth=8 deme=3
#GenSimStats gen=0 sims=12800 total=12800 durationSec=45.23 rate=283.02
1 152.234567 228.567890 498.765432    38 148 856    7 11 23
#BestOfBest gen=1 fitness=158.234567 len=38 depth=7 deme=1
#GenSimStats gen=1 sims=12800 total=12800 durationSec=44.87 rate=285.32
...
```

**Note**: The "Best" fitness in the main line is the **best deme-local fitness** (apples-to-oranges,
not globally meaningful). The `#BestOfBest` line shows the **aggregated fitness** (apples-to-apples)
which may go up OR down between generations.

#### data.dat Changes

Add best-of-best GP dump after each generation:

```
Best of generation 1 (Fitness = 152.234567, Structural Complexity = 38)
(deme-local best - may not be globally best)
(+ (* x y) (sin z) ...)

Best-of-Best generation 1 (Aggregated Fitness = 158.234567, Structural Complexity = 38, Deme = 1)
(evaluated on all scenarios)
(+ (* x y) (sin z) ...)
```

#### Logger Output Changes

Keep existing ELITE_STORE for deme-local tracking, add BESTOFBEST and HISTORICALBEST lines:

```
ELITE_STORE: gen=1 gpHash=0xabc123->0xdef456 fitness=152.234567->149.345678 len=38->36 depth=7->6 REPLACED
BESTOFBEST: gen=1 fitness=158.234567->161.456789 len=38->40 depth=7->8 deme=1->2 WORSENED
HISTORICALBEST: gen=1 fitness=158.234567 len=38 depth=7 deme=1 fromGen=1 (unchanged)
ELITE_STATUS: gen=1 reevals=0 divergences=0
```

BESTOFBEST outcomes (current gen winner, can fluctuate):
- `IMPROVED` - aggregated fitness decreased (better than last gen's winner)
- `WORSENED` - aggregated fitness increased (worse than last gen's winner, expected)
- `SAME` - same GP won again

HISTORICALBEST outcomes (monotonically improving):
- `IMPROVED` - new all-time best found this generation
- `(unchanged)` - historical best remains from earlier generation

### S3 Storage

**No change** to what gets stored - still store the best-of-best (winner of aggregated bakeoff).
This is the production candidate for flight code.

---

## Implementation Changes

### 1. Remove Elite Propagation in Demetic Mode

In `MyPopulation::endOfEvaluation()`, remove the demetic propagation block (lines 604-645).
The bakeoff winner should NOT be propagated to other demes.

```cpp
// REMOVE THIS BLOCK:
// Demetic propagation (only in demetic mode)
if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0 &&
    std::isfinite(bestAggregatedFitness)) {
  // ... propagation code ...
}
```

### 2. Don't Overwrite Bakeoff Winner's Fitness

Currently (line 601):
```cpp
best->setFitness(bestAggregatedFitness);
```

This corrupts the deme-local fitness. Instead, track aggregated fitness separately:

```cpp
// NEW: Track best-of-best separately, don't modify stdFitness
static double gBestOfBestFitness = std::numeric_limits<double>::infinity();
static int gBestOfBestIndex = -1;
static int gBestOfBestDeme = -1;

// In bakeoff winner determination:
if (averageFitness < bestAggregatedFitness) {
  bestAggregatedFitness = averageFitness;
  bestOfEvalResults = aggregatedEvalResults;
  best = candidate;
  gBestOfBestIndex = candidateIndex;
  gBestOfBestDeme = candidateIndex / gpCfg.DemeSize;
}

// DO NOT: best->setFitness(bestAggregatedFitness);
// DO NOT: bestOfPopulation = bestIndex;
```

### 3. Let Library Pick Deme-Local Best

Remove the line that sets `bestOfPopulation` to bakeoff winner (line 602):
```cpp
// REMOVE: bestOfPopulation = bestIndex;
```

Let `calculateStatistics()` pick whatever it picks. This value is meaningless
in demetic mode (comparing apples to oranges) but harmless.

### 4. Add Best-of-Best and Historical-Best Logging

Add global tracking variables:

```cpp
// Best-of-Best tracking (current generation's aggregated winner - can fluctuate)
static double gBestOfBestFitness = std::numeric_limits<double>::infinity();
static int gBestOfBestIndex = -1;
static int gBestOfBestDeme = -1;
static int gBestOfBestLength = 0;
static int gBestOfBestDepth = 0;

// Historical Best tracking (all-time best - monotonically improves)
static double gHistoricalBestFitness = std::numeric_limits<double>::infinity();
static int gHistoricalBestGen = 0;
static int gHistoricalBestDeme = -1;
static int gHistoricalBestLength = 0;
static int gHistoricalBestDepth = 0;
static std::vector<char> gHistoricalBestGP;  // Serialized GP for S3 storage

// Previous generation's best-of-best (for comparison logging)
static double gLastBestOfBestFitness = std::numeric_limits<double>::infinity();
static int gLastBestOfBestDeme = -1;
static int gLastBestOfBestLength = 0;
static int gLastBestOfBestDepth = 0;
```

After `createGenerationReport()`, add best-of-best and historical-best output:

```cpp
// After existing createGenerationReport() call:
pop->createGenerationReport(0, gen, fout, bout, *logStream);

// NEW: Best-of-best and historical-best logging for demetic mode
if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0 && gBestOfBestIndex >= 0) {
  MyGP* bob = pop->NthMyGP(gBestOfBestIndex);

  // Determine BestOfBest outcome (can go up or down)
  const char* bobOutcome;
  if (gen == 1) {
    bobOutcome = "INITIAL";
  } else if (gBestOfBestFitness < gLastBestOfBestFitness) {
    bobOutcome = "IMPROVED";
  } else if (gBestOfBestFitness > gLastBestOfBestFitness) {
    bobOutcome = "WORSENED";
  } else {
    bobOutcome = "SAME";
  }

  // data.stc line for BestOfBest
  bout << "#BestOfBest gen=" << gen
       << " fitness=" << gBestOfBestFitness
       << " len=" << bob->length()
       << " depth=" << bob->depth()
       << " deme=" << gBestOfBestDeme
       << " " << bobOutcome
       << endl;

  // Logger line for BestOfBest
  *logger.info() << "BESTOFBEST: gen=" << gen
                 << " fitness=" << gLastBestOfBestFitness << "->" << gBestOfBestFitness
                 << " len=" << gLastBestOfBestLength << "->" << bob->length()
                 << " depth=" << gLastBestOfBestDepth << "->" << bob->depth()
                 << " deme=" << gLastBestOfBestDeme << "->" << gBestOfBestDeme
                 << " " << bobOutcome << endl;

  // Check if this is a new historical best
  bool newHistoricalBest = (gBestOfBestFitness < gHistoricalBestFitness);
  if (newHistoricalBest) {
    gHistoricalBestFitness = gBestOfBestFitness;
    gHistoricalBestGen = gen;
    gHistoricalBestDeme = gBestOfBestDeme;
    gHistoricalBestLength = bob->length();
    gHistoricalBestDepth = bob->depth();

    // Serialize GP for S3 storage
    gHistoricalBestGP.clear();
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>>
      outStream(gHistoricalBestGP);
    bob->save(outStream);
    outStream.flush();
  }

  // data.stc line for HistoricalBest
  bout << "#HistoricalBest gen=" << gen
       << " fitness=" << gHistoricalBestFitness
       << " len=" << gHistoricalBestLength
       << " depth=" << gHistoricalBestDepth
       << " deme=" << gHistoricalBestDeme
       << " fromGen=" << gHistoricalBestGen
       << (newHistoricalBest ? " IMPROVED" : "")
       << endl;

  // Logger line for HistoricalBest
  *logger.info() << "HISTORICALBEST: gen=" << gen
                 << " fitness=" << gHistoricalBestFitness
                 << " len=" << gHistoricalBestLength
                 << " depth=" << gHistoricalBestDepth
                 << " deme=" << gHistoricalBestDeme
                 << " fromGen=" << gHistoricalBestGen
                 << (newHistoricalBest ? " IMPROVED" : " (unchanged)") << endl;

  // data.dat entry for BestOfBest
  fout << "Best-of-Best generation " << gen
       << " (Aggregated Fitness = " << gBestOfBestFitness
       << ", Structural Complexity = " << bob->length()
       << ", Deme = " << gBestOfBestDeme
       << ")" << endl
       << *bob << endl;

  // Update last-gen tracking for next comparison
  gLastBestOfBestFitness = gBestOfBestFitness;
  gLastBestOfBestLength = bob->length();
  gLastBestOfBestDepth = bob->depth();
  gLastBestOfBestDeme = gBestOfBestDeme;
}
```

### 5. Adjust ELITE_STORE for Demetic Mode

In demetic mode, ELITE_STORE should track the **deme-local best** (what library picks),
not the bakeoff winner. This may show fitness going up/down which is expected
when comparing across demes.

Add a note to the log:
```cpp
if (gpCfg.DemeticGrouping) {
  *logger.info() << "ELITE_STORE: (deme-local, cross-deme comparison not meaningful)" << endl;
}
```

### 6. S3 Storage: Both BestOfBest and HistoricalBest

Store two files per generation:

```cpp
// Always store current gen's BestOfBest (for analysis/debugging)
// Key: {timestamp}/gen{N}.dmp
storeToS3(computedKeyName, bestOfEvalResults);

// Store HistoricalBest only when it improves (for production)
// Key: {timestamp}/best.dmp
if (newHistoricalBest) {
  EvalResults historicalResults = bestOfEvalResults;  // Copy current eval results
  historicalResults.gp = gHistoricalBestGP;           // But use historical best GP
  storeToS3(startTime + "/best.dmp", historicalResults);

  *logger.info() << "New historical best stored to S3: " << startTime << "/best.dmp" << endl;
}
```

This gives us:
- **gen{N}.dmp**: Full history of each generation's aggregated winner (can fluctuate)
- **best.dmp**: Always contains the best generalist found so far (monotonically improves)

---

## Expected Behavior

### Deme-Local (ELITE_STORE line)

- Fitness may go up or down between generations
- Comparing apples to oranges (different scenarios)
- Used internally by library, not meaningful for production selection

### Best-of-Best vs Historical Best

**BestOfBest** (current generation's aggregated winner):
- Aggregated fitness (all scenarios) - apples-to-apples comparison
- **May go up or down** as different demes produce winners
- This is normal and expected:
  - Deme 1's specialist may beat Deme 2's specialist on aggregated fitness in gen N
  - But Deme 2's specialist may beat Deme 1's in gen N+1

**HistoricalBest** (best we've ever seen):
- Tracks the **monotonically improving** best across all generations
- Only updates when current BestOfBest beats previous HistoricalBest
- This is what we actually care about for production

Example progression:
```
Gen 1: BestOfBest=162.3 (deme 2) HistoricalBest=162.3 fromGen=1 INITIAL
Gen 2: BestOfBest=158.1 (deme 4) HistoricalBest=158.1 fromGen=2 IMPROVED
Gen 3: BestOfBest=165.7 (deme 1) HistoricalBest=158.1 fromGen=2 (unchanged)
Gen 4: BestOfBest=155.2 (deme 4) HistoricalBest=155.2 fromGen=4 IMPROVED
```

### Migration: Does It Even Help Here?

The GP library implements **pairwise swap** between adjacent demes on a linear chain:

```cpp
// Library's demeticMigration():
for each adjacent pair (deme N, deme N+1):
  if (random < DemeticMigProbability):
    swap(bestOf(deme N), bestOf(deme N+1))
```

**Problem**: In our model, each deme evaluates on a **different scenario** (path geometry).
A specialist optimized for path-1 will likely perform poorly when dropped into path-2's deme.
The swap may just inject noise rather than spreading useful genetic building blocks.

**Recommendation**: Consider setting `DemeticMigProbability = 0` to disable library migration.
Let demes evolve as pure specialists. The aggregated bakeoff identifies which specialist
happens to generalize best, but doesn't corrupt the specialist populations.

If migration is desired later, consider a custom approach:
- Migrate only individuals that show good **aggregated** fitness (generalists)
- Or migrate genetic subtrees rather than whole individuals

### S3 Storage Strategy

**Store both**:
- `{timestamp}/gen{N}.dmp` - current generation's BestOfBest (for analysis)
- `{timestamp}/best.dmp` - overwritten only when HistoricalBest improves (for production)

This gives us:
1. Full history of how specialists evolved
2. Always-current best generalist for flight code

---

## Summary of Changes

| File | Change |
|------|--------|
| autoc.cc | Remove elite propagation block (lines 604-645) |
| autoc.cc | Don't overwrite `stdFitness` with aggregated fitness (line 601) |
| autoc.cc | Don't set `bestOfPopulation` to bakeoff winner (line 602) |
| autoc.cc | Add global tracking: `gBestOfBest*` and `gHistoricalBest*` variables |
| autoc.cc | Add `#BestOfBest` and `#HistoricalBest` lines to data.stc |
| autoc.cc | Add `BESTOFBEST:` and `HISTORICALBEST:` lines to logger output |
| autoc.cc | Add best-of-best entry to data.dat |
| autoc.cc | Store `best.dmp` to S3 when historical best improves |
| config | Consider setting `DemeticMigProbability=0` (disable library migration) |

---

## Verification Checklist

After implementation:

1. [ ] Deme-local fitness in ELITE_STORE may go up/down (expected)
2. [ ] Best-of-best fitness in BESTOFBEST may go up/down (expected)
3. [ ] Historical-best fitness in HISTORICALBEST only improves (monotonic)
4. [ ] S3 `gen{N}.dmp` contains current gen's best-of-best
5. [ ] S3 `best.dmp` contains historical best (updated only on improvement)
6. [ ] No elite propagation across demes
7. [ ] Library migration disabled or tested with re-evaluation
8. [ ] Non-demetic mode unchanged (regression test)

---

## Appendix: Library Migration vs Our Model

### What the GP Library Implements

The library's `demeticMigration()` implements a **linear chain** (not full torus) with
**pairwise swap** between adjacent demes:

```cpp
// From GP/src/generate.cc lines 331-364
void GPPopulation::demeticMigration() {
  // For each adjacent pair of demes (0-1, 1-2, 2-3, ...)
  // Note: Loop stops before last deme, so no wrap-around
  for (int demeStart=0; demeStart < containerSize() - DemeSize; demeStart += DemeSize) {
    if (GPRandomPercent(DemeticMigProbability)) {
      // Select best from deme N and deme N+1
      int r1 = selectBest(deme N);
      int r2 = selectBest(deme N+1);
      // SWAP them (bidirectional exchange)
      swap(population[r1], population[r2]);
    }
  }
}
```

### Why This May Not Help Our Model

| Library Assumption | Our Reality |
|--------------------|-------------|
| All demes use same fitness landscape | Each deme has different path geometry |
| "Good" in deme A = "Good" in deme B | Specialist for path-1 may fail on path-2 |
| Migration spreads useful building blocks | Migration may just inject noise |

### Recommendation

**Option 1: Disable migration** (`DemeticMigProbability=0`)
- Let demes evolve as pure specialists
- Bakeoff identifies best generalist without corrupting specialists
- Cleanest separation of concerns

**Option 2: Custom migration** (future enhancement)
- Only migrate individuals with good **aggregated** fitness
- Or migrate genetic subtrees rather than whole individuals
- Requires custom implementation outside library

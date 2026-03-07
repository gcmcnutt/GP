# Single PRNG Architecture

## Overview

Refactor autoc to use a single PRNG source (GPrand) for all randomness, eliminating multiple independent PRNGs that create correlation issues and make determinism harder to verify.

## Problem Statement

### Current State: Multiple PRNGs

| PRNG | Seeded From | Location | Used For |
|------|-------------|----------|----------|
| GPrand() | gpSeed | GP library | GP evolution (mutation, crossover, selection) |
| Inline LCG | windSeed | variation_generator.h | Entry variations, rabbit speed profiles |
| std::mt19937 | random_device | autoc.cc:628 | Demic migration (**NON-DETERMINISTIC BUG**) |
| CRRC_Random | windSeed | crrcsim | Wind, thermals, gust |

### Issues

1. **Multiple PRNGs from same seed**: Entry variations and crrcsim wind both start from `windSeed`, creating hidden correlations
2. **Non-determinism bug**: Demic migration uses `std::random_device{}()` - not reproducible
3. **Complex seed management**: WindSeedBase, WindSeedStride computed separately from GPSeed
4. **Verification difficulty**: Hard to prove determinism with multiple PRNG instances
5. **PRNG continuation bug**: `generateVariations()` takes seed by value, so rabbit speed doesn't continue from entry variations

## Proposed Architecture

### Single Source of Truth

```
GPInit(1, gpSeed)   // ONE seed controls EVERYTHING

At startup (before GP evolution):
  1. pathSeed = GPrand()                    // For random path methods
  2. For each scenario i = 0..N-1:
       windSeeds[i] = GPrand()              // Seed for crrcsim CRRC_Random
       entryVariations[i] = generate()      // Consumes ~10 GPrand values
       rabbitSpeedProfiles[i] = generate()  // Consumes ~50-200 GPrand values

GP evolution then continues consuming GPrand()
Pre-fetched values REUSED every generation
```

### Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│                         GPSeed                               │
│                           │                                  │
│                      GPInit(gpSeed)                          │
│                           │                                  │
│                       GPrand()                               │
│                           │                                  │
│    ┌──────────────────────┼──────────────────────┐          │
│    │                      │                      │          │
│    ▼                      ▼                      ▼          │
│ pathSeed            windSeeds[N]           GP Evolution     │
│    │               entryVars[N]            (continues       │
│    │               rabbitProfiles[N]        consuming)      │
│    │                      │                                  │
│    ▼                      ▼                                  │
│ Path geometry        Passed to crrcsim                      │
│ (std::mt19937)       (CRRC_Random seeded                    │
│                       from windSeeds[i])                    │
└─────────────────────────────────────────────────────────────┘
```

### Special Case: RandomPathSeedB

`RandomPathSeedB` is **retained** as a separate config parameter. It allows reproducing specific path geometries independent of GPSeed.

When `RandomPathSeedB != -1`:
- Use configured value for path generation (override GPrand-derived pathSeed)
- Useful for debugging specific path behaviors

When `RandomPathSeedB == -1`:
- Derive pathSeed from GPrand() at startup

## Implementation

### Phase 1: Pre-computation Infrastructure

#### New Data Structures

```cpp
// autoc.h or new prefetch.h
struct ScenarioVariations {
    unsigned int windSeed;           // For crrcsim CRRC_Random
    VariationOffsets entryOffsets;   // Heading, roll, pitch, speed
    std::vector<RabbitSpeedPoint> rabbitProfile;  // Time-tagged speeds
};

// Global pre-computed table
static std::vector<ScenarioVariations> gScenarioVariations;
static unsigned int gPathSeed;  // Derived from GPrand() or RandomPathSeedB override
```

#### Startup Pre-fetch

```cpp
// Called once at startup, after GPInit()
void prefetchAllVariations(int numScenarios, const VariationSigmas& sigmas,
                           const RabbitSpeedConfig& rabbitCfg) {
    gScenarioVariations.clear();
    gScenarioVariations.reserve(numScenarios);

    // Derive pathSeed (unless overridden by config)
    if (ConfigManager::getExtraConfig().randomPathSeedB == -1) {
        gPathSeed = static_cast<unsigned int>(GPrand());
    } else {
        gPathSeed = static_cast<unsigned int>(ConfigManager::getExtraConfig().randomPathSeedB);
    }

    for (int i = 0; i < numScenarios; i++) {
        ScenarioVariations sv;

        // Wind seed for crrcsim
        sv.windSeed = static_cast<unsigned int>(GPrand());

        // Entry variations (consume ~10 GPrand values)
        sv.entryOffsets = generateVariationsFromGPrand(sigmas);

        // Rabbit speed profile (consume ~50-200 GPrand values)
        sv.rabbitProfile = generateSpeedProfileFromGPrand(rabbitCfg, SIM_TOTAL_TIME_MSEC / 1000.0);

        gScenarioVariations.push_back(std::move(sv));
    }

    *logger.info() << "Pre-fetched " << numScenarios << " scenario variations from GPrand()" << std::endl;
}
```

### Phase 2: Modify variation_generator.h

Replace inline LCG with GPrand() consumption:

```cpp
// OLD (inline LCG):
inline VariationOffsets generateVariations(unsigned int seed, const VariationSigmas& sigmas) {
    auto nextDouble = [&seed]() -> double {
        seed = seed * 1103515245 + 12345;
        return static_cast<double>((seed >> 16) & 0x7FFF) / 32768.0;
    };
    // ...
}

// NEW (consume from GPrand):
inline VariationOffsets generateVariationsFromGPrand(const VariationSigmas& sigmas) {
    auto nextDouble = []() -> double {
        return static_cast<double>(GPrand()) / static_cast<double>(RAND_MAX);
    };

    auto gaussian = [&nextDouble](double sigma) -> double {
        double u1 = nextDouble() * 0.999 + 0.001;
        double u2 = nextDouble();
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        return z * sigma;
    };

    VariationOffsets v;
    v.entryHeadingOffset = gaussian(sigmas.headingSigma);
    v.entryRollOffset = gaussian(sigmas.rollSigma);
    v.entryPitchOffset = gaussian(sigmas.pitchSigma);
    v.entrySpeedFactor = 1.0 + gaussian(sigmas.speedSigma);
    v.windDirectionOffset = gaussian(sigmas.windDirectionSigma);
    return v;
}

// Similar for generateSpeedProfileFromGPrand()
```

### Phase 3: Remove Old Infrastructure

#### Config Changes

```ini
# REMOVE these (derive from GPSeed instead):
# WindSeedBase = 12121
# WindSeedStride = 1000

# KEEP:
GPSeed = 12345              # Master seed for everything
RandomPathSeedB = -1        # -1 = derive from GPrand, else override

# KEEP variation sigmas (control magnitude, not seeds):
EntryHeadingSigma = 5.0
EntryRollSigma = 3.0
# ... etc
```

#### Code Removals

1. **autoc.cc**: Remove `seedBase`, `seedStride` computation in `buildGenerationScenarios()`
2. **autoc.cc**: Remove `populateVariationOffsets()` - use pre-fetched table instead
3. **autoc.cc**: Remove `applySpeedProfileToPath()` - use pre-fetched profiles
4. **config_manager.cc**: Remove WindSeedBase, WindSeedStride parsing
5. **pathgen.cc**: Remove dead `randomPointInCylinder()` function

### Phase 4: Fix Demic Migration Bug

```cpp
// OLD (non-deterministic):
static thread_local std::mt19937 rng(std::random_device{}());

// NEW (deterministic):
// Use GPrand() for migration probability
if (static_cast<double>(GPrand()) / RAND_MAX > migProb) {
    continue;
}
```

### Phase 5: Update Scenario Application

```cpp
// When setting up a scenario for evaluation:
void applyScenarioVariations(int scenarioIndex, ScenarioMetadata& meta,
                             std::vector<Path>& path) {
    const auto& sv = gScenarioVariations[scenarioIndex];

    // Wind seed for crrcsim
    meta.windSeed = sv.windSeed;

    // Entry variations
    meta.entryHeadingOffset = sv.entryOffsets.entryHeadingOffset;
    meta.entryRollOffset = sv.entryOffsets.entryRollOffset;
    meta.entryPitchOffset = sv.entryOffsets.entryPitchOffset;
    meta.entrySpeedFactor = sv.entryOffsets.entrySpeedFactor;
    meta.windDirectionOffset = sv.entryOffsets.windDirectionOffset;

    // Apply rabbit speed profile to path timing
    applyRabbitProfileToPath(path, sv.rabbitProfile);
}
```

## Startup Logging

At startup, after pre-fetching all variations, log a consistent summary for verification and debugging.

### Log Format

```
=== Pre-fetched Scenario Variations (GPSeed=12345) ===
PathSeed: 98765432 (override: no)
Scenarios: 49

Scenario  WindSeed    Heading°   Roll°    Pitch°   Speed%   WindDir°  RabbitSpd(min/avg/max)
--------  ----------  --------   ------   ------   ------   --------  ----------------------
   0      0x12AB34CD    +2.31    -0.45    +0.12   +3.2%     -4.56     14.2 / 16.1 / 18.3
   1      0x45DE67FA    -1.87    +0.89    -0.34   -1.8%     +2.34     15.1 / 16.0 / 17.2
   ...
  48      0xABCD1234    +0.56    -1.23    +0.67   +0.4%     +1.23     13.8 / 15.9 / 19.1

Total GPrand() consumed for pre-fetch: 4732
```

### Implementation

```cpp
void logPrefetchedVariations(int numScenarios) {
    *logger.info() << "\n=== Pre-fetched Scenario Variations (GPSeed="
                   << ConfigManager::getExtraConfig().gpSeed << ") ===" << std::endl;
    *logger.info() << "PathSeed: " << gPathSeed
                   << " (override: " << (ConfigManager::getExtraConfig().randomPathSeedB != -1 ? "yes" : "no")
                   << ")" << std::endl;
    *logger.info() << "Scenarios: " << numScenarios << std::endl;
    *logger.info() << std::endl;

    *logger.info() << "Scenario  WindSeed    Heading°   Roll°    Pitch°   Speed%   WindDir°  RabbitSpd" << std::endl;
    *logger.info() << "--------  ----------  --------   ------   ------   ------   --------  ---------" << std::endl;

    for (int i = 0; i < numScenarios; i++) {
        const auto& sv = gScenarioVariations[i];
        double minSpd = sv.rabbitProfile.empty() ? 0 : sv.rabbitProfile[0].speed;
        double maxSpd = minSpd, sumSpd = 0;
        for (const auto& pt : sv.rabbitProfile) {
            minSpd = std::min(minSpd, pt.speed);
            maxSpd = std::max(maxSpd, pt.speed);
            sumSpd += pt.speed;
        }
        double avgSpd = sv.rabbitProfile.empty() ? 0 : sumSpd / sv.rabbitProfile.size();

        *logger.info() << std::setw(4) << i << "      "
                       << "0x" << std::hex << std::setw(8) << std::setfill('0') << sv.windSeed << std::dec << std::setfill(' ')
                       << "  " << std::setw(+7) << std::fixed << std::setprecision(2) << radToDeg(sv.entryOffsets.entryHeadingOffset)
                       << "  " << std::setw(6) << radToDeg(sv.entryOffsets.entryRollOffset)
                       << "  " << std::setw(6) << radToDeg(sv.entryOffsets.entryPitchOffset)
                       << "  " << std::setw(+5) << std::setprecision(1) << ((sv.entryOffsets.entrySpeedFactor - 1.0) * 100)
                       << "%  " << std::setw(6) << std::setprecision(2) << radToDeg(sv.entryOffsets.windDirectionOffset)
                       << "  " << std::setprecision(1) << minSpd << "/" << avgSpd << "/" << maxSpd
                       << std::endl;
    }
    *logger.info() << std::endl;
}
```

## Verification

### Determinism Test

```bash
# Run twice with same seed, compare outputs
./autoc --config test.ini --seed 12345 > run1.log 2>&1
./autoc --config test.ini --seed 12345 > run2.log 2>&1
diff run1.log run2.log  # Should be empty
```

### GPrand Consumption Audit

Add optional logging to track GPrand() consumption order:

```cpp
#ifdef DEBUG_PRNG
static int gGPrandCallCount = 0;
inline int TrackedGPrand() {
    int val = GPrand();
    *logger.debug() << "GPrand[" << gGPrandCallCount++ << "] = " << val << std::endl;
    return val;
}
#define GPrand() TrackedGPrand()
#endif
```

## Migration Path

1. Implement pre-fetch infrastructure (can coexist with old code)
2. Add feature flag: `UsePreFetchedVariations = 1`
3. When enabled, use pre-fetched values; when disabled, use old path
4. Verify identical behavior with flag disabled
5. Enable by default, remove old code

## PRNG Inventory After Refactor

| PRNG | Location | Purpose | Notes |
|------|----------|---------|-------|
| GPrand() | GP library | ALL autoc randomness | Single source |
| std::mt19937(pathSeed) | pathgen.h | Path geometry | Seeded from GPrand() or RandomPathSeedB |
| CRRC_Random | crrcsim | Wind/thermal/gust | Seeded from pre-fetched windSeeds[] |

**Rule**: autoc creates NO PRNGs. It only consumes from GPrand(). Seeds passed to external systems (pathgen, crrcsim) are pre-computed values derived from GPrand().

## Files Changed

| File | Changes |
|------|---------|
| autoc.cc | Add prefetch, remove seed computation, fix demic |
| autoc.h | Add ScenarioVariations struct |
| variation_generator.h | Add GPrand-based generators, deprecate LCG versions |
| config_manager.cc | Remove WindSeedBase/Stride parsing |
| autoc.ini | Remove WindSeedBase/Stride |
| pathgen.cc | Remove dead randomPointInCylinder |

## Risks

1. **GPrand() thread safety**: Verify GP library's rand is safe for multi-threaded eval
2. **Consumption order**: Pre-fetch must happen in deterministic order
3. **Scenario count changes**: If numScenarios changes between runs, variation assignments shift

## Success Criteria

1. Single `GPSeed` config controls all randomness
2. Identical runs with same GPSeed
3. No `std::random_device` or `time(NULL)` seeding anywhere in autoc
4. Demic migration is deterministic
5. Rabbit speed properly continues PRNG sequence (no restart bug)

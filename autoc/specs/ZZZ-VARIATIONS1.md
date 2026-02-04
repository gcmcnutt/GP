# VARIATIONS1: Entry and Wind Direction Variations

## Status: IN PROGRESS (autoc complete, crrcsim pending)

## Problem Statement

Current training uses a single fixed configuration:
- **Launch heading**: 180° (South) - hardcoded in `autoc_config.xml`
- **Wind direction**: 330° (from NNW) - hardcoded in `autoc_config.xml`
- **Entry attitude**: Level flight (0° roll, 0° pitch)
- **Entry speed**: 5.5 m/s (0.5 × V_ref)

Flight 25a/b revealed that real-world activations occur with:
- Heading ±50° off from path
- Roll ±45° (sometimes inverted)
- Varying speeds

The GP controller, trained only on ideal entry conditions, immediately saturates controls and diverges when faced with these off-nominal entries.

## Goal

Expand training scenarios to include:
1. **Entry heading variations** - Aircraft heading offset from path at activation
2. **Entry attitude variations** - Roll and pitch perturbations
3. **Entry speed variations** - Faster/slower than reference
4. **Wind direction variations** - Gaussian around configured base direction

## Design Principles

1. **Deterministic from seed** - All variations derived from a single PRNG seed per scenario
2. **Simple configuration** - Extend existing `WindScenarios` / `WindSeedBase` pattern
3. **Fast switching** - Avoid config file reloads; use runtime parameter overrides
4. **Backward compatible** - Empty/zero seed = legacy behavior

## Implementation

### Architecture: What Goes Where

Autoc generates variation offsets from sigma parameters and sends them to crrcsim via ScenarioMetadata. This is cleaner than having crrcsim regenerate them because:
- Single source of truth for sigma values (autoc.ini)
- Variations logged at startup match what's sent to sim
- crrcsim just applies values, no generation logic needed

```
┌─────────────────────────────────────────────────────────────────────┐
│                          AUTOC (GP trainer)                          │
├─────────────────────────────────────────────────────────────────────┤
│  autoc.ini:                                                          │
│    - WindScenarios = 36                                              │
│    - WindSeedBase = 12121                                            │
│    - EntryHeadingSigma = 45  (degrees)                               │
│    - EntryRollSigma = 22.5                                           │
│    - ... (other sigma values)                                        │
│                                                                      │
│  variation_generator.h:                                              │
│    - generateVariations(seed, sigmas) → VariationOffsets             │
│    - VariationSigmas::fromDegrees() converts config → radians        │
│                                                                      │
│  autoc.cc:                                                           │
│    - Startup: read sigmas, log variation table                       │
│    - Scenario creation: populateVariationOffsets(meta)               │
│      fills in meta.entryHeadingOffset, etc.                          │
│                                                                      │
│  ScenarioMetadata (minisim.h):                                       │
│    - windSeed (for turbulence PRNG)                                  │
│    - entryHeadingOffset, entryRollOffset, etc. (radians)             │
│    - entrySpeedFactor (multiplier)                                   │
│    - windDirectionOffset (radians)                                   │
└─────────────────────────────────────────────────────────────────────┘
                              │
                              │ ScenarioMetadata (with variation offsets)
                              ▼
┌─────────────────────────────────────────────────────────────────────┐
│                         CRRCSIM (simulator)                          │
├─────────────────────────────────────────────────────────────────────┤
│  autoc_interface.cpp - receives ScenarioMetadata with offsets        │
│    - No generation needed, just read the values                      │
│                                                                      │
│  aircraft.cpp - apply entry variations at launch:                    │
│    - launchHeading += meta.entryHeadingOffset                        │
│    - initialRoll = meta.entryRollOffset                              │
│    - initialPitch = meta.entryPitchOffset                            │
│    - launchVelocity *= meta.entrySpeedFactor                         │
│                                                                      │
│  windfield.cpp - apply wind direction offset:                        │
│    - windDir = baseWindDir + meta.windDirectionOffset                │
│                                                                      │
│  These offsets REPLACE defaults from autoc_config.xml                │
└─────────────────────────────────────────────────────────────────────┘
```

### ScenarioMetadata - Extended for VARIATIONS1

The `ScenarioMetadata` struct now includes variation offset fields (version 5):

```cpp
// minisim.h
struct ScenarioMetadata {
  int pathVariantIndex = -1;
  int windVariantIndex = -1;
  unsigned int windSeed = 0;
  uint64_t scenarioSequence = 0;
  uint64_t bakeoffSequence = 0;
  bool enableDeterministicLogging = false;

  // VARIATIONS1: Entry and wind direction offsets (radians, computed by autoc)
  double entryHeadingOffset = 0.0;   // offset from path tangent
  double entryRollOffset = 0.0;      // initial roll attitude
  double entryPitchOffset = 0.0;     // initial pitch attitude
  double entrySpeedFactor = 1.0;     // multiplier on reference speed
  double windDirectionOffset = 0.0;  // offset from base wind direction
  // ... serialization (version 5)
};
```

Values are computed by autoc and serialized to crrcsim. crrcsim just reads and applies them.

### Variation Generator (in autoc)

The `variation_generator.h` is used by autoc only. It generates variations from a seed and configurable sigma parameters:

```cpp
#ifndef VARIATION_GENERATOR_H
#define VARIATION_GENERATOR_H

#include <cmath>

// Sigma parameters (in radians, except speed which is a fraction)
struct VariationSigmas {
    double headingSigma;       // radians
    double rollSigma;          // radians
    double pitchSigma;         // radians
    double speedSigma;         // fraction (0.1 = ±10% at 1σ)
    double windDirectionSigma; // radians

    // Convert from degrees (config) to radians (internal)
    static VariationSigmas fromDegrees(double headingDeg, double rollDeg, double pitchDeg,
                                        double speedFrac, double windDirDeg) {
        return VariationSigmas{
            headingDeg * M_PI / 180.0,
            rollDeg * M_PI / 180.0,
            pitchDeg * M_PI / 180.0,
            speedFrac,
            windDirDeg * M_PI / 180.0
        };
    }
};

// Derived variation values - computed from windSeed
struct VariationOffsets {
    double entryHeadingOffset;  // radians
    double entryRollOffset;     // radians
    double entryPitchOffset;    // radians
    double entrySpeedFactor;    // multiplier (1.0 = nominal)
    double windDirectionOffset; // radians
};

// Generate all variations from seed and sigma parameters
// Uses LCG PRNG with Box-Muller for Gaussian sampling
inline VariationOffsets generateVariations(unsigned int seed, const VariationSigmas& sigmas) {
    VariationOffsets v;

    // LCG PRNG (deterministic from seed)
    auto nextDouble = [&seed]() -> double {
        seed = seed * 1103515245 + 12345;
        return (double)((seed >> 16) & 0x7FFF) / 32768.0;
    };

    // Box-Muller transform for Gaussian sampling
    auto gaussian = [&nextDouble](double sigma) -> double {
        double u1 = nextDouble() * 0.999 + 0.001;  // avoid log(0)
        double u2 = nextDouble();
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        return z * sigma;
    };

    v.entryHeadingOffset = gaussian(sigmas.headingSigma);
    v.entryRollOffset = gaussian(sigmas.rollSigma);
    v.entryPitchOffset = gaussian(sigmas.pitchSigma);
    v.entrySpeedFactor = 1.0 + gaussian(sigmas.speedSigma);
    v.windDirectionOffset = gaussian(sigmas.windDirectionSigma);

    return v;
}

#endif // VARIATION_GENERATOR_H
```

**Gaussian distribution rationale:**
- ~68% of scenarios are within ±1σ of nominal (most training near expected conditions)
- ~95% are within ±2σ (good coverage of off-nominal)
- ~0.3% are beyond ±3σ (rare extreme cases for robustness)

### Phase 3: crrcsim Integration

In `autoc_interface.cpp`, when receiving a new scenario, the variation offsets are already
computed and present in ScenarioMetadata. crrcsim just reads and applies them:

```cpp
void AutocInterface::applyScenario(const ScenarioMetadata& meta) {
    // No generation needed - values already in metadata
    // Log for debugging (optional, controlled by verbosity)
    if (verbose) {
        printf("Scenario seed=%u: heading=%+.1f° roll=%+.1f° pitch=%+.1f° speed=%.2fx wind=%+.1f°\n",
               meta.windSeed,
               meta.entryHeadingOffset * 180.0/M_PI,
               meta.entryRollOffset * 180.0/M_PI,
               meta.entryPitchOffset * 180.0/M_PI,
               meta.entrySpeedFactor,
               meta.windDirectionOffset * 180.0/M_PI);
    }

    // Store metadata for use by aircraft launch and wind init
    currentMeta = meta;
}
```

**Launch parameters** (in `aircraft.cpp` or launch reset):
```cpp
// Apply entry variations from ScenarioMetadata (all in radians)
launchHeading += meta.entryHeadingOffset;
initialRoll = meta.entryRollOffset;
initialPitch = meta.entryPitchOffset;
launchVelocity *= meta.entrySpeedFactor;
```

**Wind direction** (in `windfield.cpp` initialization):
```cpp
// Apply wind direction offset from ScenarioMetadata
double baseWindDir = cfg->wind->getDirection();  // from autoc_config.xml
double newWindDir = baseWindDir + meta.windDirectionOffset;
// Normalize to [0, 2π)
while (newWindDir < 0) newWindDir += 2.0 * M_PI;
while (newWindDir >= 2.0 * M_PI) newWindDir -= 2.0 * M_PI;
cfg->wind->setDirection(newWindDir);
```

### Phase 4: Configuration

See [autoc.ini](../autoc.ini) for complete configuration. Key changes from baseline:

```ini
# 36 scenarios with Gaussian-distributed variations
WindScenarios = 36
WindSeedBase = 12121
WindSeedStride = 1          # was 11, now 1 for denser seed coverage

# Enable variations (requires crrcsim support, 1=enabled, 0=disabled)
EnableEntryVariations = 1
EnableWindVariations = 1

# Gaussian σ values (degrees in config, converted to radians in variation_generator.h)
EntryHeadingSigma = 45      # 68% within ±45°, 95% within ±90°
EntryRollSigma = 22.5       # 68% within ±22.5°, 95% within ±45°
EntryPitchSigma = 7.5       # 68% within ±7.5°, 95% within ±15°
EntrySpeedSigma = 0.1       # 68% within ±10%, 95% within ±20%
WindDirectionSigma = 45     # around base (330° from autoc_config.xml)

# Expanded node set (19 nodes, adds state feedback)
TrainingNodes = ...,GETROLL_RAD,GETPITCH_RAD,GETVEL
```

**Note**: The sigma values in autoc.ini are read by autoc, converted to radians, and passed to `generateVariations(seed, sigmas)`. The generated offsets are then stored in ScenarioMetadata and sent to crrcsim. This ensures a single source of truth for sigma configuration.

## Scenario Matrix

With **36 scenarios** (up from 9), each seed produces unique Gaussian-distributed variations:

| Seed | Heading (σ=45°) | Roll (σ=22.5°) | Pitch (σ=7.5°) | Speed (σ=10%) | Wind (σ=45°) |
|------|-----------------|----------------|----------------|---------------|--------------|
| 0    | +0.4 rad (+23°) | -0.2 rad (-12°)| +0.05 rad (+3°)| 1.05×         | +0.3 rad     |
| 1    | -0.8 rad (-46°) | +0.5 rad (+29°)| -0.1 rad (-6°) | 0.93×         | -0.9 rad     |
| 2    | +0.1 rad (+6°)  | -0.6 rad (-34°)| +0.2 rad (+11°)| 1.12×         | +1.3 rad     |
| ...  | (Gaussian)      | (Gaussian)     | (Gaussian)     | (Gaussian)    | (Gaussian)   |

**Gaussian distribution properties:**
- ~68% of scenarios within ±1σ of nominal (near-ideal conditions)
- ~95% of scenarios within ±2σ (moderate off-nominal)
- ~5% of scenarios beyond ±2σ (challenging edge cases)

This focuses training on typical conditions while ensuring robustness to extremes.

**Why 36?** Square numbers (25, 36, 49) provide good coverage. 36 = 6² gives 4× more scenarios than current 9. Can scale to 49 if compute allows, or drop to 25 for faster iteration.

## Fitness Function Considerations

Current fitness aggregates across scenarios. With more extreme entry conditions:

1. **Early divergence expected** - First 2-3 seconds may have high errors during recovery
2. **Recovery bonus** - Consider rewarding successful recovery from bad entry
3. **Weighted scenarios** - Could weight "easy" entries higher initially, then anneal

For first pass, keep fitness unchanged. The GP should naturally evolve recovery behaviors if the scenarios demand them.

## Training Implications

### Search Space Impact

Adding entry variations dramatically expands the state space the GP must handle:
- Current: ~1 entry condition × 9 wind patterns = 9 effective scenarios
- Proposed: ~∞ entry conditions × ~∞ wind directions = 36 seed-sampled scenarios

Node set grows from 16 to 19 nodes (+19% tree building options).

This argues for:
- More generations (current 100 may be insufficient, try 200-300)
- Larger population (current 256 may need 512+)
- Or: staged training (first master calm entries, then add perturbations)

### Node Set Expansion

The minimal node set (`GETDPHI, GETDTHETA, GETDTARGET`) is insufficient for recovery from off-nominal entries:
- No direct attitude feedback (`GETROLL_RAD`, `GETPITCH_RAD` excluded)
- Relies on target-relative angles which conflate attitude and position error
- No velocity information for energy management

**New node set for VARIATIONS1:**
```ini
# autoc-variations1.ini
TrainingNodes = SETPITCH,SETROLL,SETTHROTTLE,GETDPHI,GETDTHETA,GETDTARGET,ADD,SUB,MUL,DIV,IF,GT,0,1,2,PROGN,GETROLL_RAD,GETPITCH_RAD,GETVEL
```

Added nodes:
- `GETROLL_RAD` - Aircraft roll angle in radians (from quaternion)
- `GETPITCH_RAD` - Aircraft pitch angle in radians (from quaternion)
- `GETVEL` - Airspeed magnitude (m/s) - relative to air, accounts for wind

This brings the training set from 16 to 19 nodes. The search space grows but attitude feedback is essential for recovery behaviors.

## Migration Path

### Step 1: Update autoc (COMPLETE)
- autoc.ini: WindScenarios=36, sigma params, TrainingNodes ✓
- autoc.h: ExtraConfig fields for variations ✓
- config_manager.cc: Parse sigma params ✓
- variation_generator.h: NEW header-only variation generator ✓
- autoc.cc: Include variation_generator.h, log table, populate offsets ✓
- minisim.h: ScenarioMetadata version 5 with offset fields ✓
- CMakeLists.txt: Add variation_generator.h ✓

### Step 2: Wire up crrcsim launch variations (TODO)
- In aircraft launch code, read and apply heading/roll/pitch/speed offsets from ScenarioMetadata
- Add console logging for debugging

### Step 3: Wire up crrcsim wind variations (TODO)
- In windfield init, read and apply wind direction offset from ScenarioMetadata
- Add console logging for debugging

### Step 4: Test with short run
- Run 10 generations
- Verify renderer shows varied entry attitudes
- Verify varied wind directions in sim

### Step 5: Full training
- Run 200+ generations with 36 scenarios
- Compare to baseline (9 scenarios, no variations)
- Analyze evolved programs for recovery patterns

## Archiving Variation Details

Currently best-of-generation writes to `data.dat` with fitness and GP tree. We need to also record which scenarios (and their variation values) were used.

### Option A: Log variation table once at start

At autoc startup, compute and log all 36 scenario variations:

```
# data.dat header (or separate variations.log)
# SCENARIO VARIATIONS (seed base=12121, stride=1)
# Seed     HeadingOff  RollOff   PitchOff  SpeedFac  WindOff
  12121    +0.412      -0.198    +0.052    1.048     +0.287
  12122    -0.834      +0.512    -0.103    0.927     -0.912
  12123    +0.087      -0.623    +0.189    1.118     +1.341
  ...
```

This is logged once and referenced by windVariantIndex in per-generation output.

### Option B: Log per best-of-generation

Each time we write best-of-generation, also write the scenario breakdown:

```
Best of generation 5 (Fitness = 148326.5, Structural Complexity = 42)
  Scenario breakdown:
    wind[0] seed=12121 heading=+23.6° roll=-11.3° pitch=+3.0° speed=1.05x wind=+16.4°: fitness=12453.2
    wind[1] seed=12122 heading=-47.8° roll=+29.3° pitch=-5.9° speed=0.93x wind=-52.3°: fitness=18921.1
    ...
```

### Recommendation

Use **Option A** (log once at start) for efficiency. The variation values are deterministic from seed, so we only need to record the mapping once. Per-generation output just references windVariantIndex.

Add to autoc.cc at startup:
```cpp
void logVariationTable(const ExtraConfig& cfg, const VariationSigmas& sigmas) {
    cout << "SCENARIO VARIATIONS (base=" << cfg.windSeedBase
         << " count=" << cfg.windScenarioCount << ")" << endl;
    cout << "  Sigmas: heading=" << (sigmas.headingSigma * 180.0/M_PI) << "°"
         << " roll=" << (sigmas.rollSigma * 180.0/M_PI) << "°"
         << " pitch=" << (sigmas.pitchSigma * 180.0/M_PI) << "°"
         << " speed=" << (sigmas.speedSigma * 100.0) << "%"
         << " wind=" << (sigmas.windDirectionSigma * 180.0/M_PI) << "°" << endl;
    cout << "# Idx  Seed     Heading   Roll      Pitch     Speed   WindDir" << endl;

    for (int i = 0; i < cfg.windScenarioCount; i++) {
        unsigned int seed = cfg.windSeedBase + i * cfg.windSeedStride;
        VariationOffsets v = generateVariations(seed, sigmas);
        cout << setw(4) << i << "  " << seed
             << "  " << fixed << setprecision(1)
             << setw(8) << (v.entryHeadingOffset * 180.0/M_PI) << "°"
             << setw(8) << (v.entryRollOffset * 180.0/M_PI) << "°"
             << setw(8) << (v.entryPitchOffset * 180.0/M_PI) << "°"
             << setw(8) << v.entrySpeedFactor << "x"
             << setw(8) << (v.windDirectionOffset * 180.0/M_PI) << "°"
             << endl;
    }
}
```

Sigma values are configured in autoc.ini and passed to `generateVariations()`. The computed offsets are stored in ScenarioMetadata and sent to crrcsim - no generation happens on the crrcsim side.

## Files to Modify

### autoc (trainer) - COMPLETE

| File | Changes | Status |
|------|---------|--------|
| `autoc/autoc.ini` | WindScenarios=36, sigma params, TrainingNodes | ✓ |
| `autoc/autoc.h` | ExtraConfig fields for variations | ✓ |
| `autoc/config_manager.cc` | Parse EnableEntryVariations, sigma params | ✓ |
| `autoc/variation_generator.h` | NEW: header-only variation generator | ✓ |
| `autoc/autoc.cc` | Include variation_generator.h, log table, populate offsets | ✓ |
| `autoc/minisim.h` | ScenarioMetadata version 5 with offset fields | ✓ |
| `autoc/CMakeLists.txt` | Add variation_generator.h to autoc_common | ✓ |

### crrcsim (simulator) - TODO

| File | Changes | Status |
|------|---------|--------|
| `crrcsim/src/aircraft.cpp` | Apply entry variations from ScenarioMetadata at launch | |
| `crrcsim/src/windfield.cpp` | Apply wind direction offset from ScenarioMetadata | |
| `crrcsim/src/mod_cntrl_autoc/autoc_interface.cpp` | Log received variation values (debug) | |

Note: crrcsim does NOT need variation_generator.h - it just reads pre-computed values from ScenarioMetadata.

## Success Criteria

1. **Training completes** - 200 generations without crashes
2. **Fitness improves** - Downward trend despite harder scenarios
3. **Recovery behaviors emerge** - Evolved programs handle 45° entry rolls
4. **Flight test** - Next real flight shows improved robustness

## Risks

1. **Search space explosion** - Too many dimensions may slow convergence
   - Mitigation: 36 scenarios is 4× current; can drop to 25 if too slow

2. **Conflicting objectives** - Optimal for headwind ≠ optimal for tailwind
   - Mitigation: Aggregated fitness should find compromise

3. **Simulation fidelity** - If sim doesn't match reality, variations won't help
   - Mitigation: Flight 25a analysis shows sim dynamics are reasonable (10.1 vs 11.0 m/s cruise)

## Timeline

Most autoc work is complete. Remaining work is in crrcsim (reading and applying offset values).

1. **Phase 1**: Config ready (autoc.ini updated) ✓
2. **Phase 2**: Add variation_generator.h to autoc ✓
3. **Phase 3**: Update ScenarioMetadata with offset fields ✓
4. **Phase 4**: Populate offsets in autoc scenario creation ✓
5. **Phase 5**: Wire up launch variations in crrcsim aircraft.cpp
6. **Phase 6**: Wire up wind direction in crrcsim windfield.cpp
7. **Phase 7**: Test run (10 generations), verify in renderer
8. **Phase 8**: Full training (200+ generations)
9. **Phase 9**: Flight test

## References

- [ZZZ-NODE_MASKING.md](ZZZ-NODE_MASKING.md) - Minimal node set rationale
- [LAYERED_CONTROLLER.md](LAYERED_CONTROLLER.md) - Safety layer design (handles OOB, not recovery)
- [TODO](../TODO) - Variations backlog items (lines 29-58)
- [FLIGHT25_ANALYSIS.md](/home/gmcnutt/xiao-gp/postflight/FLIGHT25_ANALYSIS.md) - Entry condition analysis

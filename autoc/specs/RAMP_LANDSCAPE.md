# RAMP_LANDSCAPE: Gradual Variation Scaling

## Problem

Entry and rabbit speed variations create a fundamentally different fitness landscape than wind variations. With entry variations enabled (even small sigmas like 5°), GP struggles to make progress because:

1. Entry variations create transient recovery problems that dominate fitness
2. The fitness landscape becomes corrupted - selection pressure shifts to "don't crash on worst entry" rather than "track better"
3. Average fitness is dominated by hardest scenarios

Wind works because it's a continuous perturbation with preserved gradients. Entry/rabbit change the initial problem.

## Solution

Gradually ramp variation sigmas from 0 to full value over the course of training. This allows:

1. Early generations to establish good tracking behavior (easy problem)
2. Progressive introduction of harder scenarios
3. Existing tracking skills to adapt incrementally to entry recovery

## New Configuration Parameter

```ini
# Variation ramp: generations per step (0 = disabled, use full sigmas immediately)
# Linear ramp from 0.0 to 1.0 scale over NumberOfGenerations
# Example: 100 generations, VariationRampStep=5 -> 20 steps, each +5% scale
VariationRampStep               = 5
```

## Behavior

### Scale Calculation

```
current_generation = 0..NumberOfGenerations-1
scale = min(1.0, current_generation / NumberOfGenerations)
```

The scale is applied at the start of each generation (or every `VariationRampStep` generations for stepped behavior).

### Stepped vs Continuous

With `VariationRampStep=5`:
- Gen 0-4: scale = 0.00
- Gen 5-9: scale = 0.05
- Gen 10-14: scale = 0.10
- ...
- Gen 95-99: scale = 0.95

With `VariationRampStep=1` (continuous):
- Gen 0: scale = 0.00
- Gen 1: scale = 0.01
- Gen 2: scale = 0.02
- ...
- Gen 99: scale = 0.99

With `VariationRampStep=0` (disabled):
- All generations: scale = 1.0 (current behavior)

### Applied Sigmas

The scale multiplies all variation sigmas before scenario generation:

```cpp
double scale = computeVariationScale(currentGeneration, totalGenerations, rampStep);

// Applied values sent to simulator
double appliedHeadingSigma = EntryHeadingSigma * scale;
double appliedRollSigma = EntryRollSigma * scale;
double appliedPitchSigma = EntryPitchSigma * scale;
double appliedSpeedSigma = EntrySpeedSigma * scale;
double appliedRabbitSpeedSigma = RabbitSpeedSigma * scale;
// Wind can optionally be ramped too, or left at full
```

## Implementation

### 1. Config Loading (autoc.cc)

Add to config parsing:
```cpp
int variationRampStep = cfg.value("VariationRampStep", 0);
```

### 2. Scale Computation

```cpp
double computeVariationScale(int currentGen, int totalGens, int rampStep) {
    if (rampStep <= 0) return 1.0f;  // Disabled

    // Stepped: quantize to step boundaries
    int stepIndex = currentGen / rampStep;
    int totalSteps = totalGens / rampStep;
    if (totalSteps == 0) return 1.0f;

    return std::min(1.0f, static_cast<double>(stepIndex) / static_cast<double>(totalSteps));
}
```

### 3. Apply to Scenario Generation

In the scenario prefetch or per-generation scenario setup, multiply sigmas by scale before generating Gaussian samples.

### 4. Logging

Log the current scale at generation boundaries:
```
<info> Gen 10: VariationScale=0.10 (step 2/20)
```

## Example Configuration

```ini
# Train with entry variations, ramped over 100 generations
NumberOfGenerations             = 100
VariationRampStep               = 5

EnableEntryVariations           = 1
EntryHeadingSigma               = 45
EntryRollSigma                  = 22.5
EntryPitchSigma                 = 7.5
EntrySpeedSigma                 = 10

RabbitSpeedSigma                = 2.0
```

At gen 0: all sigmas effectively 0 (nominal entries)
At gen 50: sigmas at 50% (heading ±22.5°, roll ±11.25°, etc.)
At gen 100: sigmas at full value

## Future Extensions

1. **Per-parameter ramps**: Different ramp rates for heading vs roll vs rabbit speed
2. **Adaptive ramp**: Increase difficulty only when fitness plateaus
3. **Non-linear ramps**: Sigmoid or exponential curves
4. **Wind ramp control**: Separate flag to include/exclude wind from ramping

## Testing

1. Run with `VariationRampStep=0` - should match current behavior exactly
2. Run with `VariationRampStep=5` on entry variations - should show smooth progression
3. Compare gen-50 fitness to standalone run at 50% sigmas
4. Verify final controller handles full variations

## Files to Modify

- `autoc.ini` - add VariationRampStep parameter
- `autoc.cc` - config loading, scale computation, apply to scenario generation
- Possibly `autoc.h` - if scale needs to be passed to evaluation context

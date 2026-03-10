# Quickstart: Variations Redux

**Date**: 2026-03-09
**Feature**: 003-variations-redux

## Prerequisites

- GP repo on branch `003-variations-redux`
- crrcsim repo on branch `003-variations-redux`
- Both repos build successfully (perf mode recommended for training runs)

## Build

```bash
# GP/autoc (perf mode for training)
cd ~/GP/autoc && bash rebuild-perf.sh

# crrcsim (perf mode)
cd ~/crsim/crrcsim-0.9.13
mkdir -p build && cd build
cmake -DPERFORMANCE_BUILD=ON ..
make
```

## Phase 1: Tightened Distance (no variations)

1. Edit `autoc/autoc.h`: DISTANCE_NORM=2.0, DISTANCE_POWER=2.0
2. Ensure `autoc.ini` has: EnableEntryVariations=0, EnableWindVariations=0, WindScenarios=1
3. Build and run:
   ```bash
   cd ~/GP && make && ./build/autoc 2>&1 | tee log/autoc-003-p1.log
   ```
4. Verify: fitness converges, tracking distance < 10m on straights

## Phase 2: Wind Variations

1. Edit `autoc.ini`: EnableWindVariations=1, WindScenarios=9
2. Run training:
   ```bash
   ./build/autoc 2>&1 | tee log/autoc-003-p2.log
   ```
3. Verify: crrcsim logs show varied wind directions, fitness trends downward

## Phase 3: Entry + Wind Variations

1. Edit `autoc.ini`: EnableEntryVariations=1, EnableWindVariations=1, WindScenarios=36
2. Run training:
   ```bash
   ./build/autoc 2>&1 | tee log/autoc-003-p3.log
   ```
3. Verify: aircraft launches with varied attitudes, recovery behaviors emerge

## Phase 4: Variable Rabbit Speed

1. Edit `autoc.ini`: RabbitSpeedSigma=2.0
2. Run training and verify speed-adaptive throttle control

## Phase 5: Ramp Tuning

1. Edit `autoc.ini`: VariationRampStep=5 (or tune as needed)
2. Compare ramped vs non-ramped convergence

## Validation

- Check `data.dat` for fitness trend
- Check `data.stc` for GP tree structure (look for GETROLL_RAD/GETPITCH_RAD/GETVEL usage)
- Use renderer to visualize tracking: `./build/renderer -k <keyname>`

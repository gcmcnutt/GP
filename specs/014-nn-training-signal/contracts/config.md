# Contract: Configuration (autoc.ini)

## Overview

Flat key=value INI format parsed by inih. No sections. Comments with `#` and `;`.

## Keys (surviving from GP era + new)

### Evolution Core
```ini
PopulationSize = 50          # For CMA-ES; was 5000 for GA
NumberOfGenerations = 500
EvalThreads = 4
```

### NN Configuration
```ini
NNMutationSigma = 0.3        # Initial sigma (also CMA-ES sigma_0)
NNSigmaFloor = 0.05          # NEW: minimum sigma (0 = disabled)
NNCrossoverAlpha = -1         # GA only, ignored by CMA-ES
NNWeightFile = nn_weights.dat
NNInitMethod = xavier
```

### Optimizer
```ini
OptimizerType = sep-cma-es   # NEW: "ga" or "sep-cma-es"
```

### Curriculum
```ini
CurriculumEnabled = 0                    # NEW: 0/1
CurriculumSchedule = 1:50,7:150,49:0    # NEW: scenarios:until_gen,...  (0=forever)
```

### Fitness
```ini
FitnessAggregation = sum     # NEW: "sum", "minimax", or "percentile"
FitnessPercentile = 0.95     # NEW: only used with percentile mode
```

### Simulation
```ini
SimNumPathsPerGeneration = 6
PathGeneratorMethod = random
```

### S3
```ini
S3Bucket = my-bucket
S3Profile = default
```

### Output
```ini
OutputDir = runs/             # NEW: auto-created run subdirectory base
```

## Constraints

- All keys are flat (no sections, no nesting)
- Types: int, double, string only
- Unknown keys are silently ignored (forward compatibility)
- Missing keys use compiled defaults

# NODE_MASKING: Restricting Training Node Set While Preserving Evaluation

## Status: IMPLEMENTED ✓

The TrainingNodes feature is complete and verified working.

## Goal

- **Evaluation**: All nodes always available (backward compatibility for existing programs)
- **Training**: Only a subset of nodes used for tree creation and mutation
- **Configuration**: Mask defined in `autoc.ini`, easy to experiment with different subsets

## Implementation

### Configuration

```ini
# autoc.ini
# Comma-separated list of node names to INCLUDE in training
# If empty or missing, all nodes are included (backward compatible)
TrainingNodes = ADD,SUB,MUL,DIV,IF,GT,SETPITCH,SETROLL,SETTHROTTLE,GETDPHI,GETDTHETA,GETDTARGET,0,1,2,PROGN
```

### Key Files Modified

- **autoc-eval.cc**: `createNodeSet()` refactored to use NodeDef table and dynamic GPNodeSet sizing
- **autoc.cc**: Calls `setTrainingNodesMask()` before `initializeSimGP()`
- **autoc.h**: Added `trainingNodes` field to ExtraConfig, declared `setTrainingNodesMask()`
- **config_manager.cc**: Added TrainingNodes config parameter

### How It Works

1. Config loading reads `TrainingNodes` string from ini file
2. Before GP initialization, `setTrainingNodesMask(mask)` is called
3. `createNodeSet()` parses the mask and creates GPNodeSet with exact size needed
4. Only masked nodes are added to the node set
5. Empty mask = all 38 nodes (backward compatible)

### Architecture

```
autoc.cc::main()
    ↓
setTrainingNodesMask(extraCfg.trainingNodes)
    ↓
initializeSimGP() → createNodeSet(adfNs)
    ↓
GPNodeSet with N nodes (N = mask size or 38 if empty)
    ↓
Used by: Tree creation, Mutation, Evaluation
```

### GPNodeSet Sizing Fix

The gpc++ library iterates through all `containerSize()` slots expecting valid nodes.
When filtering to a subset, we must allocate GPNodeSet with the exact count needed:

```cpp
// Two-pass approach in createNodeSet():
// Pass 1: Count filtered nodes
int nodeCount = 0;
for (int i = 0; i < allNodesCount; i++) {
    if (!filterEnabled || mask.count(allNodes[i].name) > 0) {
        nodeCount++;
    }
}

// Allocate exact size
GPNodeSet& ns = *new GPNodeSet(nodeCount);

// Pass 2: Add filtered nodes
for (int i = 0; i < allNodesCount; i++) {
    if (!filterEnabled || mask.count(allNodes[i].name) > 0) {
        ns.putNode(...);
    }
}
```

## Tested Node Subsets

### Minimal (16 nodes) - Verified Working

```ini
# autoc-minimal.ini
TrainingNodes = SETPITCH,SETROLL,SETTHROTTLE,GETDPHI,GETDTHETA,GETDTARGET,ADD,SUB,MUL,DIV,IF,GT,0,1,2,PROGN
```

**Results (autoc_202minimal.log):**
- Gen 0 → Gen 98: Fitness improved from 6,217,507 to 148,326
- Outperformed full 38-node run (~154,000 plateau) with smaller search space
- Elite evolved interesting patterns like nested GETDPHI for adaptive lookahead

**Parameters tuned for minimal set:**
- TournamentSize = 7 (classic GP literature)
- MaximumDepthForCreation = 8
- MaximumDepthForCrossover = 10
- SwapMutationProbability = 10
- ShrinkMutationProbability = 10

**Observations:**
- Reduced search space finds competitive solutions faster
- Programs evolve adaptive gain via nested GETDPHI(GETDPHI(...))
- Cross-axis coupling emerges naturally (throttle depends on heading error)
- No altitude awareness by design - safety layer would handle that separately

### Full (38 nodes) - Baseline

```ini
TrainingNodes = ADD,SUB,MUL,DIV,IF,EQ,GT,SETPITCH,SETROLL,SETTHROTTLE,GETPITCH,GETROLL,GETTHROTTLE,SIN,COS,OP_PI,0,1,2,PROGN,GETDPHI,GETDTHETA,GETDTARGET,GETVEL,GETDHOME,GETALPHA,GETBETA,GETVELX,GETVELY,GETVELZ,GETROLL_RAD,GETPITCH_RAD,CLAMP,ATAN2,ABS,SQRT,MIN,MAX
```

### Excluded (rarely useful)

These nodes have never appeared in successful controllers:

- GETVELX, GETVELY, GETVELZ (redundant with GETVEL, GETDPHI, GETDTHETA)
- GETALPHA, GETBETA (aerodynamic angles - too low level)
- GETPITCH, GETROLL, GETTHROTTLE (readback of own commands - rarely useful)

## Verification

Verified with identical GP seeds:
1. Empty TrainingNodes → full 38-node evolution (autoc_200.log)
2. Explicit full list → identical evolution (autoc_201.log)

Both runs produced identical fitness values, confirming the feature works correctly.

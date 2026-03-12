# 013: Neuroevolution Controller

## Status: DRAFT
## Priority: NEXT
## Dependencies: 012-distance-temporal-nodes (complete), 005-entry-fitness-ramp (complete)
## Absorbs: 007-unify-eval (eval pipeline unification is prerequisite work within this feature)

---

## Motivation

### The GP Representation Ceiling

After implementing distance temporal nodes (012), training runs consistently show the GP
converging to **bang-bang control** — using `IF GETDIST_RATE 0 1` switching rather than
proportional distance regulation. At gen 85/200, GETDIST_RATE appears 13 times in the
elite tree, all in GT/IF switching patterns. Mean distance settles at ~22m vs the 7.5m
target, with fitness plateauing.

This is a structural limitation of tree-based GP:

1. **Node competition**: Each new sensor (GETDIST, GETDIST_PREV, GETDIST_RATE) competes
   for tree positions with all other nodes. The current 27-node training set means each
   tree position has a 1/27 chance of being any particular sensor. Adding cross-track
   sensors would push this to 33+, diluting search efficiency for all nodes.

2. **Cross-term discovery**: A proportional throttle controller needs `distance * gain`,
   but GP must evolve `MUL(GETDIST, constant)` through random tree construction. Bang-bang
   (`IF GETDIST_RATE 0 1`) is structurally simpler and dominates early.

3. **Temporal state explosion**: Every sensor wants `_PREV` and `_RATE` variants. We went
   from 2 sensors → 6 nodes (dPhi, dTheta + PREV/RATE each), then distance added 3 more.
   Cross-track would add 3-6 more. Each is an O(n) tax on search space.

4. **No simultaneous perception**: A GP tree evaluates one subtree at a time. It cannot
   naturally "see all sensors at once" the way a neural network hidden layer does. The GP
   must evolve complex nested structures to combine multiple signals.

### What a Neural Network Solves

A fixed-topology NN with weight evolution:

| Problem | GP Tree | Neural Network |
|---------|---------|----------------|
| Sensor access | 1 node per tree position, competing | All inputs always available |
| Cross-terms | Must evolve MUL(A, B) explicitly | Learned in hidden layer weights |
| Temporal state | Explicit _PREV/_RATE nodes per sensor | Recurrent hidden state (LSTM/GRU) |
| Adding sensors | Combinatorial search space expansion | +1 input neuron, minimal impact |
| Proportional control | Rare — bang-bang dominates | Natural — weighted sum |
| Deployment size | Variable tree depth, ~100-300 nodes | Fixed weight count, predictable |

### Relationship to 4D Fitness Surface

The backlog item `[DEFERRED] 4D Positional Fitness Surface` proposed decomposing distance
into along-track, cross-track, and vertical components with directional penalties. Analysis
during 012 implementation revealed that the GP **already has** the sensor information for
this (GETDPHI = lateral bearing, GETDTHETA = vertical bearing, GETDIST = range), but cannot
efficiently combine angle × distance to get cross-track meters.

Rather than encoding directional assumptions into the fitness function (which risks
over-constraining the search), providing a NN controller that can naturally compute
cross-terms from existing sensors is a cleaner solution. The NN can discover that
`lateral_crosstrack ≈ sin(dPhi) * dist` without us encoding it.

The fitness function stays simple (scalar distance + attitude), and the controller
representation becomes powerful enough to exploit directional information.

### Relationship to Existing Backlog Items

| Backlog Item | Disposition |
|---|---|
| `[DEFERRED] 4D Positional Fitness Surface` | Superseded — NN can discover directional control without fitness encoding |
| `[DEFERRED] Control Command Smoothness Fitness Objective` | May become unnecessary — NN outputs tend to be smoother than GP trees due to continuous weight space |
| `[DEFERRED] Co-evolve Random Paths` | Orthogonal — compatible with neuroevolution |
| `[DEFERRED] LLM Fitness Function Evolution` | Complementary — could guide neuroevolution fitness shaping |
| `[DEFERRED] Upper-Level Intercept Director` | Could be a separate NN head or mode — NN architecture can accommodate |

---

## Approach: Weight Evolution with Fixed Topology

Use the existing GP evolution infrastructure (population, tournament selection, crossover,
mutation) to evolve the **weights** of a fixed-topology neural network, rather than evolving
tree structure. This is sometimes called "neuroevolution of fixed topologies" (NEFT) or
simply weight-space evolutionary search.

### Why Not NEAT/HyperNEAT

NEAT (NeuroEvolution of Augmenting Topologies) co-evolves topology and weights. This adds
significant implementation complexity (speciation, innovation numbers, compatibility
distance) for uncertain benefit. Our problem has:

- **Known input dimensionality** (sensor count is fixed)
- **Known output dimensionality** (3 control commands)
- **Existing evaluation infrastructure** we want to reuse

A fixed topology with evolved weights is simpler, faster to evaluate, and easier to deploy
on embedded hardware. If the topology is too small, we can increase it. If too large, we
prune. This is much simpler than NEAT's dynamic topology management.

### Why Not Gradient-Based Training

Standard NN training (backpropagation) requires differentiable loss functions and
deterministic input/output pairs. Our fitness function involves:

- A full 3D flight simulation (non-differentiable physics)
- Stochastic wind and path variations
- Crash/timeout conditions (discontinuous)
- Multi-path aggregation

Evolutionary optimization handles all of these naturally — it only needs a scalar fitness
value per individual, which we already compute.

---

## Network Architecture

### Input Vector (14 sensors)

All values already computed by existing sensor functions in `gp_evaluator_portable.cc`:

```
Input[0]  = GETDPHI(0)          // Roll angle to rabbit (rad, body frame)
Input[1]  = GETDTHETA(0)        // Pitch angle to rabbit (rad, body frame)
Input[2]  = GETDIST             // Distance to rabbit (meters)
Input[3]  = GETDIST_RATE        // Rate of distance change (m/s, clamped [-10,10])
Input[4]  = GETDPHI_RATE        // Rate of dPhi change (rad/s)
Input[5]  = GETDTHETA_RATE      // Rate of dTheta change (rad/s)
Input[6]  = GETROLL_RAD         // Aircraft roll angle (rad)
Input[7]  = GETPITCH_RAD        // Aircraft pitch angle (rad)
Input[8]  = GETVEL              // Speed magnitude (m/s)
Input[9]  = GETALPHA            // Angle of attack (rad)
Input[10] = GETBETA             // Sideslip angle (rad)
Input[11] = GETROLL             // Current roll command (-1 to 1)
Input[12] = GETPITCH            // Current pitch command (-1 to 1)
Input[13] = GETTHROTTLE         // Current throttle command (-1 to 1)
```

**Note on gravity**: Gravity direction relative to the aircraft body frame is fully
encoded by `GETROLL_RAD` + `GETPITCH_RAD`. When the aircraft is banked 30 degrees,
`GETROLL_RAD = 0.52` tells the NN that gravity has a lateral component in body frame.
When pitched up, `GETPITCH_RAD > 0` tells it gravity pulls "behind and below." No
explicit gravity sensor is needed.

**Note on cross-track**: `GETDPHI` and `GETDTHETA` are the angular bearing to the
rabbit in body frame — exactly what a forward-facing camera would see. Combined with
`GETDIST`, the NN can compute cross-track distance as `sin(dPhi) * dist` via a single
hidden layer neuron. This is the "visual sensor" equivalent discussed in design review.

### Candidate Topologies

**Minimal (14-8-3)**: 14 inputs → 8 hidden → 3 outputs
- Parameters: 14×8 + 8 + 8×3 + 3 = 139 weights
- RAM: 556 bytes (float32)
- Compute: ~139 MACs + 11 activations

**Baseline (14-16-8-3)**: 14 inputs → 16 hidden → 8 hidden → 3 outputs
- Parameters: 14×16 + 16 + 16×8 + 8 + 8×3 + 3 = 403 weights
- RAM: 1,612 bytes (float32)
- Compute: ~403 MACs + 27 activations

**Large (14-32-16-3)**: 14 inputs → 32 hidden → 16 hidden → 3 outputs
- Parameters: 14×32 + 32 + 32×16 + 16 + 16×3 + 3 = 1,011 weights
- RAM: 4,044 bytes (float32)
- Compute: ~1,011 MACs + 51 activations

### Output Layer

```
Output[0] = tanh(...)  →  pitchCommand    // [-1, 1]
Output[1] = tanh(...)  →  rollCommand     // [-1, 1]
Output[2] = tanh(...)  →  throttleCommand // [-1, 1]
```

The tanh activation naturally maps to the [-1, 1] control command range, matching the
existing GP SETPITCH/SETROLL/SETTHROTTLE clamping.

### Activation Functions

- **Hidden layers**: tanh (matches existing sensor value ranges, which are roughly centered
  around zero). Alternative: ReLU is faster to compute but loses negative gradients.
- **Output layer**: tanh (maps directly to [-1, 1] control commands).

### Recurrent Option (Future)

For temporal state without explicit _PREV/_RATE inputs, add a recurrent connection:

```
hidden_t = tanh(W_input × input + W_recurrent × hidden_{t-1} + bias)
```

This replaces the 6 temporal inputs (GETDIST_RATE, GETDPHI_RATE, GETDTHETA_RATE and
_PREV variants) with learned temporal dynamics. Increases parameter count by
`hidden_size²` but eliminates the need for ring buffer history on embedded.

**Recommendation**: Start with feedforward + explicit rate/prev inputs. Add recurrence
only if the feedforward plateau is reached and temporal modeling is identified as the
bottleneck.

---

## Embedded Platform Analysis: Seeed XIAO nRF52840

### Hardware Specs

| Aspect | Value |
|--------|-------|
| MCU | Nordic nRF52840 (ARM Cortex-M4F) |
| Clock | 64 MHz |
| FPU | VFPv4-SP-D16 (single-precision, 16 double-word registers) |
| RAM | 256 KB total, ~18 KB used by current GP system |
| Flash | 1 MB program, 2 MB external QSPI (logging) |
| Control loop | 100 ms (10 Hz) = 6,400,000 cycles budget |
| Current GP eval | ~1,000 cycles (<1% budget) |

### FPU Operation Timing (Cortex-M4F VFPv4)

| Operation | Cycles | Notes |
|-----------|--------|-------|
| VADD.F32 (add) | 1 | Single-cycle throughput |
| VMUL.F32 (multiply) | 1 | Single-cycle throughput |
| VMLA.F32 (multiply-accumulate) | 1 | Fused MAC, key for NN |
| VDIV.F32 (divide) | 14 | Avoid in hot path |
| VSQRT.F32 | 14 | Avoid in hot path |
| VCMP.F32 (compare) | 1 | For activation functions |
| Register file | 32 × 32-bit | S0-S31 (or D0-D15 as pairs) |

**Key**: VMLA.F32 (multiply-accumulate) is the workhorse for NN inference. One MAC per
cycle means a single neuron with N inputs costs N cycles for the weighted sum.

### FP16 Considerations

The Cortex-M4F **does not have native FP16 compute**. Options:

1. **FP32 compute, FP32 storage** (baseline): Simplest, uses native FPU. 4 bytes/weight.
2. **FP32 compute, FP16 storage**: Store weights in half-precision, convert to FP32 at
   load time or per-inference. Halves weight storage. ARM `__fp16` type supports this
   with `VCVT.F32.F16` instruction (1 cycle). Good for flash-constrained scenarios.
3. **INT8 quantized weights**: Store as int8, multiply by scale factor. 1 byte/weight
   but needs int→float conversion per MAC. Only worthwhile if RAM is severely constrained.
4. **CMSIS-NN library**: ARM's optimized NN inference library for Cortex-M. Supports
   int8/int16 quantized inference with SIMD (though M4 SIMD is integer-only, not FP).

**CMSIS-NN INT8 SIMD alternative**: ARM's CMSIS-NN library deliberately avoids the FPU.
Instead it uses the Cortex-M4's DSP SIMD instructions (SXTB16, SMLAD) to pack two int16
MACs into a single cycle — **2x throughput vs FP32**. The original paper reports 4.6x
overall speedup (including memory access benefits from smaller data). Accuracy loss from
int8 quantization averages ~2% on classification tasks. However, for a 1,600-parameter
network the total compute is so small that this optimization is premature.

**Recommendation**: FP32 throughout. With 256 KB RAM and ~4 KB for the largest candidate
network, memory is not the constraint. FP16 storage is a future optimization if we need
to fit multiple networks or much larger topologies. INT8 quantization via CMSIS-NN is
available if we ever scale to larger networks where the 2x throughput matters.

### Inference Timing Estimates

| Topology | MACs | Activations | Est. Cycles | Time @ 64 MHz | Budget % |
|----------|------|-------------|-------------|---------------|----------|
| 14-8-3 | 139 | 11 | ~360 | 5.6 µs | 0.006% |
| 14-16-8-3 | 403 | 27 | ~940 | 14.7 µs | 0.015% |
| 14-32-16-3 | 1,011 | 51 | ~2,030 | 31.7 µs | 0.032% |
| 14-64-32-16-3 | 3,504 | 115 | ~5,800 | 90.6 µs | 0.091% |
| 14-128-64-32-3 | 12,768 | 227 | ~17,300 | 270 µs | 0.27% |

**Activation cost estimate**: tanh approximation via polynomial (Padé or lookup) ≈ 20
cycles each. Even a 5th-order polynomial approximation is <30 cycles.

**Conclusion**: We have **enormous headroom**. Even a 14-128-64-32-3 network (12,768
parameters, ~51 KB) runs in 270 µs — 0.27% of the 100ms budget. The constraint is not
compute but whether evolution can search the weight space effectively. Smaller networks
(14-16-8-3, 403 params) are better starting points for evolutionary search.

### Sensor Evaluation Cost

The NN inputs require calling the same sensor functions the GP currently uses. These costs
are **identical** to current GP evaluation — the sensors don't get more expensive.

```
Sensor eval: ~14 calls × ~50 cycles each ≈ 700 cycles (11 µs)
NN inference: ~1,000 cycles (16 µs) for baseline topology
Total: ~1,700 cycles (27 µs) per control step
```

This is comparable to or less than the current 91-instruction bytecode evaluation.

---

## Prior Art: Neuroflight

**Neuroflight** (Koch et al., 2019) is the closest real-world analogue to this project.
It is the first open-source neural network flight controller firmware, running on real
FPV racing drones performing aerobatic maneuvers.

| Aspect | Neuroflight | Our Proposal |
|--------|-------------|--------------|
| Network | 6→32→32→4 (tanh) | 14→16→8→3 (tanh) |
| Parameters | 1,412 | 403 (baseline) |
| MCU | Cortex-M7 @ 216 MHz | Cortex-M4F @ 64 MHz |
| Inference rate | 2.67 kHz (<0.375 ms) | ~500-800 Hz est. (1-2 ms) |
| Training | RL (PPO) in Gazebo | Evolution in minisim |
| Deployment | tfcompile AOT → C | nn2cpp code gen → C |
| Control outputs | 4 motor commands | pitch, roll, throttle |

**Key validation**: A network almost identical in size to ours has been proven in real
drone flight. Their 2x32 hidden architecture with tanh activation is the same pattern
we propose. On Cortex-M4 at 64 MHz (3.4x slower clock), our smaller network (403 vs
1,412 params) would comfortably run at 500+ Hz — far exceeding our 10 Hz control loop.

---

## Embedded Deployment: Framework vs Hand-Written

### Framework Options

| Framework | Flash Overhead | RAM Overhead | INT8 SIMD | PlatformIO | Notes |
|-----------|---------------|-------------|-----------|------------|-------|
| Hand-written C | ~0 | ~0 | No | Yes | ~50 lines for full forward pass |
| NNoM | ~10 KB | ~2 KB | Optional (CMSIS-NN backend) | Yes | Pure C, Keras import |
| nn4mc | ~5 KB | ~1 KB | No | Yes | Generates standalone C from Keras |
| TFLM | 10-26 KB | ~5 KB | Yes (CMSIS-NN kernels) | Partial | Heavy for tiny nets |
| Edge Impulse | ~15 KB | ~2 KB | Yes (EON compiler) | Yes | 1ms inference on nRF52840 proven |
| CMSIS-NN direct | ~5 KB | ~1 KB | Yes | Manual | Low-level, max performance |

### Recommendation: Hand-Written Forward Pass

For a network with 403-1,635 parameters, a hand-written forward pass is the best fit:

1. **Zero dependencies** — no framework overhead, no version management
2. **~50 lines of C** — three nested loops (one per layer) + tanh LUT
3. **Matches existing pattern** — `bytecode2cpp` already generates standalone C evaluators;
   `nn2cpp` would follow the same pattern
4. **Existing tanh LUT** — `gp_evaluator_portable.cc` already has a 512-entry sin LUT
   with linear interpolation. A tanh LUT follows the same pattern.
5. **Full control** — no quantization assumptions, no allocator requirements
6. **Proven at this scale** — Edge Impulse reports 1 ms inference for similar-sized
   networks on nRF52840; hand-written will be faster.

If we later scale to larger networks (10K+ params) or need RNN/GRU layers, NNoM provides
a clean upgrade path with optional CMSIS-NN SIMD acceleration.

---

## Evolution Strategy

### Genome Representation

Each individual in the population is a **flat float array** of NN weights:

```cpp
struct NNGenome {
    std::vector<float> weights;  // All weights + biases, flattened
    int topology[MAX_LAYERS];    // Layer sizes (fixed for population)
    int num_layers;
};
```

For the baseline 14-16-8-3 topology: genome = 403 floats.

### Evolution Algorithm Options

Research identified three practical approaches for evolving 403-1,635 weights:

**Option A: Simple GA with Gaussian Mutation (recommended for Phase 1)**

Uber's Deep Neuroevolution (2017) showed that a plain truncation-selection GA with
Gaussian mutation can train networks with 4M+ parameters on Atari, competitive with
DQN and A3C. For 403 parameters, this is entirely practical and dead simple.

Pros: Simple to implement, reuses existing GP infrastructure, proven at scale.
Cons: Less sample-efficient than CMA-ES.

**Option B: Separable CMA-ES (recommended if Phase 1 plateaus)**

CMA-ES is the gold standard for continuous optimization. Full CMA-ES is O(n^2) memory
/ O(n^3) per generation, impractical at 403 dimensions. But **sep-CMA-ES** restricts
the covariance to diagonal form, scaling to 1000+ dimensions. Default population:
`lambda = 4 + floor(3 * ln(N))` = ~22 for N=403. IPOP-CMA-ES doubles population on
each restart for principled multimodality handling.

Pros: Most sample-efficient (critical with 3-min evals), self-adapting step size.
Cons: More complex to implement, may need custom integration with existing infra.

**Option C: OpenAI ES / Natural Evolution Strategies**

Designed for massive parallelism (720+ cores). With limited parallelism (8 threads),
the gradient estimate is too noisy for 403 dimensions. Not recommended at our scale.

### Population Sizing with 3-Minute Fitness Evals

| Method | Population | Serial Time/Gen | 8-Thread Time/Gen |
|--------|-----------|----------------|-------------------|
| Simple GA | 200 | 10 hours | ~1.3 hours |
| Simple GA | 500 | 25 hours | ~3.1 hours |
| sep-CMA-ES | 22 | 66 min | ~8 min |
| sep-CMA-ES | 50 | 2.5 hours | ~19 min |

CMA-ES's tiny population is a major advantage when evaluations are expensive.
Note: current GP runs use pop=500 with 8 threads at ~160s/gen — the fitness eval
is per-individual, so NN evolution has identical per-eval cost.

### Genetic Operators (Simple GA)

**Crossover**: Arithmetic crossover (blend).
```
child[i] = alpha * parent1[i] + (1-alpha) * parent2[i]
```
where alpha ~ uniform(0, 1). Preserves co-adapted weight groups better than uniform
crossover on the weight vector.

**Mutation**: Gaussian perturbation.
```
weight[i] += N(0, sigma)
```
Starting sigma = 0.1. Consider self-adaptive sigma (each individual carries its own
sigma that also mutates) for automatic step-size control.

**Selection**: Tournament selection (reuse existing GPc++ tournament infrastructure).

**Population**: 200-500 for GA; 22-50 for CMA-ES.

### Initialization

Xavier/Glorot initialization: `weight[i] = uniform(-1/sqrt(fan_in), +1/sqrt(fan_in))`.
This prevents saturation of tanh activations at initialization — critical for
neuroevolution where there is no gradient signal to recover from saturation.

---

## Integration Points

### Where the NN Evaluator Plugs Into Existing Code

The current evaluation chain has three clear integration boundaries:

#### 1. Training Mode: Replace MyGene::evaluate()

**File**: `autoc-eval.cc` line 184

Currently: GP tree recursion calling `evaluateGPOperator()` per node.

NN replacement: Call all 14 sensor functions once, run forward pass, set control commands.

```cpp
// Conceptual — replaces recursive tree evaluation
gp_scalar NNGene::evaluate(std::vector<Path>& path, MyGP& run, gp_scalar arg) {
    VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());

    // Gather inputs (reuse existing sensor functions)
    float inputs[14] = {
        executeGetDPhi(pathProvider, aircraftState, 0.0f),
        executeGetDTheta(pathProvider, aircraftState, 0.0f),
        executeGetDist(pathProvider, aircraftState),
        executeGetDistRate(aircraftState),
        executeGetDPhiRate(aircraftState),
        executeGetDThetaRate(aircraftState),
        aircraftState.getRollRad(),
        aircraftState.getPitchRad(),
        aircraftState.getRelVel(),
        executeGetAlpha(aircraftState),
        executeGetBeta(aircraftState),
        aircraftState.getRollCommand(),
        aircraftState.getPitchCommand(),
        aircraftState.getThrottleCommand(),
    };

    // Forward pass
    float outputs[3];
    nn_forward(inputs, 14, weights_, topology_, num_layers_, outputs);

    // Set control commands
    aircraftState.setPitchCommand(outputs[0]);
    aircraftState.setRollCommand(outputs[1]);
    aircraftState.setThrottleCommand(outputs[2]);

    return 0.0f;  // NN doesn't return a "tree value"
}
```

#### 2. Deployment Mode: Replace GPBytecodeInterpreter

**File**: `gp_bytecode.h` line 92

Create `NNInterpreter` class with same interface contract:

```cpp
class NNInterpreter {
    std::vector<float> weights_;
    int topology_[MAX_LAYERS];
    int num_layers_;
public:
    bool loadModel(const std::string& filename);
    gp_scalar evaluate(AircraftState& state, std::vector<Path>& path, gp_scalar arg);
};
```

Minisim detection logic (line 128-150) adds a third format check for NN weight files.

#### 3. Embedded Deployment: Replace generatedGPProgram()

**File**: `xiao-gp/generated/gp_program_generated.cpp`

Create `nn2cpp` tool (sibling to `bytecode2cpp`) that generates:

```cpp
// Auto-generated NN inference — topology 14-16-8-3, 403 weights
static const float weights[] = { /* ... 403 floats ... */ };

gp_scalar generatedNNProgram(PathProvider& pathProvider,
                             AircraftState& aircraftState,
                             gp_scalar arg) {
    float inputs[14] = { /* sensor calls */ };
    float hidden1[16], hidden2[8], outputs[3];

    // Layer 1: 14 → 16
    for (int j = 0; j < 16; j++) {
        float sum = weights[bias_offset_1 + j];
        for (int i = 0; i < 14; i++)
            sum += inputs[i] * weights[w_offset_1 + j*14 + i];
        hidden1[j] = fast_tanh(sum);
    }

    // Layer 2: 16 → 8
    // ... same pattern ...

    // Output: 8 → 3
    // ... same pattern ...

    aircraftState.setPitchCommand(outputs[0]);
    aircraftState.setRollCommand(outputs[1]);
    aircraftState.setThrottleCommand(outputs[2]);
    return 0.0f;
}
```

This compiles to ~2 KB of ARM code + ~1.6 KB weights. The existing `controller.cpp`
dispatch calls it identically to the current `generatedGPProgram()`.

#### 4. RPC Serialization

The existing `EvalData.gp` field (a `std::vector<char>` blob) can carry NN weights
using the same Boost binary serialization. Add a format magic number to distinguish
GP tree / bytecode / NN weight data.

### Fitness Computation: No Changes

The NN controller sets `pitchCommand`, `rollCommand`, `throttleCommand` on
`AircraftState` via the same setter functions the GP uses. The fitness calculation
in `autoc.cc` lines 1278-1290 reads aircraft position from the simulation trajectory,
which is identical regardless of what produced the control commands. No fitness code
changes needed.

---

## Implementation Phases

### Phase 1: Unify Evaluation Pipeline (absorbs 007-unify-eval)

This is prerequisite infrastructure. The current codebase has ~350 lines of duplicated
evaluation logic between `MyGP::evalTask()` (GP tree mode) and
`BytecodeEvaluationGP::evalTask()` (bytecode mode). Adding a third controller type (NN)
without unification would create a third copy. Instead, unify first.

- Extract `fitness_computer` — shared fitness calculation from both eval paths
- Extract `eval_data_builder` — shared scenario metadata construction
- Extract `eval_logger` — shared logging/diagnostics (data.dat, data.stc output)
- Unify core eval loop — single `evalTask()` parameterized by an evaluation backend
  interface (GP tree, bytecode, or NN)
- Define `ControllerBackend` interface: `evaluate(AircraftState&, PathProvider&) → void`
  (sets control commands on AircraftState as side effect)
- Verify: both GP tree and bytecode modes produce identical results after refactor

### Phase 2: NN Evaluator Core

- Implement `nn_forward()` — portable forward pass in `nn_evaluator_portable.cc/h`
  (same file works desktop + embedded, same pattern as `gp_evaluator_portable`)
- Implement `NNControllerBackend` conforming to the Phase 1 backend interface
- Sensor gathering: call existing executeGetDPhi, executeGetDist, etc. to build
  input vector, then run forward pass, then call setPitch/setRoll/setThrottle
- Add tanh LUT to portable math (alongside existing sin/cos/atan LUTs)
- Unit tests: forward pass correctness, output range [-1,1], known-weight regression

### Phase 3: NN Serialization & Archive

This is an incompatible change to the archive format — a fork in the road. No backward
compatibility with GP tree or bytecode archives is required.

- Define `NNGenome` serialization format: magic number + topology + flat weight array
- Implement Boost binary serialization for NNGenome (for RPC to minisim)
- Add NN format detection in minisim alongside GP tree and bytecode detection
- Define S3 archive format for NN populations: topology metadata + best weights per gen
- Update `data.dat` output: replace GP tree S-expression dump with weight vector dump
  (or summary statistics: weight mean/stdev per layer, activation distributions)
- Update `data.stc` output: same fitness/distance/attitude columns, add NN-specific
  metrics (weight magnitude, gradient proxy = fitness delta / weight delta)
- Implement `nnextractor` tool (sibling to `gpextractor`): extract best weights from
  S3 archive → standalone weight file for deployment

### Phase 4: Evolution Integration

- Implement NNPopulation: flat float arrays as genomes, tournament selection
- Implement arithmetic crossover (BLX-alpha) on weight vectors
- Implement Gaussian mutation with configurable sigma
- Xavier/Glorot initialization for initial population
- Add `ControllerType = GP | NN` config option in autoc.ini
- Add `NNTopology = 14,16,8,3` config option (comma-separated layer sizes)
- Add `NNMutationSigma = 0.1` config option
- Wire into existing elite store (re-evaluation, divergence tracking)
- First training runs: compare NN vs GP convergence on same scenario set

### Phase 5: Embedded Code Generation & Deployment

- Implement `nn2cpp` tool (sibling to `bytecode2cpp`): generates standalone C++
  inference code from weight file
  - Output: `nn_program_generated.cpp` with `generatedNNProgram()` function
  - Same function signature as `generatedGPProgram()` — drop-in replacement
  - Embeds weights as `static const float[]` in flash
  - Embeds topology, layer loop unrolling, tanh LUT call
  - Generates stack allocation sized to max layer width
- Update xiao-gp build: swap `gp_program_generated.cpp` → `nn_program_generated.cpp`
  via PlatformIO build flag or file rename
- Add `nn_evaluator_portable.cc/h` to xiao-gp includes (same shared-source pattern
  as `gp_evaluator_portable`)
- Validate timing on real nRF52840 hardware (target: <1ms inference)
- Validate control output matches desktop inference (bit-exact with same LUT)
- Flight test: compare NN controller vs GP bytecode controller on same paths

### Phase 6: Architecture Search (Optional)
- Try multiple topologies (14-8-3 vs 14-16-8-3 vs 14-32-16-3)
- Experiment with recurrent connections (simple Elman or GRU)
- Evaluate whether _PREV/_RATE inputs can be dropped with recurrence
- If GPU-native (011) is available: scale to pop=10K+ and larger topologies

---

## Functional Requirements

### FR-001: Forward Pass Evaluator
The system shall implement a portable feedforward neural network forward pass that
operates on float32 values and supports configurable layer topology.

### FR-002: Weight Evolution
The system shall evolve NN weights using the existing population management, tournament
selection, and fitness evaluation infrastructure.

### FR-003: Genetic Operators
The system shall implement arithmetic crossover (BLX-alpha) and Gaussian mutation
operators on flat weight vectors.

### FR-004: Sensor Integration
The NN evaluator shall call the same sensor functions (executeGetDPhi, executeGetDist,
etc.) as the GP evaluator, ensuring identical perception.

### FR-005: Control Output
The NN shall output three control commands (pitch, roll, throttle) in [-1, 1] via
tanh activation, set through the existing AircraftState setter functions.

### FR-006: Serialization
NN weight vectors shall be serializable via Boost binary archive for RPC transport to
minisim. The format uses a magic number prefix to distinguish from GP tree and bytecode
data. This is an incompatible format change — no backward compatibility with GP/bytecode
archives is required.

### FR-007: Format Detection
Minisim shall auto-detect NN weight format alongside GP tree and bytecode formats,
using a magic number prefix.

### FR-008: Code Generation (nn2cpp)
A `nn2cpp` tool shall generate standalone C++ inference code from trained weights,
suitable for embedded deployment (no dynamic allocation, fixed topology). Output follows
the same pattern as `bytecode2cpp`: a `generatedNNProgram()` function with identical
signature to `generatedGPProgram()`.

### FR-008a: Weight Extraction (nnextractor)
An `nnextractor` tool shall extract the best NN weights from an S3 archive and write
a standalone weight file for deployment, analogous to `gpextractor` for bytecode.

### FR-008b: Unified Eval Pipeline
The evaluation loop shall be refactored into a shared core with a pluggable controller
backend interface. GP tree, bytecode, and NN evaluation shall all use the same fitness
computation, scenario handling, and logging paths. (Absorbs 007-unify-eval.)

### FR-009: Configuration
`autoc.ini` shall support a `ControllerType` option selecting between GP tree evolution
and NN weight evolution.

### FR-010: Topology Configuration
The NN topology (layer sizes) shall be configurable in `autoc.ini`, not hardcoded.

### FR-011: Activation Function Portability
Activation functions shall use the existing LUT-based fast math (or new tanh LUT)
for platform portability between desktop and embedded.

### FR-012: Fitness Compatibility
NN-evolved controllers shall be evaluated using the identical fitness function as GP
controllers, with no changes to the fitness computation.

### FR-013: Embedded Timing
NN inference on the XIAO nRF52840 (64 MHz Cortex-M4F) shall complete within 1 ms for
the baseline topology (14-16-8-3), using less than 1% of the 100 ms control budget.

### FR-014: Weight Initialization
Initial populations shall use Xavier/Glorot initialization to prevent activation
saturation.

---

## Non-Functional Requirements

### NFR-001: GP Tree Preservation
GP tree evolution shall remain functional for ongoing flight test data generation.
However, NN serialization, archive format, and diagnostic output are intentionally
incompatible with GP formats — this is a fork in the road, not a backward-compatible
extension.

### NFR-002: Code Reuse
Sensor evaluation functions, fitness computation, RPC transport, and simulation physics
shall be shared between GP and NN evaluation paths with no duplication.

### NFR-003: Minimal New Dependencies
The NN forward pass shall be implemented in plain C++ with no external ML framework
dependencies (no TensorFlow, PyTorch, ONNX Runtime). The entire inference is a few
nested loops.

---

## GPU-Native Evaluation Synergy (011-gpu-native)

Neuroevolution becomes dramatically more practical under GPU-native simulation. This
section analyzes how 013 and 011 interact.

### Current CPU Bottleneck

Each fitness evaluation runs 36 scenarios × ~100 timesteps × physics sim = ~160s per
individual on 8 CPU threads. This makes population sizing the binding constraint for
neuroevolution — large populations are prohibitively slow.

### GPU Changes the Economics

The GB10 can batch 64-256 physics simulations simultaneously. If the physics sim
moves to CUDA (011), the time/generation drops by 10-100x:

| Platform | Evals/sec (est.) | Pop=500 Time/Gen | Pop=10K Time/Gen |
|----------|-----------------|------------------|-------------------|
| CPU 8-thread (current) | ~50 | ~160s | ~53 min |
| GPU 64-batch | ~1,000-5,000 | ~4-18s | ~1-6 min |
| GPU 256-batch | ~5,000-20,000 | ~1-4s | ~15-60s |

At GPU speed, pop=10,000 takes about the same time as current pop=500. This unlocks:

- **Large populations**: 5K-20K individuals with Gaussian mutation, matching the regime
  where Uber's Deep Neuroevolution was proven (they used 1K for 4M params; we'd use
  10K for 403 params — massively over-sampled, very reliable convergence).
- **CMA-ES becomes less attractive**: Its advantage is sample efficiency with expensive
  evals. When evals are cheap, brute-force population search wins.
- **Topology search**: Run multiple sub-populations with different architectures
  (14-8-3, 14-16-8-3, 14-32-16-3) in parallel and compare within a single run.
- **Full scenario coverage**: With 36+ wind/entry/speed scenarios per eval and fast
  sim, every individual sees the full variation space every generation.

### NN Eval Is Naturally GPU-Friendly

A GP tree has variable structure and recursive evaluation — hard to batch efficiently
on GPU (each tree follows a different execution path). An NN forward pass is a
**fixed-size matrix multiply** — literally what GPUs are built for.

Per-generation GPU workload for NN eval:
```
Pop=500 × 36 scenarios × 100 timesteps = 1.8M forward passes per generation
Each forward pass: 14→16→8→3 = ~400 MACs

Total: 1.8M × 400 = 720M MACs per generation
GPU throughput: ~10 TFLOPS (GB10 FP32)
Time: 720M / 10T = 0.072 ms (effectively free)
```

The NN eval overhead is negligible. Physics simulation dominates entirely. This means
the controller type (GP vs NN) doesn't affect GPU sim throughput at all — the
bottleneck is always the physics integration, aerodynamics, and state management.

### GPU Integration Architecture

The NN weights for each individual are a small constant buffer uploaded once per
individual. During the physics sim loop on GPU, each timestep:

1. Gather sensor inputs from GPU-side aircraft state (already in GPU memory)
2. Run NN forward pass (batched GEMM across all concurrent sims)
3. Write control outputs to GPU-side aircraft state
4. Advance physics (the expensive part)

Steps 1-3 are a single batched matrix multiply. Step 4 is the existing GPU physics
kernel from 011. No CPU↔GPU round-trips during simulation — weights upload once,
results download once.

For GP tree eval on GPU, the bytecode interpreter could also work (uniform instruction
dispatch), but would be less efficient than NN due to variable-length programs and
branch-heavy execution. NN's fixed compute graph is inherently better suited.

### Dependency Analysis

| Spec | Disposition |
|------|-----------|
| 007 (unify eval) | **Absorbed into 013 Phase 1** — eval pipeline unification is the first phase |
| 004 (vacate boost) | Independent, but simplifies serialization for GPU path later |
| 011 (GPU native) | Future — depends on 013 Phase 1 (unified eval backend interface) |

**Sequencing**: 013 Phase 1-5 (CPU) → 004 (vacate boost) → 011 (GPU) → scale pop to 10K+

013 Phase 1 (unify eval) creates the clean backend interface that 011 needs anyway.
When 011 lands, scaling up population is a config change (`PopulationSize=10000`),
not a code change.

---

## Open Questions

1. **Population size**: Weight-space evolution on 403 dimensions may need larger
   populations than tree GP. Start with existing 500-20K range and measure.

2. **Mutation rate / sigma**: Self-adaptive sigma (CMA-ES style) vs fixed sigma?
   Start simple with fixed sigma=0.1, tune empirically.

3. **Recurrence**: Should Phase 1 include simple recurrence (Elman network), or
   start feedforward-only? Feedforward is simpler and may be sufficient given
   the explicit rate inputs.

4. **Input normalization**: Sensor values have different ranges (dPhi in [-pi, pi],
   dist in [0, 100], vel in [0, 30]). Normalize inputs to [-1, 1]? Or let evolution
   discover appropriate scaling in the first layer weights?

5. **Elitism**: Current GP uses elite store with re-evaluation. Same strategy for NN,
   or different? NN fitness should be more stable (no structural changes between
   generations) so elite divergence may be less of an issue.

6. **Hybrid approach**: Could GP trees evolve alongside NN individuals in the same
   population? Probably not — the genome representations are incompatible. But
   could run GP and NN populations in parallel and compare.

---

## Research References

### Neuroevolution for Flight Control
- **Neuroflight** (Koch et al., 2019): First open-source NN flight controller firmware.
  6→32→32→4 tanh, 1412 params, 2.67 kHz on Cortex-M7. Trained via PPO in Gazebo.
  [arXiv:1901.06553](https://arxiv.org/abs/1901.06553)
- **Neuroevolutionary soaring** (Aerospace 2021): NEAT-evolved controllers for sustained
  dynamic and thermal soaring. Domain randomization for robustness.
- **Quadcopter neuroevolution** (MDPI Aerospace 2023): NEAT-based neurocontrollers
  producing interpretable, compact networks.

### Evolution Strategies
- **Deep Neuroevolution** (Uber, 2017): Plain GA with Gaussian mutation trains 4M+
  parameter networks on Atari, competitive with DQN/A3C.
  [arXiv:1712.06567](https://arxiv.org/abs/1712.06567)
- **CMA-ES Tutorial** (Hansen): Gold standard for continuous optimization.
  sep-CMA-ES scales to 1000+ dimensions with diagonal covariance.
  [arXiv:1604.00772](https://arxiv.org/abs/1604.00772)
- **OpenAI ES** (Salimans et al., 2017): Scalable alternative to RL, designed for
  massive parallelism (720+ cores). Less practical at 8 threads.

### Embedded NN Inference
- **CMSIS-NN** (ARM, 2018): 4.6x speedup via int8/int16 SIMD on Cortex-M. ~2% accuracy
  loss. [arXiv:1801.06601](https://arxiv.org/abs/1801.06601)
- **TFLM** (Google, 2021): <2KB interpreter core, zero dynamic allocation.
  Tested on nRF52840 DK with Nordic support.
- **Edge Impulse on nRF52840**: 1 ms inference, 1.5KB RAM for motion recognition.
- **nn4mc** (Correll Lab, CU Boulder): Keras→C code generation for MCUs.
  [github.com/correlllab/nn4mc](https://github.com/correlllab/nn4mc_cpp)
- **NNoM**: Pure C NN library for MCUs, optional CMSIS-NN backend, Keras import.
  [github.com/majianjia/nnom](https://github.com/majianjia/nnom)

### GP vs NN for Control
- NN weight spaces are inherently smooth — small weight changes produce small behavior
  changes. GP tree spaces are rugged — single crossover can drastically change behavior.
  This smoothness makes NNs more amenable to evolution strategies.
- For continuous control, evolved NNs and GP trees show comparable deployed performance,
  but NNs scale better to higher input dimensionality.

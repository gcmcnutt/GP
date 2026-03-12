# Data Model: 013-neuroevolution

## Entities

### NNGenome

The fundamental unit of evolution. A flat weight array representing a fixed-topology
feedforward neural network.

| Field | Type | Description |
|-------|------|-------------|
| weights | `std::vector<float>` | All weights + biases, flattened in layer-major order |
| topology | `std::vector<int>` | Layer sizes, e.g., {14, 16, 8, 3} |
| fitness | `double` | Aggregated fitness from evaluation |
| generation | `uint32_t` | Generation when this genome was created |
| mutation_sigma | `float` | Per-individual mutation step size (self-adaptive) |

**Invariants**:
- `weights.size()` must equal the sum of `topology[i] * topology[i+1] + topology[i+1]`
  for all layers (weights + biases)
- `topology[0]` = number of sensor inputs (14 for current config)
- `topology[last]` = 3 (pitch, roll, throttle outputs)
- All weights are finite (no NaN/Inf)

**Weight layout** (layer-major, weights then biases):
```
Layer 0→1: weights[0 .. input×hidden1-1], biases[.. hidden1-1]
Layer 1→2: weights[.. hidden1×hidden2-1], biases[.. hidden2-1]
Layer 2→3: weights[.. hidden2×output-1],  biases[.. output-1]
```

### NNPopulation

Collection of NNGenome individuals with shared topology and evolution parameters.

| Field | Type | Description |
|-------|------|-------------|
| individuals | `std::vector<NNGenome>` | Population of genomes |
| topology | `std::vector<int>` | Shared topology (all individuals same shape) |
| population_size | `int` | Number of individuals |
| generation | `uint32_t` | Current generation counter |
| best_fitness | `double` | Best fitness in current generation |
| best_index | `int` | Index of best individual |

### ControllerBackend (Interface)

Abstract interface for controller evaluation. All controller types implement this.

| Method | Signature | Description |
|--------|-----------|-------------|
| evaluate | `void evaluate(AircraftState&, PathProvider&)` | Run controller, set control commands on AircraftState |
| getName | `const char* getName()` | Return controller type name for logging |

**Implementations**:
- `GPTreeBackend` — wraps existing `MyGene::evaluate()` recursive tree traversal
- `BytecodeBackend` — wraps existing `GPBytecodeInterpreter::evaluate()`
- `NNControllerBackend` — gathers sensor inputs, runs `nn_forward()`, sets commands

### NNSerializationFormat

Binary format for RPC transport and S3 storage.

```
Offset  Size    Field
0       4       Magic bytes: "NN01"
4       4       Format version (uint32, currently 1)
8       4       num_layers (uint32)
12      4×N     layer_sizes[num_layers] (uint32 each)
12+4N   4       num_weights (uint32)
16+4N   4×W     weights[num_weights] (float32 each)
16+4N+4W 8      fitness (float64)
24+4N+4W 4      generation (uint32)
28+4N+4W varies  s3_key (null-terminated string)
```

### FitnessComputer

Extracted shared fitness computation (from Phase 1 unification).

| Field | Type | Description |
|-------|------|-------------|
| distance_target | `double` | Target following distance (7.5m) |
| distance_norm | `double` | Distance penalty normalization (5.0) |
| distance_power | `double` | Distance penalty exponent (0.75) |
| attitude_norm | `double` | Attitude delta normalization (0.349) |
| attitude_power | `double` | Attitude delta exponent (1.5) |

**Methods**:
- `computeStepPenalty(distance, attitude_delta, intercept_scale) → double`
- `computeCrashPenalty(fraction_completed) → double`
- `computeAttitudeScale(path_distance, path_turn_rad) → double`

### EvalLogger

Extracted shared logging (data.dat, data.stc output).

| Method | Description |
|--------|-------------|
| logStepHeader | Write column headers to data.dat |
| logStep | Write per-timestep data row |
| logGenerationStats | Write per-generation summary to data.stc |
| logBestController | Write best controller representation (GP tree S-expr or NN weight summary) |

## State Transitions

### NNGenome Lifecycle

```
INITIALIZED (Xavier/Glorot random weights)
    ↓ evaluate across all scenarios
EVALUATED (fitness assigned)
    ↓ tournament selection
SELECTED (chosen as parent)
    ↓ arithmetic crossover + Gaussian mutation
OFFSPRING (new weights, unevaluated)
    ↓ evaluate across all scenarios
EVALUATED (fitness assigned)
    ↓ ... repeat
```

### Controller Mode Selection

```
autoc.ini: ControllerType=GP
    → GPTreeBackend created
    → GP tree population evolved
    → gpextractor → bytecode → bytecode2cpp → xiao-gp

autoc.ini: ControllerType=NN
    → NNControllerBackend created
    → NNPopulation evolved
    → nnextractor → weight file → nn2cpp → xiao-gp
```

## Relationships

```
NNPopulation 1──* NNGenome
NNGenome uses→ NNControllerBackend (wraps genome for evaluation)
NNControllerBackend implements→ ControllerBackend
GPTreeBackend implements→ ControllerBackend
BytecodeBackend implements→ ControllerBackend
FitnessComputer ←used by── evalTask() (shared across all backends)
EvalLogger ←used by── evalTask() (shared across all backends)
```

# Data Model: 014-nn-training-signal

## Entities

### NNGenome (existing, no changes)
- `topology`: vector of layer sizes [22, 16, 8, 3]
- `weights`: flat vector of 531 floats (layer-contiguous)
- `fitness`: double (lower is better)
- `generation`: int
- Serialized as "NN01" magic + binary

### CMAESState (new)
- `mean`: vector<double> (N=531) — distribution center
- `diagonal_cov`: vector<double> (N=531) — diagonal covariance
- `p_sigma`: vector<double> (N=531) — conjugate evolution path
- `p_c`: vector<double> (N=531) — covariance evolution path
- `sigma`: double — global step size
- `generation`: int
- `best_weights`: vector<double> (N=531) — best solution found
- `best_fitness`: double
- `rng_state`: uint64_t seed (for reproducibility)

### CurriculumSchedule (new)
- `stages`: vector of CurriculumStage
- `current_stage_index`: int

### CurriculumStage
- `start_generation`: int
- `scenario_count`: int (how many scenarios to evaluate)
- `description`: string (for logging)

### FitnessAggregator (new)
- `mode`: enum { SUM, MINIMAX, PERCENTILE }
- `percentile_value`: double (e.g., 0.95 — only used in PERCENTILE mode)

### TimestepRecord (new, for per-timestep streaming)
- `time_ms`: int
- `distance_to_target`: double
- `dphi`: double (roll angle to target)
- `dtheta`: double (pitch angle to target)
- `pitch_cmd`: double
- `roll_cmd`: double
- `throttle_cmd`: double

### TrajectorySegment (new, for per-segment scoring)
- `start_time_ms`: int
- `end_time_ms`: int
- `start_distance`: double
- `end_distance`: double
- `error_reduction`: double (start_distance - end_distance)
- `difficulty_weight`: double (higher for turns/crosswind)
- `segment_score`: double

### EvolutionCheckpoint (new)
- `generation`: int
- `curriculum_stage_index`: int
- `population`: vector<NNGenome> (or CMAESState if using CMA-ES)
- `optimizer_type`: enum { GA, SEP_CMA_ES }
- `cmaes_state`: CMAESState (if optimizer_type == SEP_CMA_ES)
- `rng_state`: uint64_t

### RunConfig (replaces GPConfiguration + GPVariables)
- All config from autoc.ini parsed via inih
- Key fields: population_size, num_generations, eval_threads, s3_bucket, s3_profile, nn_mutation_sigma, nn_sigma_floor, optimizer_type, curriculum_enabled, fitness_aggregation_mode

## Relationships

```
RunConfig ──configures──> EvolutionLoop
EvolutionLoop ──uses──> CMAESState OR NNPopulation (GA mode)
EvolutionLoop ──uses──> CurriculumSchedule
EvolutionLoop ──uses──> FitnessAggregator
EvolutionLoop ──produces──> EvolutionCheckpoint (per generation)
EvolutionLoop ──evaluates──> NNGenome via Simulator
Simulator ──returns──> vector<TimestepRecord> (per scenario)
FitnessAggregator ──scores──> vector<TrajectorySegment> OR scalar
```

## State Transitions

### Evolution Loop
```
INIT → EVALUATE → AGGREGATE_FITNESS → SELECT → REPRODUCE → [CHECKPOINT] → EVALUATE → ...
  │                                                              │
  └── curriculum stage check on generation boundary ─────────────┘
```

### CMA-ES Loop (replaces SELECT/REPRODUCE)
```
ASK (sample lambda candidates from distribution) →
EVALUATE (all candidates across scenarios) →
TELL (rank, update mean/sigma/covariance) →
[CHECKPOINT] → ASK → ...
```

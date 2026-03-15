# Feature Specification: NN Training Signal Improvement

**Feature Branch**: `014-nn-training-signal`
**Created**: 2026-03-14
**Status**: Clarified
**Input**: User description: "Improve NN training signal for path-following under variations. The current simple GA with scalar fitness (49 scenarios x 600 timesteps -> one number) hit a wall in nn13 (sigma collapsed 0.20->0.024, fitness stalled at 3.99M). This feature bundles related efforts to get better gradient signal from the simulator."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Prevent Search Freeze via Sigma Floor (Priority: P1)

As an operator running NN evolution, I want the mutation step size to never collapse below a useful threshold, so that evolution continues exploring the weight space instead of stalling.

**Why this priority**: This is the simplest fix for the immediate problem observed in nn13. Self-adaptive sigma decayed from 0.20 to 0.024, effectively freezing the search. A sigma floor directly prevents this failure mode with minimal code change.

**Independent Test**: Can be tested by running a short evolution (50-100 generations) and verifying that mutation sigma never drops below the configured floor, and that fitness continues improving past the point where nn13 stalled.

**Acceptance Scenarios**:

1. **Given** an NN evolution run with sigma floor configured, **When** self-adaptive sigma would decay below the floor value, **Then** the system clamps sigma to the floor value and logs the clamping event.
2. **Given** a sigma floor of 0.05, **When** evolution runs for 500 generations, **Then** sigma never falls below 0.05 at any point in the run.
3. **Given** sigma floor is set to 0, **When** evolution runs, **Then** the system behaves identically to the current unclamped behavior (backward compatible).

---

### User Story 2 - Curriculum Scenario Ramp for NN (Priority: P1)

As an operator, I want NN evolution to start with easy scenarios (few wind variations, simple paths) and progressively ramp to the full 49-scenario suite, so that the NN can learn basic flight control before being challenged with hard conditions.

**Why this priority**: The existing VariationRampStep infrastructure (from features 003/005) already supports scenario ramping for GP evolution but is not wired up for NN. Connecting it gives the NN a progressive learning curriculum, which research shows dramatically improves convergence for hard optimization problems.

**Independent Test**: Can be tested by running NN evolution and observing that early generations evaluate against fewer/easier scenarios and later generations evaluate against the full suite.

**Acceptance Scenarios**:

1. **Given** curriculum ramping is enabled with 3 stages, **When** NN evolution starts, **Then** generation 0 evaluates against only the first stage's scenario subset.
2. **Given** a curriculum schedule (e.g., stage 1: 1 scenario for gens 0-50, stage 2: 7 scenarios for gens 51-150, stage 3: 49 scenarios for gens 151+), **When** evolution crosses a stage boundary, **Then** the scenario count increases and a log message records the transition.
3. **Given** curriculum ramping is disabled, **When** NN evolution runs, **Then** all 49 scenarios are used from generation 0 (backward compatible).

---

### User Story 3 - Per-Scenario Fitness Decomposition (Priority: P2)

As an operator, I want fitness aggregation to penalize individuals that fail catastrophically on any single scenario, rather than allowing a good average to mask a bad worst case.

**Why this priority**: The current sum-over-all-scenarios aggregation lets an individual score well on 48 easy scenarios while failing on 1 hard scenario. Switching to worst-case (minimax) or 95th-percentile aggregation forces evolution to find controllers that are robust across all conditions.

**Independent Test**: Can be tested by comparing fitness rankings of a population under sum vs minimax aggregation and verifying that minimax penalizes individuals with high variance across scenarios.

**Acceptance Scenarios**:

1. **Given** minimax aggregation is selected, **When** an individual scores {100, 100, 100, 9000} across 4 scenarios, **Then** its fitness is 9000 (worst case), not 2325 (mean).
2. **Given** 95th-percentile aggregation is selected, **When** fitness is computed across 49 scenarios, **Then** the reported fitness equals the value at the 95th percentile of per-scenario scores.
3. **Given** sum aggregation is selected (default), **When** fitness is computed, **Then** behavior is identical to the current system.

---

### User Story 4 - sep-CMA-ES Optimizer (Priority: P2)

As an operator, I want to replace Gaussian mutation with sep-CMA-ES, so that evolution learns the covariance structure of the 531-weight space and converges more efficiently with a much smaller population.

**Why this priority**: sep-CMA-ES is the gold standard for continuous optimization in the 100-1000 dimension range. It replaces both self-adaptive sigma and arithmetic crossover with a principled covariance-adapted search. Population drops from 5000 to ~50-100, dramatically reducing per-generation simulation cost.

**Independent Test**: Can be tested by running sep-CMA-ES on a known benchmark function (e.g., Rosenbrock) to verify convergence, then on the actual NN flight controller problem to compare against baseline GA.

**Acceptance Scenarios**:

1. **Given** sep-CMA-ES is selected as the optimizer, **When** evolution runs, **Then** the system uses CMA-ES sampling (mean + sigma * sqrt(diagonal covariance) * normal) instead of Gaussian mutation.
2. **Given** sep-CMA-ES with population size 50, **When** evolution runs for 500 generations, **Then** the per-generation wall-clock time is significantly less than the current 5000-individual GA.
3. **Given** sep-CMA-ES evolution, **When** the operator inspects logs, **Then** step size (sigma), best fitness, and mean fitness are reported per generation.

---

### User Story 5 - Per-Segment Credit Assignment (Priority: P3)

As an operator, I want trajectory scoring to break 60-second flights into overlapping segments and score each segment by error reduction relative to its starting condition, so that the fitness signal rewards local improvements rather than only end-to-end performance.

**Why this priority**: A single scalar for 600 timesteps provides poor gradient signal. Segment-level scoring amplifies the signal by rewarding controllers that locally reduce error, even if the overall trajectory is mediocre. This requires per-timestep fitness streaming (Story 7).

**Independent Test**: Can be tested by computing segment scores on recorded trajectory data and verifying that a controller with good local corrections scores higher than one with constant error.

**Acceptance Scenarios**:

1. **Given** segment scoring is enabled with 5-second segments, **When** a 60-second trajectory is evaluated, **Then** the system produces 12 segment scores (or more if overlapping).
2. **Given** a segment where distance error decreases from 15m to 5m, **When** segment score is computed, **Then** the score reflects the error reduction (delta-based) weighted by segment difficulty.
3. **Given** segments are weighted by difficulty, **When** a segment occurs during a sharp turn with crosswind, **Then** its weight is higher than a straight-line calm segment.

---

### User Story 6 - Behavioral Cloning Bootstrap (Priority: P3)

As an operator, I want to initialize NN weights from supervised learning on recorded GP flight data, so that evolution starts from a reasonable controller rather than random Xavier weights.

**Why this priority**: The best GP controllers already achieve good path following. Training the NN via supervised MSE on GP state-action pairs provides a warm start that should be dramatically better than random initialization. Evolution then refines rather than discovers.

**Independent Test**: Can be tested by training an NN on recorded GP flight data, deploying it in the simulator, and verifying it produces reasonable (if imperfect) path following before any evolutionary refinement.

**Acceptance Scenarios**:

1. **Given** recorded GP flight data (state-action pairs at 10Hz), **When** behavioral cloning training runs, **Then** the MSE loss converges and the resulting weights produce a controller that follows paths (even if poorly).
2. **Given** a behaviorally-cloned NN, **When** it is used as the initial population seed for evolution, **Then** generation-0 fitness is significantly better than Xavier-initialized random weights.
3. **Given** no GP flight data is available, **When** behavioral cloning is requested, **Then** the system reports an error and falls back to Xavier initialization.

---

### User Story 7 - Per-Timestep Fitness Streaming (Priority: P3)

As an operator, I want the simulation to return per-timestep fitness data (distance, attitude, commands) rather than just a single aggregate scalar, so that advanced scoring methods (segment scoring, behavioral cloning data collection, policy gradient) have the raw data they need.

**Why this priority**: This is infrastructure that enables Stories 5, 6, and future policy gradient work. Without per-timestep data, advanced scoring is impossible.

**Independent Test**: Can be tested by running a single simulation evaluation and verifying that the returned data includes per-timestep distance, attitude delta, and command values.

**Acceptance Scenarios**:

1. **Given** per-timestep streaming is enabled, **When** a 60-second simulation completes, **Then** the system returns ~600 timestep records each containing distance, attitude delta, and 3 command values.
2. **Given** per-timestep streaming, **When** the aggregate fitness is computed from the streamed data, **Then** it matches the legacy single-scalar fitness value (backward compatibility).
3. **Given** per-timestep streaming is disabled, **When** simulation runs, **Then** behavior is identical to the current system (only aggregate scalar returned).

---

### User Story 8 - Checkpoint/Resume (Priority: P3)

As an operator, I want evolution state to be checkpointed at each generation so that crashed runs can be resumed and curriculum stage transitions are seamless.

**Why this priority**: Long evolution runs (500+ generations) are vulnerable to crashes, network interruptions, and power failures. Checkpointing also enables curriculum stage transitions where the optimizer state is preserved across difficulty ramps.

**Independent Test**: Can be tested by running evolution for 10 generations, killing the process, resuming, and verifying that generation 11 continues correctly with the same population state.

**Acceptance Scenarios**:

1. **Given** checkpointing is enabled, **When** generation N completes, **Then** the full evolution state (population, optimizer state, generation counter, curriculum stage) is written to a recoverable format.
2. **Given** a checkpoint exists at generation N, **When** the operator resumes evolution, **Then** generation N+1 begins with the exact population and optimizer state from the checkpoint.
3. **Given** no checkpoint exists, **When** resume is requested, **Then** the system starts fresh from generation 0.

---

### Edge Cases

- What happens when sigma floor is set higher than the initial sigma? The system should log a warning and use the floor as the effective sigma.
- What happens when curriculum ramp reaches the final stage mid-run and the operator increases total generations? The system should stay at the final stage for all remaining generations.
- What happens when per-timestep data streaming increases RPC payload size beyond available memory? The system should enforce a maximum trajectory length or use streaming/chunked transfer.
- What happens when a checkpoint file is corrupted? The system should detect corruption (via magic bytes or checksum) and fall back to starting fresh with a warning.
- What happens when sep-CMA-ES population size is set below the minimum recommended (e.g., <4+3*ln(N))? The system should warn and clamp to the recommended minimum.

## Requirements *(mandatory)*

### Functional Requirements

**Tier 0 — Prerequisite Cleanup:**

- **FR-012**: System MUST remove GP controller code (tree evaluation, bytecode interpreter, gpextractor) from the NN fork, preserving only the interface contracts required by NN evolution: RPC transport protocol, S3 archive layer, and ini-file configuration system.
- **FR-013**: System MUST replace the current multi-format Boost-serialized RPC transport (GP tree, bytecode, NN) with a single unified wire format for NN weights + per-timestep response data. Train and eval modes MUST use the same transport from the evaluator side. Boost serialization dependency MUST be removed from the transport layer.
- **FR-014**: System MUST remove all dependencies on the parent GP repo (libgp.a, GP headers, parent Makefile) so that autoc/ is fully self-contained and can be extracted into its own standalone repository as the final step.

**Tier 1 — Quick Wins:**

- **FR-001**: System MUST enforce a configurable minimum mutation sigma (sigma floor) that prevents self-adaptive sigma from decaying below a specified threshold.
- **FR-002**: System MUST support curriculum scenario ramping for NN evolution, starting with a subset of scenarios and progressively adding more according to a configurable schedule.
- **FR-003**: System MUST support configurable fitness aggregation across scenarios: sum (default), minimax (worst-case), or percentile-based (e.g., 95th percentile).

**Tier 2 — Better Evolution:**

- **FR-004**: System MUST support sep-CMA-ES as an alternative optimizer, replacing Gaussian mutation and arithmetic crossover with covariance-adapted sampling.
- **FR-005**: System MUST support per-segment credit assignment that scores trajectory segments by error reduction relative to starting conditions, weighted by segment difficulty.
- **FR-006**: *(DEFERRED to backlog)* Behavioral cloning from GP flight data. Only relevant if direct NN training is exhausted and GP→NN weight transfer is attempted.

**Tier 3 — Infrastructure:**

- **FR-007**: System MUST support per-timestep fitness data streaming from the simulator, providing distance, attitude, and command data at each timestep rather than only aggregate fitness.
- **FR-008**: System MUST support checkpointing of full evolution state (population, optimizer state, generation, curriculum stage) at configurable intervals for crash recovery and resume.

**Cross-Cutting:**

- **FR-009**: The surviving contract is sensor-in/control-out: NN training (autoc → sim evaluation → fitness), NN export (nnextractor → nn2cpp → xiao), and the NN01 serialization format. Backward compatibility with GP-era code, transport versioning, and modal branching (if GP/if NN) is explicitly NOT required — simplify aggressively to enable fast iteration.
- **FR-010**: All new features MUST be configurable via the ini-file configuration system. The config format may change (e.g., add sections) — no backward compatibility with the current flat format is required, only the same data. Inline comments should be supported.
- **FR-011**: sep-CMA-ES optimizer state MUST be included in checkpoint/resume data.

### Key Entities

- **Sigma Floor**: A minimum bound on the per-individual mutation step size, preventing search freeze. Configured per-run, applied during mutation.
- **Curriculum Schedule**: A sequence of stages, each specifying which scenarios to evaluate and at which generation to transition. Drives progressive difficulty ramping.
- **Fitness Aggregator**: A strategy for combining per-scenario fitness values into a single individual fitness. Supports sum, minimax, and percentile modes.
- **sep-CMA-ES State**: The optimizer's internal state including mean vector, step size, diagonal covariance, and evolution path. Replaces per-individual sigma and crossover.
- **Trajectory Segment**: A contiguous slice of a simulation trajectory (e.g., 5-10 seconds), scored by error reduction from its starting condition, weighted by difficulty metrics (turn rate, crosswind).
- **Behavioral Cloning Dataset**: A collection of (state, action) pairs recorded from GP controller flights at 10Hz, used for supervised MSE training of NN weights.
- **Timestep Record**: Per-timestep simulation output containing distance to target, attitude delta, and 3 command values. Foundation for segment scoring and future policy gradient.
- **Evolution Checkpoint**: A snapshot of complete evolution state at a given generation, sufficient to resume evolution exactly where it left off.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: NN evolution with sigma floor enabled MUST show continued fitness improvement beyond generation 500 (the point where nn13 stalled), achieving at least 20% lower fitness than the 3.99M plateau.
- **SC-002**: Curriculum-trained NN MUST achieve equivalent or better final fitness compared to non-curriculum training, in fewer total scenario evaluations.
- **SC-003**: Minimax or percentile fitness aggregation MUST produce controllers where the worst-case scenario fitness is at least 30% better than controllers evolved with sum aggregation.
- **SC-004**: sep-CMA-ES MUST achieve equivalent fitness to the current GA in at most 1/10th the total number of individual evaluations (e.g., 50 population x 500 gens = 25K evals vs 5000 x 500 = 2.5M evals).
- **SC-005**: Behaviorally-cloned NN weights MUST produce generation-0 fitness at least 50% better than Xavier-initialized random weights.
- **SC-006**: Checkpoint/resume MUST produce bit-identical evolution trajectories — a run interrupted and resumed at generation N MUST produce the same results as an uninterrupted run through generation N+K.

## Clarifications

### Session 2026-03-14

- Q: Is GP removal in-scope for this branch? → A: Yes, as prerequisite cleanup. Spec defines surviving contracts. Deep dependency unwinding expected.
- Q: What survives in the RPC transport after GP removal? → A: Strip GP/bytecode formats AND replace Boost serialization with simpler wire format. Single unified transport for train and eval (evaluator side is identical in both modes).
- Q: Should S3 archive layer be abstracted or NN-only? → A: NN-only, strip GP archive code. No backward compatibility required — all codebases are forking. The surviving contract is sensor-in/control-out (training loop + xiao export). Simplify aggressively: drop transport versioning, modal if/else branches, and fancy interfaces to enable fast iteration on fitness, sensors, and training styles.
- Q: Does autoc still depend on libgp.a? → A: No. Drop libgp.a entirely — NNPopulation/NNGenome are self-contained. Clean out all dependencies on the parent GP repo so autoc/ can be extracted into its own standalone repo as the final step of this feature.
- Q: What goes into the standalone repo? → A: Start with autoc/ contents only (CMake already builds independently). After extraction, re-examine and migrate relevant specs/docs/backlog items from the GP repo as follow-up work in the new repo.
- Q: What operational interfaces survive refactoring? → A: S3 archiving machinery (profiles, bucket config, key name synthesis) stays as-is. Diagnostic outputs (data.dat, data.stc) and console logging structure stay in the same general scheme. During refactoring, when encountering logging/diagnostics/archiving, reason about what stays, what gets simplified, and where it moves — don't just delete without thought.
- Q: How do we manage risk during the large refactor? → A:
  - **Contract tests first**: Add contract tests around key modules (evaluator, config, basic evolution engine) before ripping things out. These tests define the surviving behavior.
  - **Config format**: The ini parser is part of libgp and will be replaced. Any contemporary format is fine (TOML, JSON, YAML, etc). Command-line args remain similarly.
  - **Happy path only**: Strip all defensive hack-arounds, fallback guesses for missing data, if-this-or-that-isn't-there paths. Expect all-or-nothing — the complexity those branches add is itself a problem.
  - **Analysis before surgery**: Before incremental build-fix-crash cycles, do thorough upfront analysis: enumerate all definitions exposed by libgp.a (`nm`), find all transitive includes of gp.h, scan for every symbol from libgp used in autoc code. Understand the full dependency surface, then plan cuts. Less "mushing around", more informed removal.
  - **Wire format**: The versioned S3 archive format from autoc trainer is fine, keep it simple. Renderer reads from the same archives so it stays aligned.
  - **Source reorganization**: As GP code is removed, organize autoc/ into a cleaner directory structure: `util/` (logging, data output, config, S3), `eval/` (evaluator, simulation interface), `nn/` (NNGenome, NNPopulation, nn_forward, topology), with main programs at top level. The newer NN code is already better factored — extend that pattern. Can lean slightly more into modern C++ idioms now that the C-style GP library constraint is going away.
  - **Source layout**: Adopt conventional C++ project structure: `include/autoc/` for public headers (with subdirs `nn/`, `eval/`, `util/`), `src/` for implementation (mirroring include subdirs), `tests/` for all test sources, `tools/` for secondary executables (nnextractor, nn2cpp, renderer). Main executables (autoc, minisim) at `src/` top level. Headers use `#include "autoc/nn/nn_topology.h"` style for clear provenance. All file moves via `git mv` to preserve history.
  - **Test organization**: Test binaries should build into a subdirectory (e.g., `build/tests/`) rather than cluttering the top-level build dir alongside main executables. Use conventional CMake patterns (`add_subdirectory(tests)`, `ctest`). Tests always build and run — no conditional compilation or `BUILD_TESTS` flags.
  - **Run output management**: Now that we own data.dat/data.stc (no longer constrained by libgp), clean up the output story: output file paths should be configurable params, no more hardcoded `eval-` prefixes, and all run artifacts (data.dat, data.stc, logs, detailed logs, weight files) should go into an auto-created run subdirectory (via startup script or similar). This prevents test runs from overwriting production data and sets up for eventual S3 archiving of full run artifacts. The current `stdbuf` piping and primitive logging can be improved as part of this.
  - **Zero parent dependency**: When done, autoc's CMake is fully self-contained — no references to parent GP repo files, headers, or build artifacts. This also means the shared code included in crrcsim and xiao-gp becomes cleaner: the `#ifdef GP_BUILD` / `#ifdef GP_TEST` guards exist today to prevent GP infrastructure from leaking into those systems. Once GP infra is removed from autoc itself, the shared evaluator code becomes genuinely portable without conditional compilation masks. The xiao-gp Arduino engine and crrcsim's own build system each have their own concerns, but the shared NN evaluator contract (sensor-in/control-out) should be includable without bringing GP baggage along.
- Q: What is the implementation order? → A: Refactor-first approach:
  1. Remove GP dependencies from autoc (tree eval, bytecode, gpextractor). Replace GP config parser with something standalone. Prove with minisim, renderer, autoc using simplified autoc.ini.
  2. Get all remaining programs running (nnextractor, nn2cpp, etc).
  3. Update crrcsim (simplified RPC transport).
  4. Update xiao-gp (verify export pipeline).
  5. Repo cleanup — extract autoc/ to standalone repo.
  6. Tear out Boost (flexible ordering, may interleave with above).
  7. Resume neuro tuning work (sigma floor, curriculum, CMA-ES, etc).
  Keep the order flexible — Boost removal can shift around as convenient.

## Assumptions

- All three codebases (GP, CRRCSim, xiao-gp) are forking — no backward compatibility with GP-era code is required. The fork enables aggressive simplification and fast iteration.
- The surviving contract is sensor-in/control-out: the NN receives normalized sensor inputs and produces control commands. This contract spans training (autoc→sim), export (nn2cpp→xiao), and the NN01 serialization format.
- The existing VariationRampStep infrastructure from features 003/005 can be extended to support NN evolution without major architectural changes.
- GP flight data for behavioral cloning can be collected from existing S3 archives by replaying best-generation GP controllers through the simulator.
- The CRRCSim RPC protocol will be replaced with a simplified NN-only wire format; no transport versioning or multi-format detection needed.
- sep-CMA-ES with diagonal covariance (O(n) memory) is sufficient for the 531-dimensional weight space; full CMA-ES (O(n^2)) is not needed.
- The 10Hz state-action recording rate is sufficient for behavioral cloning; higher rates would increase data volume without proportional benefit.

## Scope Exclusions

- **NEAT / topology evolution**: Fixed topology (22,16,8,3) is correct for this problem; topology search is out of scope.
- **OpenAI-ES**: Requires 720+ distributed cores; not feasible with current infrastructure.
- **Full CMA-ES**: O(n^2) memory is wasteful at 531 dimensions; sep-CMA-ES provides the same benefit at O(n).
- **GP removal**: In-scope as prerequisite cleanup (deep dependency unwinding). Surviving interface contracts (RPC transport, S3 archiving, config system) are defined in this spec.
- **Policy gradient (REINFORCE)**: Deferred to future work; requires per-timestep streaming (FR-007) and possibly CRRCSim state injection as prerequisites.
- **CRRCSim state injection**: Restarting simulation from arbitrary mid-trajectory states is future research, not part of this feature.

## Dependencies

- **013-neuroevolution**: All NN infrastructure (NNGenome, nn_forward, NNPopulation, serialization, S3 archiving) from the current feature.
- **003-variations-redux / 005-entry-fitness-ramp**: VariationRampStep infrastructure for curriculum scenario ramping.
- **CRRCSim minisim**: RPC protocol changes needed for per-timestep fitness streaming (FR-007).

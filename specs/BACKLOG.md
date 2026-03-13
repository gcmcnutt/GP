# Project Backlog

**Location**: `~/GP/specs/BACKLOG.md`
**Source**: Migrated from `~/GP/autoc/TODO`
**Last Updated**: 2026-03-10 (migrated items to speckit features)

## Legend

- `[ACTIVE]` - Currently being worked on
- `[NEXT]` - High priority, ready to start
- `[BLOCKED]` - Waiting on dependencies
- `[DEFERRED]` - Lower priority, will revisit
- `[ABANDONED]` - Investigated and rejected
- `[DONE]` - Completed (reference only)
- `[SPEC]` - Transferred to speckit feature spec

---

## Infrastructure & Build

### [SPEC] Migrate Away from Boost Serialization → `specs/004-vacate-boost`

### [DEFERRED] Make pathgen.h Portable for Embedded
- Current: Separate `embedded_pathgen_selector.h` duplicates AeroStandard logic
- Issue: Code duplication, maintenance burden, version drift risk
- Goal: Single pathgen.h that works on both desktop and embedded
- Method: Conditional compilation for std::vector vs fixed arrays
- Defer: Wait until path system stabilizes

### [SPEC] Refactor Opcode Enum → `specs/008-opcode-enum`

### [SPEC] GP Library Fitness Serialization Precision → `specs/006-fitness-precision`

### [DONE] Separate data.dat/data.stc Filenames for Train vs Eval Mode → `012-distance-temporal-nodes`
- Eval mode now writes to `eval-data.dat` / `eval-data.stc`; training mode unchanged

### [DEFERRED] Unify S3 Archive Format for GP and NN Controllers
- Current: GP archives use Boost-serialized GP trees (`.dmp`); NN archives use custom binary "NN01" format (`.dmp`)
- S3 prefixes diverge: `autoc-*` (GP) vs `nn-*` (NN); both use INT64_MAX-epoch trick for descending time order
- Issue: Two distinct serialization formats creates a fork — renderer, extractor, CRRCSim all need format-aware branching
- Direction: GP path likely to be deprecated in favor of NN; consider a single envelope format (e.g., magic + type tag + payload) or converge on NN format only
- Defer: Wait until NN pipeline is validated end-to-end and GP deprecation decision is made

---

## Evolution & GP Engine

### [DEFERRED] Demetic Mode Elite Preservation
- Current: Best fitness jumps wildly (135→191→192→191→174→180→173→177)
- Root cause: Elite's aggregated fitness lost when re-evaluated on single scenario
- Impact: Cannot track evolution progress, GP search is ineffective
- Deferred: Hard to ensure with variable ramp or deme; revisit when fitness ramp stabilizes

### [DEFERRED] Co-evolve Random Paths
- Consider co-evolving the random paths alongside the GP population

### [DEFERRED] LLM Fitness Function Evolution
- Progressively raise the bar on complexity - staged evolution
- Could enable co-evolution of fitness criteria

### [DEFERRED] Simplify Eval Tree Interface
- Try to simplify the eval tree by passing a class reference instead of N args

### [NEXT] Checkpoint/Resume Evolution
- Dump full state at end of each generation so a crashed/paused run can resume from last completed gen
- State includes: population (GP trees + fitness), PRNG state, generation counter, elite store, scenario variation table
- Use local storage (not S3) for checkpoint files — fast writes, no network dependency
- Must restore deterministic PRNG sequence so resumed run produces identical results to uninterrupted run
- Motivation: long runs (200 gens × ~160s/gen ≈ 9 hours) are vulnerable to crashes, OOM, machine reboots
- Future: could also enable LLM-guided fitness function evolution between checkpoints

---

## Robustness & Variations

### [SPEC] Entry Variations with Ramped Fitness → `specs/005-entry-fitness-ramp`

### [DONE] Variations Project - Wind + Speed (003-variations-redux)
- Wind direction and speed variations validated in training

---

## Controller Architecture

### [SPEC] Layered Controller Safety Layer → `specs/010-safety-layer`

### [DEFERRED] Upper-Level Intercept Director
- Current GP controller assumes "intercept and track" — it just flies toward the rabbit
- With entry position/orientation variations, the aircraft may start pointed away, inverted, etc.
- A higher-level controller (or pre-intercept phase) should *judge* the appropriate initial
  maneuver based on current position and orientation before engaging track mode
- Examples: if starting 90° off-heading, first turn toward path; if starting inverted, first
  recover wings-level; if starting far displaced, fly direct intercept rather than chasing rabbit
- This is distinct from the current intercept-budget *fitness scaling* (005), which merely
  reduces penalties during intercept — it doesn't guide the aircraft's approach strategy
- Could be a separate GP program, a hand-coded state machine, or a learned policy
- Depends on: entry variation training (005) reaching maturity, safety layer (010) for guardrails

### [DONE] Path Interpolation → `specs/002-path-interpolation`

### [DEFERRED] Control Command Smoothness Fitness Objective
- Current fitness penalizes aircraft *attitude* changes but not *control command* changes
- GP *sometimes* exploits crrcsim's airframe inertia smoothing — does pulse duty-cycle control
- However, multiple runs with pop=20K and the expanded node set (23 nodes including
  temporal RATE/PREV nodes) have evolved smooth proportional controllers without any
  smoothness penalty. The bang-bang behavior appears seed-dependent, not systematic.
- **Observation**: With enough population size and the right node set, GP finds smooth
  solutions on its own. A smoothness penalty may be unnecessary overhead that constrains
  the search space. Revisit only if bang-bang becomes dominant when variations are enabled.

### [ACTIVE] Neuroevolution Controller → `specs/013-neuroevolution`
- Phases 1-7 complete (NN population, evaluator, fitness, serialization, S3, eval-mode, minisim integration)
- Phase 8 (polish/cleanup) remaining

### [DEFERRED] GP/NN Architecture Decision: Refactor or Fork
- Current: GP and NN code paths coexist with significant duplication (evalTask vs computeNNFitness, GPrand vs local RNG, parallel config parsing, format-aware branching in renderer/extractor/CRRCSim)
- Two paths forward:
  1. **Big refactor**: Unify GP and NN behind a clean ControllerBackend polymorphic interface — shared fitness, shared eval pipeline, shared serialization envelope, pluggable controller backends
  2. **Fork to NN-only**: Rip out GP tree evolution entirely, simplify codebase around NN-only workflow — bytecode/GP tree serialization, gpextractor, bytecode2cpp all become dead code
- Decision depends on: whether GP tree evolution provides unique value not achievable by NN (e.g., interpretability, modular subtree transfer), or whether NN supersedes GP for all use cases
- Defer until: 013-neuroevolution is fully validated end-to-end and comparative training runs (T098) inform the decision

### [SPEC] Fix LongSequential Path Immelman Segment → `specs/009-immelman-path`

### [SPEC] 4D Positional Fitness Surface → `specs/013-neuroevolution`
- Superseded by neuroevolution approach: NN can discover directional control from existing
  angular sensors (GETDPHI, GETDTHETA) × distance without encoding assumptions in fitness
- See 013 spec for analysis of why sensor-side representation is preferred over fitness encoding

### [DEFERRED] Error Cone for Future Path Points
- The further ahead we look from rabbit's current position, the less accurate the target point becomes
- Current evolved programs already show adaptive lookahead via nested GETDPHI

### [DEFERRED] Target Pose Estimation/Interpolation
- Current path interpolation is position-only
- Future controllers may need predicted target orientation (slerp quaternions)
- Depends on: Path interpolation implementation (002-path-interpolation)

---

## Embedded/Hardware Integration

### [NEXT] Training Record Consistency & Provenance
- Multiple S3 profiles/buckets exist (AWS default, minio, autoc-eval-arm) — need
  consistent provenance in all outputs so artifacts can be traced back to source
- **Fitness formatting**: `bytecode2cpp.cc` generated code comments use default
  `operator<<` for `header.fitness` (may show exponent notation). Should use
  `std::fixed << std::setprecision(6)` like gpextractor and gp_bytecode.cc do.
- **S3 provenance gaps** — these locations print S3 key without bucket/profile:
  - `autoc.cc:1982`: `"Results stored to S3: " << computedKeyName` — no bucket/profile
  - `gp_bytecode.cc:62`: `"S3 Key: " << getS3Key()` — no bucket/profile
  - `bytecode2cpp.cc:195,234`: generated comments have S3 Key but no bucket/profile
  - `gpextractor.cc:273`: prints fitness and S3 key but no bucket/profile
- **Good pattern**: `config_manager.cc:48-50` logs both bucket AND profile at startup
- **Fitness in bout stream**: `autoc.cc:846` prints fitness without explicit precision
- **Fix**: Add bucket+profile to GPBytecodeHeader (or a separate provenance struct),
  embed in bytecode files and generated code, use consistent fixed-point formatting
- Relevant for xiao-gp codebase (not yet in this repo)

### [DONE] Configurable Output File Prefixes for Train vs Eval → `012-distance-temporal-nodes`
- Auto-prefix `eval-` based on EvaluateMode in `autoc.cc`

### [NEXT] Export RC Commands to Xiao Log
- Currently: GP Output (rc=[...]) only logged during autoc=Y test spans
- Need: Log RC commands throughout entire flight for full playback visualization
- Location: xiao-gp/src/msplink.cpp - add logging in non-autoc code path

### [DEFERRED] GP to Autoc in INAV via Controls
- Selector mechanism
- Activate mechanism

---

## Scale & Performance

### [DONE] Unify Evaluation Pipelines → absorbed into `specs/013-neuroevolution` Phase 1

### [SPEC] GB10 GPU Native Evaluation → `specs/011-gpu-native`

---

## Visualization & Renderer

### [DEFERRED] Blackbox Rendering Improvements
- Select a path and a blackbox log for comparisons
- User interface to allow selecting path, log file, etc
- Additional attributes on renderer playback (e.g., loss of GPS signal)
- FPV mode for playback

### [DEFERRED] CRRCSim Display Dependency
- CRRCSim is dependent on a valid DISPLAY even when in headless mode

### [DEFERRED] Clean Shutdown Method
- Now that a crrcsim batch runs long, perhaps have polling loop for keepalive
- When autoc exits, crrcsim should exit cleanly

---

## Performance

### [DEFERRED] Wind Speed Variation
- Currently only wind direction varies across scenarios (via windDirectionOffset)
- Wind speed is fixed at base value from autoc_config.xml (e.g., 12 m/s)
- crrcsim adds Dryden turbulence on top (MIL-HDBK-1797), which provides continuous
  body-axis gust perturbations (~±10% of base wind speed, σ_w = 0.1 * V_wind)
- Different windSeeds produce different turbulence sequences, so GP already sees
  *stochastic* speed variation — but not *systematic* bias (e.g., 8 vs 16 m/s base)
- **Enhancement**: Add windSpeedOffset or windSpeedFactor to ScenarioMetadata
- Would need: new field in ScenarioMetadata, application in crrcsim windfield.cpp,
  sigma config in autoc.ini, generation in variation_generator.h
- This would train robustness to *different wind regimes*, complementing the
  turbulence which trains robustness to *gusts within a single regime*

### [NEXT] Batch and Cache Deterministic Scenarios
- With 36 wind scenarios, autoc serializes the full scenario table per individual evaluation
- This causes ~4x throughput hit (observed: 400% slower with 36 scenarios vs 1)
- Scenarios ARE deterministic per generation (seed + sigmas → offsets)
- **Optimization**: Send scenario table once at generation start, then only send population
  batches (GP trees) back and forth per individual
- crrcsim caches the scenario table and applies it to each incoming individual
- Reduces per-eval serialization from O(scenarios × individual) to O(individual)
- **Files**: `autoc/autoc.cc` (eval dispatch), `autoc/minisim.h` (ScenarioMetadata),
  crrcsim `inputdev_autoc.cpp` (scenario reception)

---

## Code Cleanup

### [DEFERRED] Memory Leak Investigation
- Walk various classes for copy by reference, copy by pointer, copy by value
- Small memory leak exists in autoc

### [DEFERRED] Buffer Best Run
- Buffer the best run over a population vs re-running it

### [DEFERRED] Interactive Mode Fallback
- Cannot interactively flip to autoc input as it attempts to use the callback port

### [DEFERRED] CRRCSim Robots for Lead Plane
- Figure out how to use crrcsim robots to form lead plane

### [DEFERRED] Check for Dangling Source Files
- Check for dangling or vestigial source files in the codebase

---

## Abandoned / Not Pursuing

### [ABANDONED] BLE Download Pipelining (PRN-Style)
- 1.42x gain doesn't justify complexity/reliability issues
- **Reference**: `autoc/specs/ZZZ-FLASH_SPEEDUP.md`

---

## Completed (Reference)

### [DONE] Variable Speed Rabbit — `autoc/specs/ZZZ-VARIABLE_RABBIT.md`
### [DONE] Path Interpolation — `specs/002-path-interpolation`
### [DONE] Wind + Speed Variations — `003-variations-redux`
### [DONE] GP Eval Tests — `specs/001-gp-eval-tests`

---

## Related Documentation

| Location | Description |
|----------|-------------|
| `~/GP/CLAUDE.md` | Main project guidance for Claude Code |
| `~/GP/autoc/specs/` | Design specs (ZZZ- prefix = archived/done) |
| `~/GP/specs/` | Feature specs (speckit workflow) |
| `~/GP/.specify/` | Speckit templates and configuration |

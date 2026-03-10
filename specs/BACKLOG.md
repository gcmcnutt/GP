# Project Backlog

**Location**: `~/GP/specs/BACKLOG.md`
**Source**: Migrated from `~/GP/autoc/TODO`
**Last Updated**: 2026-03-09 (003-variations-redux)

## Legend

- `[ACTIVE]` - Currently being worked on
- `[NEXT]` - High priority, ready to start
- `[BLOCKED]` - Waiting on dependencies
- `[DEFERRED]` - Lower priority, will revisit
- `[ABANDONED]` - Investigated and rejected
- `[DONE]` - Completed (reference only)

---

## Infrastructure & Build

### [DEFERRED] Migrate Away from Boost Serialization
- Boost serialization is NOT portable across different Boost versions
- Current issue: Boost 1.74 (x86) cannot read archives from Boost 1.83 (ARM64)
- Archive version in header changes between Boost releases (v18 vs v19)
- This breaks cross-architecture artifact sharing via S3
- Consider alternatives: cereal, JSON, FlatBuffers, custom text format
- Migration required for: EvalData, EvalResults, Path, AircraftState, ScenarioMetadata
- **Reference**: `autoc/specs/REPLACE_BOOST.md`

### [DEFERRED] Make pathgen.h Portable for Embedded
- Current: Separate `embedded_pathgen_selector.h` duplicates AeroStandard logic
- Issue: Code duplication, maintenance burden, version drift risk
- Goal: Single pathgen.h that works on both desktop and embedded
- Method: Conditional compilation for std::vector vs fixed arrays
- Defer: Wait until path system stabilizes

### [NEXT] Refactor Opcode Enum
- Need to refactor the opcode enum to only have one definition
- Stop using numeric opcodes in the bytecode2cpp generator

### [DONE] GP Library Fitness Serialization Precision
- ~~Renderer shows fitness as 164073.000 instead of 164073.136503~~
- Fixed in 003-variations-redux: `src/gp.cc` now uses `std::setprecision(17)` for double round-trip
- Note: Existing S3 archives still have truncated fitness; only new runs benefit

---

## Evolution & GP Engine

### [NEXT] Demetic Mode Elite Preservation
- Current: Best fitness jumps wildly (135→191→192→191→174→180→173→177)
- Root cause: Elite's aggregated fitness lost when re-evaluated on single scenario
- Impact: Cannot track evolution progress, GP search is ineffective
- Need: Preserve elite's full aggregated fitness across generations
- Approach: Fix attempted solution in autoc.cc lines 623-661, 665-744
  - Ensure elite is NOT re-evaluated (skip in evaluate() call)
  - OR restore aggregated fitness immediately after evaluation
  - Add logging to verify elite fitness is truly preserved
  - Test that best fitness is monotonically improving

### [DEFERRED] Co-evolve Random Paths
- Consider co-evolving the random paths alongside the GP population

### [DEFERRED] LLM Fitness Function Evolution
- Progressively raise the bar on complexity - staged evolution
- Could enable co-evolution of fitness criteria

### [DEFERRED] Simplify Eval Tree Interface
- Try to simplify the eval tree by passing a class reference instead of N args

### [DEFERRED] Checkpoint/Resume Evolution
- Checkpoint/resume a run - including for LLM enhancement over time

---

## Robustness & Variations

### [ACTIVE] Variations Redux - Entry/Wind/Speed Variations
- **Moved to**: `specs/003-variations-redux/` (feature branch `003-variations-redux`)
- Covers: crrcsim integration, fitness distance tuning, progressive training ramp
- Remaining for future VARIATIONS2:
  - Craft variations: CG, wing loading, control surface rates, drag, power/thrust
- **Reference**: `autoc/specs/ZZZ-VARIATIONS1.md`

---

## Controller Architecture

### [NEXT] Layered Controller Architecture
- **Reference**: `autoc/specs/LAYERED_CONTROLLER.md` for full design document

Current tracking layer (autoc-minimal.ini) follows the rabbit but has no awareness of:
- Am I about to crash?
- Am I going out of bounds?
- Am I in recovery vs tracking mode?

**SAFETY LAYER** (next priority):
- Hardcoded constraints that override GP outputs when danger detected
- OOB detection: if approaching geofence, break off and head toward origin
- Altitude protection: if too low, force climb
- Stall prevention: if too slow, force dive + throttle
- Re-engagement: when near origin/safe, resume tracking
- Implementation: post-GP safety clamping in evaluateOutput() or similar

**Testing**:
- Create intentional OOB paths to verify safety layer behavior
- Paths that go toward boundaries, into the ground, etc.
- Verify break-off and re-engagement logic works

**STRATEGY LAYER** (future):
- Phase detection (intercept → track → recover)
- Energy management
- Path preview / anticipation

### [DEFERRED] Control Command Smoothness Fitness Objective
- Current fitness penalizes aircraft *attitude* changes but not *control command* changes
- GP *sometimes* exploits crrcsim's airframe inertia smoothing — does pulse duty-cycle control
- However, multiple runs with pop=20K and the expanded node set (23 nodes including
  temporal RATE/PREV nodes) have evolved smooth proportional controllers without any
  smoothness penalty. The bang-bang behavior appears seed-dependent, not systematic.
- **Observation**: With enough population size and the right node set, GP finds smooth
  solutions on its own. A smoothness penalty may be unnecessary overhead that constrains
  the search space. Revisit only if bang-bang becomes dominant when variations are enabled.
- **If needed**: Add per-step control delta penalty matching existing pattern:
  - `control_delta = |Δroll_cmd| + |Δpitch_cmd| + |Δthrottle_cmd|`
  - `control_sum += pow(control_delta / CONTROL_NORM, CONTROL_POWER)`
  - CONTROL_NORM ~0.2 (20% of range per tick = 2.0/sec full travel)
  - CONTROL_POWER 1.5 (matches distance/attitude convention)
- **Related tuning**: DISTANCE_NORM/POWER tuning moved to `003-variations-redux`.
- **Files**: `autoc/autoc.cc` (fitness computation), `autoc/autoc.h` (constants)

### [NEXT] Fix LongSequential Path Immelman Segment
- Last segment of `longSequentialPath` generation should be an Immelman turn
- Current implementation has discontinuities due to coordinate convention issues
- Need to correct the final path section so the Immelman connects smoothly
- **File**: `autoc/pathgen.h` (AeroStandard::longSequentialPath)

### [DEFERRED] Error Cone for Future Path Points
- The further ahead we look from rabbit's current position, the less accurate the target point becomes
- Options evaluated:
  1. Weight by distance - MUL(GETDPHI(n), DIV(1, n))
  2. New operator GETDPHI_WEIGHTED(n)
  3. Keep it simple - For pure tracking layer, uncertainty may help
- Current evolved programs already show adaptive lookahead via nested GETDPHI

### [DEFERRED] Target Pose Estimation/Interpolation
- Current path interpolation is position-only
- Future controllers may need predicted target orientation (slerp quaternions)
- Would enable pose estimation sensors (target heading, target attitude)
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

### [NEXT] Export RC Commands to Xiao Log
- Currently: GP Output (rc=[...]) only logged during autoc=Y test spans
- Need: Log RC commands throughout entire flight for full playback visualization
- Location: xiao-gp/src/msplink.cpp - add logging in non-autoc code path
- Benefit: Control HUD in renderer 'a' mode will show stick/throttle for full flight

### [DEFERRED] GP to Autoc in INAV via Controls
- Selector mechanism
- Activate mechanism

---

## Visualization & Renderer

### [DEFERRED] Blackbox Rendering Improvements
- Select a path and a blackbox log for comparisons
- User interface to allow selecting path, log file, etc
- Additional attributes on renderer playback (e.g., loss of GPS signal)
- Blackbox logs are 1/32 Hz which is still a lot of data, consider subselect (5 Hz?)
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
- Save the final eval steps

### [DEFERRED] Interactive Mode Fallback
- Cannot interactively flip to autoc input as it attempts to use the callback port
- Soften this dependency - if port not set, fall back to Mouse

### [DEFERRED] CRRCSim Robots for Lead Plane
- Figure out how to use crrcsim robots to form lead plane

### [DEFERRED] Check for Dangling Source Files
- Check for dangling or vestigial source files in the codebase

---

## Abandoned / Not Pursuing

### [ABANDONED] BLE Download Pipelining (PRN-Style)
- Baseline: 1.77 KB/s with simple request-response (reliable)
- Attempted: Nordic DFU PRN-style pipelining (browser-driven)
- Result: 2.52 KB/s on 17.9KB file (1.42x), but 171KB wedges frequently
- Root cause: ArduinoBLE written() doesn't queue - rapid writeValueWithoutResponse calls overwrite each other
- Complexity: state tracking (outstanding/pending/drain/crc flags) has race conditions
- Conclusion: 1.42x gain doesn't justify complexity/reliability issues
- **Reference**: `autoc/specs/ZZZ-FLASH_SPEEDUP.md` "Phase 3: PRN-Style Pipelining"
- Future options: compression (Heatshrink), larger MTU, different BLE stack

---

## Completed (Reference)

### [DONE] Path Interpolation
- **Reference**: `specs/002-path-interpolation/`

### [DONE] Variable Speed Rabbit
- **Reference**: `autoc/specs/ZZZ-VARIABLE_RABBIT.md`

---

## Related Documentation

| Location | Description |
|----------|-------------|
| `~/GP/CLAUDE.md` | Main project guidance for Claude Code |
| `~/GP/autoc/specs/` | Design specs (ZZZ- prefix = archived/done) |
| `~/GP/specs/` | Feature specs (speckit workflow) |
| `~/GP/.specify/` | Speckit templates and configuration |

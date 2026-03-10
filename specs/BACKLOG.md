# Project Backlog

**Location**: `~/GP/specs/BACKLOG.md`
**Source**: Migrated from `~/GP/autoc/TODO`
**Last Updated**: 2026-03-06 (added Path Interpolation)

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

### [DEFERRED] GP Library Fitness Serialization Precision
- Renderer shows fitness as 164073.000 instead of 164073.136503
- Root cause: GP::save() in `/home/gmcnutt/GP/src/gp.cc:148-157` uses default ostream precision
- Fix: Add `os << std::setprecision(17) << stdFitness` (17 = max_digits10 for double)

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

### [NEXT] Variations Project - Zero-shot Transfer
- Goal: Get closer to zero-shot transfer from simulation to real aircraft
- Variations to implement:
  - Power, drag, service response times
  - ✅ [DONE] Varying speed rabbit (see `autoc/specs/ZZZ-VARIABLE_RABBIT.md`)
  - Craft variations (for robustness to airframe differences):
    - CG position variation
    - Wing loading variation
    - Control surface response times/rates
    - Drag coefficient variation
    - Power/thrust variation
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

### [NEXT] Path Interpolation for Smooth Target Tracking
- **Problem**: Current `getPathIndex()` snaps to discrete waypoints, causing issues:
  - **Timing jitter overshoot**: Waypoint at t=98ms is skipped when looking for t=100ms
    because `98 < 100`, jumping to t=200ms instead (discovered in gp_evaluator_tests)
  - **Discontinuous sensors**: GETDPHI/GETDTHETA jump discretely between waypoints
  - **Real-time sensitivity**: On xiao-gp, eval loop jitter causes erratic sensor values
- **Solution**: Interpolate target position based on exact time instead of snapping to waypoints
  - Binary search for bracketing waypoints around goal time
  - Linear interpolation between them: `pos = lerp(p0, p1, frac)`
  - Eval at t=99ms vs t=101ms gives nearly identical positions (smooth, not discrete jump)
- **Implementation**:
  - Add `getInterpolatedTargetPosition(pathProvider, currentTimeMsec, offsetSteps)`
  - Update `executeGetDPhi/DTheta/DTarget` to use interpolated positions
  - Remove discrete `getPathIndex()` function
- **Benefits**:
  - Robust to timing jitter (no more overshoot problem)
  - GP trains on smooth continuous tracking matching real-world behavior
  - Works correctly with variable rabbit speed
- **Files**: `gp_evaluator_portable.cc`, `aircraft_state.h`, `tests/gp_evaluator_tests.cc`

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
- **Related tuning**: Current DISTANCE_NORM=5.0 / DISTANCE_POWER=1.5 creates a "good enough"
  zone at 15-20m where GP stops trying to get closer. Consider tightening alongside smoothness
  if needed: e.g., NORM=3 POWER=2.0 would push harder toward tight tracking.
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

### [NEXT] Record S3 Profile in Extracted Artifacts
- `gpextractor` and `bytecode2cpp` print a timestamp but not which S3 profile was used
- Multiple S3 profiles exist (AWS default, minio) — need to record which one the
  extract came from so artifacts can be traced back to source
- Relevant for xiao-gp codebase (not yet in this repo)
- Simple fix: embed the profile string in the generated output header/comment

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

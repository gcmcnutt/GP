# FASTMATH Scalar Rollout Plan

This document expands on `specs/FASTMATH.md` by detailing how to remove all `float`/`double` math from the GP runtime, replacing it with a single scalar representation that works across training, simulator, and flight builds.

## 1. Goals and Constraints
- **Single scalar domain**: Every opcode executed by `evaluateGPOperator` (`GP/autoc/gp_evaluator_portable.cc:12-177`) and every stack slot in `evaluateBytecodePortable` (`GP/autoc/gp_evaluator_portable.cc:230-320`) must use the same fixed-precision scalar type. Legacy `double`/`float` usage is allowed only when ingesting external sensor/simulator inputs and when logging.
- **Scalar-native AircraftState**: Internally, `AircraftState`, `Path`, and all cached sensor values live entirely in the integer domain; minisim, Xiao/INAV, and CRRCSim perform conversions at their boundaries so the evaluator never pays per-opcode casts.
- **Eigen removal in hot paths**: Frequent sensors such as `GETALPHA`, `GETBETA`, `GETROLL_RAD`, and `GETPITCH_RAD` used to depend on Eigen quaternions/vectors (e.g., `Eigen::Quaterniond` in `aircraft_state.h:211-248`). Lightweight helpers in `fastmath/orientation_math.h` now cover those opcodes so the evaluator build no longer drags Eigen; the remaining Eigen usage lives at the simulator/firmware boundaries until they finish migrating their state containers.
- **Determinism and performance**: Controllers trained in `minisim` (`GP/autoc/minisim.cc`) and validated in CRRCSim (`~/crsim/crrcsim-0.9.13/.../inputdev_autoc.cpp`) must see the same scalar math that runs on the Xiao flight controller (`~/xiao-gp/src/msplink.cpp:312-360`). Expensive ops (e.g., trig, quaternion math) must be rewritten as LUT/CORDIC routines to handle thousands of opcode calls per evaluation without falling back to libm.
- **Accept hard saturation**: Sensor opcodes can clamp to their natural ranges; training from scratch will adapt to those limits, so we bias toward deterministic saturation over preserving float headroom.

## 2. Scalar Abstraction
1. **`gp_scalar` traits**  
   - `fastmath/gp_scalar.h` defines a templated `GPScalar<T>` wrapper with:
     - Underlying storage (`int32_t` in signed Q15.15 for production, optional `float` for desktop debugging). Q15.15 gives ~4–5 decimal digits of resolution above/below zero, matching the 10.10 fidelity target while leaving margin for command stacking.
     - Compile-time scale (`1 << 15`) and helpers (`fromDouble`, `toFloat`, `mul`, `div`, `saturatingAdd`, `narrowToCommand`).
     - A build flag `GP_ENABLE_FASTMATH` selecting the Q-format vs. float.
   - Provide constexpr constructors so literal opcodes (`OP_PI`, `ONE`, etc.) fold into integers during compilation.
   - Document the expected dynamic range so edge adapters know when to pre-clamp inputs before converting to `gp_scalar`.

2. **Math helpers**  
   - `fastmath/fixed_math.h` exports LUT-based `sinApprox`, `cosApprox`, `atan2Approx`, and integer `sqrtApprox`/`rsqrtApprox`. Each function consumes/produces `gp_scalar` and uses higher-precision accumulators (e.g., 64-bit) internally.
   - Range pre-conditioning (e.g., modulo 2π before sine lookup) ensures we stay within table bounds even when opcodes feed offsets repeatedly (per `specs/FASTMATH.md:13-20` strategy note).

3. **Vector/quaternion micro-library**  
   - `fastmath/orientation_math.h` replaces Eigen usage inside hot sensor helpers with small structs (`FastVec3`, `FastQuat`) storing scalars but exposing conversion helpers to/from the float domain.  
   - Inline routines (`quat_normalize`, `quat_rotate_vec`, `quat_to_euler`) use `gp_scalar` math when `GP_ENABLE_FASTMATH` is set, keeping compatibility with the Xiao build where Eigen is unavailable. Remaining Eigen references will be removed once the Xiao/crrcsim adapters finish their conversions.

## 3. Project-Specific Rollout

### 3.1 GP/autoc (core runtime)
1. **State containers**  
   - Convert `AircraftState` members (`velocity`, `position`, `pitchCommand`, etc.) to store `gp_scalar`. External loaders (minisim log replay, CRRCSim integration) translate host floats/doubles via `GPScalar::fromFloat` when writing into the state, and back to float only for logging. This keeps the evaluation core strictly integer while letting edge modules remain in their native formats.
   - Update `Path` to store positions/orientations using `gp_scalar`. Retain float serialization for Boost archives by converting at the boundary.

2. **Evaluator templating**  
   - Template `evaluateGPOperator` and `evaluateBytecodePortable` on the scalar type. Keep thin inline wrappers so existing call sites remain unchanged while the compiler instantiates the correct precision.
   - Implement opcode handlers using the helper math:
     - Pure arithmetic uses saturating integer ops.
     - Sensor ops (`GETVELX`, `GETVELY`, etc.) read `gp_scalar` directly from `AircraftState`.
     - Trig ops invoke the LUT functions, never Eigen.

3. **Bytecode generation and embedding**  
   - Update `bytecode2cpp.cc` so generated C arrays emit `gp_scalar` constants instead of floats (e.g., `constexpr gp_scalar embedded_gp_bytecode_constant = GPScalar::fromFloat(0.5f);`).
   - Ensure `gp_bytecode.h` uses scalar types for literal storage (`int32_t value` instead of `float constant` when FASTMATH is on).

4. **Testing harness**  
   - `gp_fastmath_tests` exercises opcode semantics, getters/setters, nav math, and evaluator stack flow. `ctest` (and therefore `make`) now builds and runs the tests automatically. Future work: add a dual-mode minisim regression that compares float vs. scalar outputs (<0.01 command delta) for full-program validation.

5. **Interpreter flattening**  
   - Phase 1 landed: `evaluateBytecodePortable` now handles arithmetic/logical/sensor opcodes inline with range-checked pushes instead of bouncing through `evaluateGPOperator`. Phase 2 (post-scalar validation) will explore per-program unrolling or macro expansion to erase the remaining switch/branch overhead.

### 3.2 xiao-gp (flight firmware)
1. **State ingestion**  
   - In `src/msplink.cpp:298-360`, convert all MSP/INAV telemetry into scalar form immediately (position, velocity, quaternions). Keep cached float copies only for telemetry/log output, so the evaluator sees pre-clamped Q15.15 inputs.
   - Replace direct Eigen calls with the shared `fastmath` quaternion helpers.

2. **Evaluator usage**  
   - Ensure `SinglePathProvider` and the cached `Path` segments store scalars. When refreshing `flight_path` from host uploads, convert once and discard float buffers.
   - Gate FASTMATH via build flag so flight builds always use the fixed-point backend; provide an optional debug menu to run the float reference for comparison if needed.

3. **Instrumentation**  
   - Record evaluation duration plus the number of high-cost opcodes per frame using the scalar backend. Store stats in `MSP_SEND_LOG` so we can confirm the <5 ms goal once fixed-point math replaces float/libm calls. When `GP_FASTMATH_TRACE` is set, reuse the shared metrics output (totals + per-eval averages, suppressed when all counters are zero).

### 3.3 CRRCSim integration
1. **Boundary translation**  
   - `inputdev_autoc.cpp` currently copies Eigen quaternions/positions straight into `AircraftState` (`lines 400-520`). Rework this step so CRRCSim values remain `double` until the final assignment, where we call `GPScalar::fromDouble`.
   - Remove Eigen from the CRRCSim build of `AircraftState` as well—use the shared quaternion helpers so simulator and firmware share identical math.

2. **Desktop evaluator**  
   - `gp_evaluator_desktop.cc/h` should instantiate the templated evaluator with the scalar type so CRRCSim training/evaluation exercises the same fixed-point logic as the flight build.
   - Provide a CLI flag (or env var `GP_FASTMATH_REFERENCE=float`) to temporarily run the float version for debugging regressions.

3. **Boost serialization**  
   - Update CRRCSim’s usage of Boost archives to serialize scalar values as integers. Provide version bumps in `Path`/`AircraftState` serialization to remain backward compatible with existing logs.

## 4. Migration Steps
1. **Instrumentation phase (DONE)**
   - Optional tracing (`GP_FASTMATH_TRACE`) in `evaluateBytecodePortable` captures saturation, clamp, and zero-injection counts per opcode. Logs emit total counts plus per-evaluation averages and stay quiet when all counters are zero. Use these stats from minisim, CRRCSim, and Xiao runs to confirm Q15.15 (or fallback 10.10) covers the observed ranges.

2. **Scalar scaffolding**
   - Land the `gp_scalar` traits, math helpers, and quaternion micro-library. Keep builds in float mode initially to avoid behavior changes while we wire up the new types.

3. **Core conversion**
   - Migrate `AircraftState`, `Path`, and the evaluator stack to the scalar type. Enable Q-format mode behind `GP_ENABLE_FASTMATH`, run dual-mode regression tests, and fix any overflow/saturation issues revealed by the instrumentation logs. Confirm that every opcode now touches only integer data between entry/exit conversions.

4. **Consumer rollout**
   - Update xiao-gp (MSP firmware) and CRRCSim integration to feed/consume scalars. Ensure their build systems pick up the new headers and drop Eigen from embedded targets.

5. **Validation**
   - Run:
     - `minisim` eval suites comparing float vs. fixed-point outputs.
     - CRRCSim autopilot evals to confirm determinism between Linux and firmware builds.
     - On-device timing tests (using the MSP logging already in `src/msplink.cpp`) to verify the <5 ms evaluation loop goal.

6. **Cleanup**
   - Remove unused Eigen includes/macros from `aircraft_state.h`.
   - Delete legacy float/double helper macros once all builds rely on `gp_scalar`.
   - Update documentation (`specs/FASTMATH.md`) with the final scalar format (Q15.15) and measured performance wins.
   - Schedule the evaluator inlining work after the scalar migration to reclaim the remaining call/branch overhead.

## 5. Open Questions
- Final Q-format selection (Q15.15 assumed; verify against instrumentation—downgrade to 10.10 only if saturation is excessive).
- Should we allow per-opcode scaling (e.g., store angles in radians vs. normalized units) to reduce saturation risk?
- How do we surface scalar overflow/underflow events during training so GP evolution can adapt rather than exploit undefined behavior?

Answering these during the implementation phases above will keep all environments—desktop training, CRRCSim validation, and Xiao flight hardware—in lockstep while delivering the targeted evaluation speedups.

## 6. Prototype Rollout Log
- **2025-12-04 – Scalar stack prototype**: Introduced `fastmath::GPScalar` (Q15.15) and rewired `gp_evaluator_portable` to keep its bytecode stack, arithmetic opcodes, and sensor outputs in integer space. Sensors and Eigen-heavy helpers still compute in `double` temporarily, but values now enter/exit the evaluator via `GPScalar`, eliminating per-opcode conversions.
- **2025-12-04 – API shim updates**: Updated `autoc-eval`, the desktop/embedded/CRRCSim evaluators, and the bytecode-to-C generator so legacy entry points still accept/return `double` while internally converting to `GPScalar`. This documents the contract change: anything calling `evaluateGPOperator`/`evaluateBytecodePortable` must now pass scalars (or run through the provided shims) instead of raw doubles.
- **2025-12-04 – Fast trig helpers**: Added `fastmath/fixed_math` with LUT-based `sin`/`cos`, a polynomial `atan2`, and a Newton-based `sqrt` so `gp_evaluator_portable` no longer calls libm directly for these hot-path operations. Bytecode generators now emit scalar stacks, and CMake links the new helpers everywhere the evaluator is built.
- **2025-12-04 – Metrics scaffolding**: Added optional `GP_FASTMATH_TRACE` counters (via `fastmath::getFastMathMetrics()`) so long runs can report how often `GPScalar` saturates, range limits clamp, or outputs get zero-snapped. This keeps instrumentation zero-cost unless explicitly enabled.
- **2025-12-04 – Metrics logging & Eigen-light sensors**: Enabled FastMath metric dumps in `minisim` and CRRCSim (under `GP_FASTMATH_TRACE`) and replaced the GETALPHA/BETA/ROLL/PITCH implementations with lightweight quaternion math that avoids Eigen matrix conversions. World→body rotations now use direct quaternion formulas, reducing per-opcode overhead.
- **2025-12-05 – Metrics normalization**: FastMath logs now include per-evaluation averages (total counts divided by evaluated paths) so CRRCSim/minisim logs remain short even over thousands of sims, making it easier to correlate saturation rates with GP size.
- **2025-12-05 – Metrics noise reduction**: Instrumentation skips printing when totals are zero, keeping CRRCSim logs readable even when saturation rarely triggers.
- **2025-12-05 – Inline opcode handling**: `evaluateBytecodePortable` now handles arithmetic/logical/sensor opcodes directly (with range limiting on every push) instead of re-entering `evaluateGPOperator`, reducing per-instruction overhead while keeping a fallback for uncommon operations.
- **2025-12-05 – Cached sensor state & unit tests**: `AircraftState` now caches body-frame velocity plus roll/pitch angles so hot opcodes don’t recompute Eigen transforms. Added `gp_fastmath_tests` to exercise opcode semantics and evaluator stack ordering via CTest.
- **Pending validation**:
  1. Run `minisim` regression suite comparing legacy float vs. new scalar evaluator (expect <0.01 command delta). Capture any saturation warnings.
  2. Exercise CRRCSim autopilot path with the scalar evaluator to confirm deterministic behavior before flashing Xiao builds.
  3. Update `xiao-gp` consumer to consume the new API (currently relying on implicit conversions) and measure on-device evaluation latency.

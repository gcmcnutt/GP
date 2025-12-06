# FASTMATH: fixed-point GP evaluation

We are replacing every floating-point value that flows through the GP evaluator with a deterministic fixed-point scalar so the interpreter runs quickly on desktop simulations, CRRCSim builds, and Xiao firmware. Long-running CRRCSim jobs were spending 100 ms+ per evaluation, which makes evolved controllers unusable on embedded hardware. The FASTMATH effort aligns all environments on the same math, removes per-opcode conversions, and gives us instrumentation to validate the new precision.

## Target state

- **Single scalar domain**: All opcodes and interpreter stack slots use `fastmath::GPScalar`, a signed Q15.15 fixed-point type (1 sign + 15 integer + 15 fractional bits). External interfaces (minisim logs, CRRCSim sensors, INAV telemetry) convert from their native floats/doubles into `GPScalar` once, before entering the evaluator. Outputs convert back only when sending commands to the simulator/firmware.
- **Integer-native AircraftState/Path**: `AircraftState`, `Path`, and cached sensor values store `GPScalar` internally. Expensive quaternion/velocity derivations happen once per simulator sample, and hot opcodes read the cached scalars directly.
- **Eigen-free hot paths**: Quaternion-to-Euler and body-frame velocity helpers now live in `fastmath/orientation_math.h` so we can drop Eigen from the evaluator build, particularly for Xiao. Sensors such as `GETALPHA/BETA` pull from the cached body velocity instead of recomputing Eigen transforms.
- **Deterministic math helpers**: `fastmath/fixed_math.{h,cc}` provides LUT/polynomial `sin`, `cos`, `atan2`, and `sqrt` implementations that accept/return `GPScalar`. These replace libm calls inside opcodes like `SIN`, `COS`, `ATAN2`, and `SQRT`.
- **Instrumentation-first**: Optional counters (`fastmath/metrics.h`, enabled with `-DGP_FASTMATH_TRACE`) track saturation, clamp, and zeroing events per evaluation. Minisim and CRRCSim can log aggregate stats (total and per-eval averages) every generation to confirm the Q15.15 range is sufficient.
- **Interpreter flattening**: `evaluateBytecodePortable` now handles the core arithmetic/logical/sensor opcodes inline, with range-checked pushes. The generic `evaluateGPOperator` fallback still exists for rare opcodes, but the common path no longer re-enters libm or Eigen.

## Performance goals

- Keep per-evaluation latency under 10 ms on CRRCSim/minisim hardware; 5 ms is the stretch goal needed for Xiao flight loops.
- Maintain simulator throughput of ≥100 sims/sec even with larger (500–600 node) GP individuals.
- Record saturation/clamp rates so training can adapt instead of silently operating near overflow.

## Implementation status (Dec 2025)

1. **Scalar stack + API shim**: `fastmath::GPScalar` (Q15.15) in place, evaluator stack retyped, and callers convert to/from scalars via shim functions. All arithmetic/logical opcodes now operate entirely in integer space.
2. **Math + sensors**: LUT trig (`sinApprox`, `cosApprox`, `atan2Approx`) and integer `sqrtApprox` are linked everywhere the evaluator builds. `AircraftState` caches body-frame velocity, roll, and pitch in scalar form, removing repeated Eigen math.
3. **Instrumentation**: `fastmath::FastMathMetrics` accumulates saturation/clamp counts. CRRCSim/minisim dump metrics only when non-zero to keep logs concise.
4. **Unit tests**: `gp_fastmath_tests` exercises every opcode (getters/setters, sensors, nav, constants) plus stack flow. `ctest` runs automatically as part of `make`, so regressions show up during CI.
5. **Build toggles**: Both GP/autoc and CRRCSim expose a `GP_FASTMATH_TRACE` CMake option and link `fastmath` helpers explicitly so fixed-point builds work from a clean checkout.

## Remaining roadmap

1. **Widen caching**: Precompute DPhi/DTheta inputs or other navigation-heavy values when the state updates; only recompute when the path provider advances. This should unlock additional sims/sec.
2. **Evaluator specialization**: After confidence in the scalar backend, experiment with per-program unrolled evaluators (fully flattening opcodes into inline functions) to remove the remaining switch/branch overhead.
3. **Xiao alignment**: Xiao firmware currently uses `bytecode2cpp` to generate C; we will migrate it to the shared interpreter once throughput is verified on desktop. The scalar backend already ensures the math matches.
4. **Extended validation**: Run long CRRCSim/minisim jobs (TRACE off) to gather throughput numbers, occasionally re-enabling TRACE to inspect saturation. Confirm Xiao builds hit the <5 ms target once integrated.

## References

- Detailed rollout steps: `specs/FASTMATH_SCALAR_ROLLOUT.md`
- Instrumentation docs: `fastmath/metrics.h`
- Unit tests: `tests/gp_evaluator_tests.cpp`
- Interpreter changes: `gp_evaluator_portable.cc`

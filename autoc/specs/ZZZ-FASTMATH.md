# Accelerate math operations during GP eval
Today we use a mix of expensive double, transcendental and generic EIGEN libraries which are all relatively expensive. For
example to execute them on the Seeed xiao nrf58240 we see 100msec for a 500 instruction GP.  We need to work back to something that
runs efficiently on the flight hardware.  Some examples:

## Cheatsheet summary (nRF52840 @ 64 MHz):

Integer (Cortex-M4)
- Simple ALU (add/logic): ~1 cycle → ~64 M/s
- Multiply: 1 cycle → ~64 M/s
- Divide: ~2–12 cycles → ~5–30 M/s

Float (float / f32, hardware FPU)
- Add / Mul / MAC: 1 cycle throughput → ~64 M/s
- Divide: ~14–20 cycles → ~3–5 M/s

Double (double / f64, software)
- Add / Mul: several dozen cycles → 5–10× slower than float
- Divide: hundreds of cycles → ~0.1–0.5 M/s range

Transcendentals (software, using float)
- sinf / cosf / expf / logf / atanf:
- ~1k–5k cycles each → ~10k–60k calls/sec

## Strategy
We need to make changes to several code bases
- Change ***all*** double to float around the evaluator -- including fitness functions, path generators, transports/serialization, and Eigen types (flip Eigen to float first, not double). Always stick to `gp_scalar` and cast legacy inputs at the edge.
- Examine places where we can share a divisor to minimise the use to divide
  - in some places we may want an accerated reciprocal operation if it winds up being better than divide
- Prefer multiply-by-inverse over divide where it’s safe/clear (e.g., use `* 0.01f` instead of `/ 100.0f` in coordinate conversions).
- Change all transcendentals and perhaps sqrt to lookup tables (no enable/disable flags; just switch to LUT-backed paths — YOLO it)
- Do direct minimal conversions of the quaternion calculations
- Convert all transports to and from evaluator to be float given we really don't even need to compute double precision at all (especially all the boost serializations)
- the general accuracy here is 2 or 3 digits so we don't need to carry around any additional precision beyond float
- we don't need to compare to earlier code as the evolved programs won't be portable anyway
- Validation/measurement:
  - Use gp_scalar/gp_vec3/gp_quat everywhere (Eigen now float-only); any legacy inputs get cast on entry rather than flowing as double
  - GTest functional checks cover every GP op (math, LUT trig, side-effects) using small synthesized values; no float-vs-double comparisons
  - We don’t need micro-benchmarks of the bytecode evaluator; just watch macro-level sims/sec (minisim ~150–200/s on laptop) and crash/nan counts
  - CRRCSim AUTOC plugin now uses float/LUT GP stack; desktop runs show ~50 sims/sec on 4c laptop (4–8 threads similar)

## Components
- ~/GP/autoc/
  - gp_evaluator_portable.cc is the root of the code we want 
  to speed up -- a lot of operations wind up referencing this code
  - aircraft_state.h will likely need some simplification too
  - minisim.cc is one of the simulators we test
  - crrcsim (see below) is the other one that is key for the sim to real work we are doing here.  most training at scale uses this.
  - xiao-gp is the flight implementation which today uses an unrolled evaluator hard coded as cpp code (bytecode2cpp.cpp).  For now this is ok, but we may go back to a faster bytecode evaluator
  - autoc.ini is the execution configuration we edit often for various training runs
- ~/crsim/crrcsim-0.9.13/
  - src/modinput_dev/inputdev_autoc is the interface library from the GP to the simulator
  - src/modinput_dev/inputdev_autoc is now float-only and uses the shared LUT ops; autoc_config.xml stays as-is for scenarios
- ~/xiao-gp/src
  - msplink.cpp is the primary interface to and from inav (which we won't change but is at ~/inav)
  - bytecode2cpp-generated program is already unrolled float code; main loop uses the LUT-backed evaluator helpers

## Steps
- Identify code in GP that needs to change and strike double everywhere possible (transcendentals, Eigen, transport, fitness)
- Stepwise rollout:
  - flip scalar/Eigen types to float everywhere (evaluator, aircraft state, transports, sim/glue)
  - replace transcendentals/sqrt with LUT-backed implementations
  - tighten Eigen usage or hand-roll equivalents where it helps
- First get things working with minisim, then with crrcsim, then finally with xiao; watch sims/sec instead of micro-benching the bytecode path
- Add Gtest-based functional checks for each node/op (add, divide, trig, side-effect nodes, stack order) comparing to the canonical float operations on small synthesized values; no float-vs-double comparisons
- Opcode dispatch note: interpreter still uses a switch/dispatch per opcode; xiao-gp’s generated code is already unrolled. If we chase more CPU wins, consider a jump-table or fully unrolled emitted C for desktop/minisim too, but weigh i-cache size (1050 steps) vs branch cost on nRF52840.

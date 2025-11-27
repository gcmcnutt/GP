# Profiling autoc and crrcsim

These steps capture high-resolution CPU samples and emit Flame Graphs so we can
see where time is spent inside both `autoc` and `crrcsim`.

## Prerequisites

1. `perf` (usually from the `linux-tools` package).
2. Brendan Gregg's [FlameGraph](https://github.com/brendangregg/FlameGraph)
   scripts cloned somewhere locally, e.g.
   ```bash
   git clone https://github.com/brendangregg/FlameGraph ~/FlameGraph
   ```
3. The kernel must allow unprivileged perf sampling. Either run `sudo sysctl
   kernel.perf_event_paranoid=1` (or lower) once per boot, or be prepared to
   enter your password when the helper script invokes `perf` via `sudo`.

Set `FLAMEGRAPH_DIR` to your clone location if it is not `~/FlameGraph`.

## Build with profiling flags

Both codebases now expose a `CMake` option that keeps frame pointers and emits
debug symbols so `perf` can unwind stack traces.

### `autoc`

```bash
cd ~/GP/autoc
cmake -B build -S . \
  -DENABLE_PROFILING=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

### `crrcsim`

```bash
cd ~/crsim/crrcsim-0.9.13
./cmake.sh -DENABLE_PROFILING=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

## Recording a single-worker baseline

1. Use one evaluation thread so the traces are easier to reason about. In
   `~/GP/autoc/autoc.ini` set `EvalThreads = 1` and pick a fixed port:
   `MinisimPortOverride = 5500`.
2. Start `crrcsim` manually under the profiler:
   ```bash
   cd ~/xiao-gp
   ./tools/record_flamegraph.sh crrcsim \
     ~/crsim/crrcsim-0.9.13/build/crrcsim \
     -g ~/crsim/crrcsim-0.9.13/autoc_config.xml \
     -p 5500 \
     -i AUTOC
   ```
   This produces `crrcsim.perf.svg`.
3. In another terminal run `autoc` under the same helper, pointing it at the
   same config file:
   ```bash
   cd ~/xiao-gp
   ./tools/record_flamegraph.sh autoc \
     ~/GP/autoc/build/autoc \
     -i ~/GP/autoc/autoc.ini
   ```
   This produces `autoc.perf.svg`.
4. Revert `EvalThreads`/`MinisimPortOverride` afterwards if needed.

The helper script stores the raw `perf.data`, the folded stacks, and the SVG
flame graph next to whichever output prefix you provide. Open the SVG in any
browser or VS Code's built-in viewer.

## Tips

- You can profile any other command (e.g. `minisim`, renderer) by pointing
  `record_flamegraph.sh` at the corresponding binary.
- To sample a running multi-worker session, launch `autoc` normally and run
  `sudo perf top` or `sudo perf record -p <PID>` (still benefitting from the
  profiling build flags).

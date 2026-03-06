# Profiling Guide: autoc and crrcsim

## Overview

This guide covers flame graph profiling for identifying bottlenecks in the autoc/crrcsim simulation pipeline. Understanding where time goes is critical before attempting GPU scale-up (see SCALEUP.md).

**Goal:** Determine actual serial vs parallel fraction for Amdahl's law analysis.

---

## Prerequisites

### Install perf (Linux/WSL2)

```bash
# WSL2 requires building perf from kernel source
# First, check current kernel version
uname -r  # e.g., 6.6.87.2-microsoft-standard-WSL2

# Option 1: Try installing pre-built (may not match kernel)
sudo apt-get install linux-tools-generic linux-cloud-tools-generic

# Option 2: Build perf from WSL2 kernel source (more reliable)
sudo apt-get install build-essential flex bison libelf-dev libdw-dev
git clone --depth 1 https://github.com/microsoft/WSL2-Linux-Kernel.git
cd WSL2-Linux-Kernel/tools/perf
make
sudo cp perf /usr/local/bin/

# Verify
perf --version
```

### Install FlameGraph tools

```bash
git clone https://github.com/brendangregg/FlameGraph.git ~/FlameGraph
export PATH=$PATH:~/FlameGraph
```

### Enable kernel perf events (WSL2)

```bash
# May need to adjust paranoid level
echo 1 | sudo tee /proc/sys/kernel/perf_event_paranoid

# For full call graphs (may need root)
echo 0 | sudo tee /proc/sys/kernel/kptr_restrict
```

---

## Building for Profiling

### autoc

```bash
cd ~/GP
rm -rf build && mkdir build && cd build

# Build with profiling flags (keeps frame pointers for call graph)
cmake -DENABLE_PROFILING=ON ../autoc
make -j$(nproc)
```

### crrcsim

```bash
cd ~/crsim/crrcsim-0.9.13
rm -rf build && mkdir build && cd build

# Build with profiling flags
cmake -DENABLE_PROFILING=ON ..
make -j$(nproc)
```

---

## Test Configuration

Create a minimal test configuration for profiling:

### ~/GP/autoc/autoc_profile.ini

```ini
# Profiling configuration - single worker, moderate load
[General]
SimProgram = ../crsim/crrcsim-0.9.13/build/crrcsim

[GP]
PopulationSize = 100
NumberOfGenerations = 10
Crossover = 0.9
Mutation = 0.1

[Extra]
evalThreads = 1           # Single worker for clean profiling
simNumPathsPerGen = 3     # Moderate path count
windScenarioCount = 2     # Some variation
```

---

## Profiling Procedure

### Step 1: Profile autoc (dispatcher side)

```bash
cd ~/GP/build

# Record with call graph (dwarf for C++, fp for frame pointer)
# Use -g for call graphs, -F for frequency
perf record -g -F 99 --call-graph dwarf -o autoc.perf.data \
    ./autoc -f ../autoc/autoc_profile.ini 2>&1 | tee autoc_run.log

# If dwarf doesn't work, try frame pointer mode
perf record -g -F 99 --call-graph fp -o autoc.perf.data \
    ./autoc -f ../autoc/autoc_profile.ini

# Generate report
perf report -i autoc.perf.data --stdio > autoc_perf_report.txt

# Generate flame graph
perf script -i autoc.perf.data | \
    ~/FlameGraph/stackcollapse-perf.pl | \
    ~/FlameGraph/flamegraph.pl --title "autoc CPU Flame Graph" > autoc_flame.svg
```

### Step 2: Profile crrcsim (worker side)

Since crrcsim is spawned by autoc, we need to profile it separately:

```bash
# Option A: Profile the child process directly
# Start crrcsim manually in headless mode with a test scenario

# Option B: Use perf to follow forks
perf record -g -F 99 --call-graph dwarf -o crrcsim.perf.data \
    --inherit -f -- \
    ./crrcsim -i ../documentation/examples/autoc_headless.xml

# Option C: Attach to running crrcsim
# First, find PID of crrcsim while autoc is running
ps aux | grep crrcsim
perf record -g -F 99 --call-graph dwarf -p <PID> -o crrcsim.perf.data -- sleep 30

# Generate flame graph
perf script -i crrcsim.perf.data | \
    ~/FlameGraph/stackcollapse-perf.pl | \
    ~/FlameGraph/flamegraph.pl --title "crrcsim CPU Flame Graph" > crrcsim_flame.svg
```

### Step 3: Profile Both Together (follow forks)

```bash
# This will capture autoc and all spawned crrcsim processes
perf record -g -F 99 --call-graph dwarf \
    --inherit \
    -o combined.perf.data \
    ./autoc -f ../autoc/autoc_profile.ini

# Generate combined flame graph
perf script -i combined.perf.data | \
    ~/FlameGraph/stackcollapse-perf.pl | \
    ~/FlameGraph/flamegraph.pl --title "autoc+crrcsim Combined" > combined_flame.svg

# Separate by process name
perf script -i combined.perf.data | grep -E "^autoc|^\s" | \
    ~/FlameGraph/stackcollapse-perf.pl | \
    ~/FlameGraph/flamegraph.pl --title "autoc Only" > autoc_only_flame.svg

perf script -i combined.perf.data | grep -E "^crrcsim|^\s" | \
    ~/FlameGraph/stackcollapse-perf.pl | \
    ~/FlameGraph/flamegraph.pl --title "crrcsim Only" > crrcsim_only_flame.svg
```

---

## Expected Bottleneck Areas

### autoc (dispatcher)

| Component | Expected % | Location |
|-----------|------------|----------|
| Serialization (boost) | 10-30% | `sendRPC<>`, `receiveRPC<>` in minisim.h |
| Fitness computation | 20-40% | `MyGP::evalTask()` fitness loop |
| Thread pool dispatch | 5-15% | `ThreadPool::submit()`, condition variables |
| GP evolution | 5-15% | `GP::nextGeneration()`, crossover, mutation |
| I/O (logging) | 5-10% | `fout <<` in evalTask |

### crrcsim (worker)

| Component | Expected % | Location |
|-----------|------------|----------|
| Aero calculations | 30-50% | `CRRC_AirplaneSim_Larcsim::aero()` |
| EOM integration | 20-30% | `EOM01::ls_step()`, `ls_accel()` |
| Wind/turbulence | 10-20% | `T_Wind::CalculateWind()`, Dryden model |
| GP evaluation | 10-20% | `eval_gp()` bytecode interpreter |
| Serialization | 5-15% | boost archive encode/decode |

---

## Analyzing Results

### Reading Flame Graphs

- **Width** = time spent (wider = more time)
- **Height** = call stack depth
- **Color** = arbitrary (helps distinguish functions)
- **Click** to zoom into a subtree
- **Search** (Ctrl+F) for specific function names

### Key Questions to Answer

1. **Serial fraction**: What % of autoc time is NOT in worker communication?
   - Look for: `ThreadPool`, `sendRPC`, `receiveRPC`
   - Everything else is serial overhead

2. **Simulation hotspots**: Where does crrcsim spend time?
   - Look for: `aero`, `ls_step`, `CalculateWind`, `eval_gp`
   - These are candidates for GPU porting

3. **Serialization overhead**: How much time in boost?
   - Search for: `boost::archive`, `serialize`
   - This motivates custom protocol (REPLACE_BOOST.md)

4. **Memory allocation**: Hidden overhead?
   - Look for: `malloc`, `operator new`, `vector::push_back`
   - Pre-allocation may help

### Calculating Amdahl's Limit

From the flame graphs, estimate:

```
T_serial = time in (evolution + dispatch + fitness_calc + logging)
T_parallel = time in (worker_communication + waiting)
T_total = T_serial + T_parallel

Serial_fraction = T_serial / T_total
Amdahl_limit = 1 / Serial_fraction
Max_speedup = Amdahl_limit × current_sims_per_sec
```

---

## Alternative: Manual Instrumentation

If perf doesn't work (common in WSL2), use manual timing:

### autoc instrumentation (add to autoc.cc)

```cpp
#include <chrono>
#include <atomic>

// Global timing accumulators
static std::atomic<uint64_t> g_serialize_ns{0};
static std::atomic<uint64_t> g_deserialize_ns{0};
static std::atomic<uint64_t> g_fitness_ns{0};
static std::atomic<uint64_t> g_dispatch_ns{0};
static std::atomic<uint64_t> g_evolution_ns{0};

// Timing macro
#define TIMED_BLOCK(accumulator) \
    auto _start = std::chrono::high_resolution_clock::now(); \
    auto _guard = [&]() { \
        accumulator += (std::chrono::high_resolution_clock::now() - _start).count(); \
    }; \
    struct _ScopeGuard { decltype(_guard) g; ~_ScopeGuard() { g(); } } _sg{_guard};

// Usage in evalTask:
void MyGP::evalTask(WorkerContext& context) {
    {
        TIMED_BLOCK(g_serialize_ns)
        sendRPC(*context.socket, evalData);
    }
    {
        TIMED_BLOCK(g_deserialize_ns)
        context.evalResults = receiveRPC<EvalResults>(*context.socket);
    }
    {
        TIMED_BLOCK(g_fitness_ns)
        // ... fitness computation loop ...
    }
}

// Print at end of run:
void printTimingStats() {
    printf("Timing breakdown:\n");
    printf("  Serialize:   %8.2f ms\n", g_serialize_ns / 1e6);
    printf("  Deserialize: %8.2f ms\n", g_deserialize_ns / 1e6);
    printf("  Fitness:     %8.2f ms\n", g_fitness_ns / 1e6);
    printf("  Dispatch:    %8.2f ms\n", g_dispatch_ns / 1e6);
    printf("  Evolution:   %8.2f ms\n", g_evolution_ns / 1e6);
}
```

### crrcsim instrumentation (add to fdm_larcsim.cpp)

```cpp
#include <chrono>
#include <atomic>

static std::atomic<uint64_t> g_aero_ns{0};
static std::atomic<uint64_t> g_ls_step_ns{0};
static std::atomic<uint64_t> g_wind_ns{0};
static std::atomic<uint64_t> g_gear_ns{0};
static std::atomic<uint64_t> g_engine_ns{0};
static std::atomic<uint64_t> g_aero_calls{0};

void CRRC_AirplaneSim_Larcsim::update(...) {
    auto wind_start = std::chrono::high_resolution_clock::now();
    env->CalculateWind(...);
    g_wind_ns += (std::chrono::high_resolution_clock::now() - wind_start).count();

    for (int n = 0; n < multiloop; n++) {
        auto step_start = std::chrono::high_resolution_clock::now();
        ls_step(dt);
        g_ls_step_ns += (std::chrono::high_resolution_clock::now() - step_start).count();

        auto aero_start = std::chrono::high_resolution_clock::now();
        aero(...);
        g_aero_ns += (std::chrono::high_resolution_clock::now() - aero_start).count();
        g_aero_calls++;

        // ... etc for engine, gear
    }
}

// Print on shutdown
void printFDMStats() {
    printf("FDM timing breakdown (%lu aero calls):\n", g_aero_calls.load());
    printf("  aero():    %8.2f ms (%.2f us/call)\n",
           g_aero_ns / 1e6, g_aero_ns / 1000.0 / g_aero_calls);
    printf("  ls_step(): %8.2f ms\n", g_ls_step_ns / 1e6);
    printf("  wind():    %8.2f ms\n", g_wind_ns / 1e6);
}
```

---

## Quick Start Commands

```bash
# 1. Build both projects with profiling
cd ~/GP && rm -rf build && mkdir build && cd build && cmake -DENABLE_PROFILING=ON ../autoc && make -j$(nproc)
cd ~/crsim/crrcsim-0.9.13 && rm -rf build && mkdir build && cd build && cmake -DENABLE_PROFILING=ON .. && make -j$(nproc)

# 2. Run profiled test (if perf works)
cd ~/GP/build
perf record -g -F 99 --call-graph dwarf -o profile.data ./autoc -f ../autoc/autoc.ini

# 3. Generate flame graph
perf script -i profile.data | ~/FlameGraph/stackcollapse-perf.pl | ~/FlameGraph/flamegraph.pl > flame.svg

# 4. Open in browser
# WSL2: explorer.exe flame.svg
# Linux: firefox flame.svg
```

---

## Next Steps After Profiling

1. **Document findings** in SCALEUP.md with actual percentages
2. **Calculate true Amdahl's limit** from measured serial fraction
3. **Prioritize optimizations**:
   - If serialization >20%: Implement custom protocol (REPLACE_BOOST.md)
   - If fitness >30%: Move to worker side
   - If aero/ls_step >50% of crrcsim: Primary GPU port targets
4. **Re-profile after each optimization** to verify improvement

---

## References

- [Brendan Gregg's Flame Graphs](https://www.brendangregg.com/flamegraphs.html)
- [perf wiki](https://perf.wiki.kernel.org/)
- [WSL2 perf setup](https://github.com/microsoft/WSL/issues/4297)

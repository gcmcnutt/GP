# Replace boost with a platform independent faster/bespoke serialization strategy
Right now we are using boost to communicate between the main autoc program and the simulators; ~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc and ~/GP/autoc/minisim.cc.  This library is relatively expensive and not compatible between ubuntu22/24 versions. We should replace this with a bespoke serialization strategy that is faster and platform independent.

---

## UPDATE LOG (Jan 2026)

**Current Status:** Boost binary_archive is working but with known limitations:
- Binary format is deterministic (text_archive was not)
- Cross-architecture (x86/ARM64) works if both are little-endian
- Boost version incompatibility (1.74 vs 1.83) breaks archive portability
- Performance is adequate but not optimal for scale-up

**Key Insight for Fast Strategy:** Given the expected scale-up (5,000-50,000 sims/sec target), the strategy must:
1. Minimize per-message overhead (zero-copy where possible)
2. Support caching of repeated data (paths, scenarios)
3. Be simple enough to maintain across x86/ARM64 only

**Recommended Path Forward:** Custom binary protocol with optional FlatBuffers migration path

---

## Strategy
Examine ~/GP/autoc/TODO to find the details of a boost removal.  From there study how it is used in autoc and the simulators.  Also notice how it is used to store build results used by renderer.cc and gpextrator.cc.  We want to remove boost completely from the current code and move on to another serialization format we can control over versions.  Also, we note that the serialization is a bottleneck in performance, so we want to make sure the new strategy is faster.

## Steps
- study use of boost
- propose alternate strategies
- first implement changes in ~/GP/autoc directory; autoc, renderer, minisim, gpextractor, etc if there are any others
- we will test all these end to end first
- then we can port the change to crrcsim
- secondary checks
  - what is being sent, given we want to greatly speed up the evals, we will likely need to reduce the routine data sent back and forth by caching what we can; paths that don't change, etc
  - perhaps a redistribution of what is a task to help group and reduce round trips to fewer, perhaps bulkier task specs

## verification
- cmake is used for GP
  - a clean build of GP is (cd ~/GP; make clean; rm -rf build; mkdir build; cd buildl cmake -DCMAKE_BUILD_TYPE=Debug ../autoc; cd ~/GP; make)
  - a simpler incremental build is (cd ~/GP; make) <- this is when the build system doesn't change
- cmake is used for crrcsim 
  - clean build (cd ~/crsim/crrcsim-0.9.13; rm -rf build; mkdir build; cd build; cmake -DCMAKE_BUILD_TYPE=Debug ..; make)
  - a simpler incremental build is (cd ~/crsim/crrcsim-0.9.13/build; make)
- actual testing is best done manually at this point.
  - modify ~/GP/autoc/autoc.ini - to use minisim, small runs, etc
  - run short runs
  - verify that the various reader tools work
- then later we will bring in crrcsim

## profiling
- today a typical autoc run (8 threads of execution) on this laptop shows
  - 250% cpu for autoc core
  - 8 crrcsim runs at 70% each
- Basically a significant amount of work is stuck in the autoc tasks, here's where maybe we do a bit more refactor or just send larger jobs to the task workers.  Sure the end of generation work should be done in autoc, and evolution, etc.  But right now, a lot of dispatch (or more likely needless extra copies of repeating data) are adding up to the scaling limit.
- Think ahead to when we have ported the simulator to big GPU systems.  Instead of laptop at 50 sim/second or DGX Spark at 500 sim/second, we want to get to GP10 at perhaps 5000-50000 sims/sec.  This is because we will greatly increase variations, population size, etc...

## Determinism Investigation Findings (2026-01-04)

### Root Cause Identified
Non-deterministic fitness values in GP evolution were caused by **boost::archive::text_archive** having non-deterministic float precision/formatting. Different text representations of the same floating-point values led to different deserialized results after network transmission.

### Investigation Process
1. **Verified identical inputs**: GP bytecode, initial simulation state, PRNG seeds all identical across workers
2. **Ruled out FPU configuration**: All workers had identical rounding modes (FE_TONEAREST), MXCSR=8099
3. **Ruled out compiler optimizations**: Issue persisted even with -O0
4. **Ruled out uninitialized heap memory**: MALLOC_PERTURB_=255 didn't crash
5. **Key breakthrough**: Hash verification showed:
   - Worker-side (crrcsim): 3 unique result hashes sent
   - Autoc-side: 36-50 unique result hashes received
   - **Serialization/deserialization was corrupting data**

### Solution Implemented (Jan 2026)
Switched from `text_archive` to `binary_archive` for deterministic serialization:

**Files Modified:**
- `~/GP/autoc/minisim.h`: Changed from text_archive to binary_archive with size-prefixed framing
- `~/GP/autoc/aircraft_state.h`: Added zero-initialization to AircraftState default constructor
- `~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp`:
  - Added `serializeForHash()` helper to hash serialized data instead of raw memory
  - Fixed padding issues by serializing before hashing

**Test Results:**
```
SUCCESS: All individuals have identical fitness (deterministic)
Group 1: fitness = 193.081543 (count=50)
```

### Cross-Platform Portability Notes
Current `binary_archive` implementation:
- ✅ **Works**: x86_64 ↔ x86_64, ARM64 ↔ ARM64, x86_64 ↔ ARM64 (both little-endian)
- ✅ **Assumption**: Modern systems use little-endian + IEEE 754 floats (true for x86/ARM cloud instances)
- ❌ **Won't work**: Little-endian ↔ big-endian cross-platform (rare in modern deployments)
- ✅ **Eigen types**: Already portable via element-wise serialization

**For guaranteed cross-endian support**, consider:
- Option 1: Implement explicit byte-swapping wrapper (htonl/ntohl for all primitives)
- Option 2: Use `eos::portable_binary_archive` (external library)
- Option 3: Custom portable format (protobuf, flatbuffers, capnproto)

### Remaining Cleanup Tasks
1. Remove debug logging from autoc.cc (AUTOC_RECV_HASH lines)
2. Consider disabling dtest_tracker in production builds
3. Remove MALLOC_PERTURB_=255 from crrcsim.sh (was for debugging)
4. Restore -O2 optimization in CMakeLists.txt (currently -O0 for testing)
5. Remove FPU state logging from dtest_tracker.h (not needed in production)

### Performance Considerations for Boost Replacement
When replacing boost serialization:
- Binary format is ~3x smaller than text (108 bytes vs 36 bytes for AircraftState vector)
- Size-prefixed framing adds 4 bytes overhead per message (negligible)
- Binary serialization is faster (no float→string conversion)
- Custom format could be even faster by removing boost overhead entirely
- Consider FlatBuffers/Cap'n Proto for zero-copy deserialization at scale

---

## Data Structure Analysis (Jan 2026)

### Current Serialized Types and Sizes

| Type | Fields | Est. Binary Size | Frequency | Cacheable? |
|------|--------|------------------|-----------|------------|
| **EvalData** (request) | | | Per-eval | Partially |
| → gp (bytecode) | vector<char> | ~200-2000 B | Per-eval | No |
| → gpHash | uint64 | 8 B | Per-eval | No |
| → isEliteReeval | bool | 1 B | Per-eval | No |
| → pathList | vector<vector<Path>> | ~50-200 KB | Per-gen | **Yes** |
| → scenarioList | vector<ScenarioMetadata> | ~50-100 B | Per-gen | **Yes** |
| **EvalResults** (response) | | | Per-eval | No |
| → crashReasonList | vector<CrashReason> | ~4-16 B | Per-eval | No |
| → aircraftStateList | vector<vector<AircraftState>> | ~100-500 KB | Per-eval | No |
| → debugSamples | vector<DebugSample> | 0 (usually) | Elite only | No |
| → physicsTrace | vector<PhysicsTraceEntry> | 0 (usually) | Elite only | No |
| **Path** | start, orient, dist, rad, time | 36 B | ~1000/path | Cacheable |
| **AircraftState** | idx, vel, quat, pos, cmds, time | ~80 B | ~1000/path | No |
| **DebugSample** | 40+ fields | ~400 B | Debug only | No |
| **PhysicsTraceEntry** | 60+ fields | ~600 B | Elite only | No |

### Key Optimization Opportunities

1. **Path caching is the big win** - pathList is 50-200KB per eval but changes only per-generation
   - With 100+ evals/gen, caching saves 5-20 MB/gen of serialization
   - Worker can cache paths by generation ID or path hash

2. **AircraftState dominates response size** - ~100-500KB per eval
   - 1000 steps × ~80 bytes × N paths
   - Can't cache, but can minimize fields or compress

3. **Debug/physics trace rarely sent** - only for elite re-eval
   - Already conditional, no change needed

### Proposed Wire Protocol Optimization

**Phase 1: Caching Protocol**
```
EvalRequest {
    uint32 pathGeneration;    // Generation number (for cache lookup)
    uint64 pathHash;          // Hash of pathList (cache key)
    bool includesPaths;       // If true, pathList follows
    PathList? pathList;       // Only sent if cache miss

    uint64 gpHash;
    bytes gpBytecode;
    ScenarioList scenarioList;
    bool isEliteReeval;
}
```

Worker maintains LRU cache of recent pathList by (pathGeneration, pathHash).
On cache miss, autoc sends full paths. On hit, worker uses cached paths.

**Expected savings:** 90%+ reduction in request size after first eval per generation.

---

## Serialization Library Comparison (Updated Jan 2026)

### Benchmark Data (from cpp-serializers and FlatBuffers benchmarks)

| Library | Serialize | Deserialize | Size | Zero-copy? | Header-only? |
|---------|-----------|-------------|------|------------|--------------|
| **FlatBuffers** | ~40 ns | ~80 ns | 1.0x | Yes | Yes |
| **Cap'n Proto** | ~30 ns | ~100 ns | 1.1x | Yes | No |
| **Protobuf** | ~200 ns | ~350 ns | 0.9x | No | No |
| **Boost binary** | ~150 ns | ~200 ns | 1.0x | No | No |
| **Custom memcpy** | ~10 ns | ~10 ns | 1.0x | Yes | Yes |

*Note: Benchmarks are indicative; actual performance depends on data structure.*

### Recommendation: Two-Phase Approach

**Phase 1: Custom Binary (Immediate)**
- Simple, fast, no dependencies
- Direct memcpy for POD structs (Path, AircraftState, etc.)
- Explicit versioning via magic number + version byte
- Works across x86/ARM64 little-endian (our only target platforms)

**Phase 2: FlatBuffers (Future)**
- If we need schema evolution, cross-language support, or stricter portability
- FlatBuffers is ~4x faster than Protobuf, has zero-copy deserialization
- Good tooling, active development, Google-backed

### Why NOT Cap'n Proto
- More complex to integrate (code generator, not header-only)
- Cross-platform support weaker than FlatBuffers
- Similar performance to FlatBuffers, less ecosystem

### Why NOT Protobuf
- Significantly slower than zero-copy alternatives
- Overkill for our x86/ARM64-only use case
- Good for interop, but we control both ends

---

## Fast Custom Protocol Design

### Design Principles

1. **Fixed-size headers, variable payloads** - parse header without reading entire message
2. **POD structs where possible** - direct memcpy, no per-field serialization
3. **Length-prefixed arrays** - uint32 count, then N × element
4. **Explicit version field** - allows schema evolution
5. **Little-endian everywhere** - native on x86/ARM64, no byte-swapping

### Wire Format

```
Message Header (12 bytes):
  uint32 magic;        // 0x47504556 "GPEV" for EvalData, 0x47505245 "GPRE" for EvalResults
  uint32 version;      // Protocol version (start at 1)
  uint32 totalSize;    // Total message size including header

EvalData Payload:
  uint32 gpSize;
  bytes gp[gpSize];
  uint64 gpHash;
  uint8 isEliteReeval;
  uint8 flags;          // bit 0: includesPaths
  uint32 pathGeneration;
  uint64 pathHash;
  [if includesPaths]:
    uint32 numPathGroups;
    for each group:
      uint32 numPaths;
      Path paths[numPaths];   // Direct memcpy of POD
  uint32 numScenarios;
  ScenarioMetadata scenarios[numScenarios];  // Direct memcpy of POD

EvalResults Payload:
  uint64 gpHash;
  uint32 numPaths;
  CrashReason crashReasons[numPaths];  // enum as uint8
  uint32 numStateGroups;
  for each group:
    uint32 numStates;
    AircraftState states[numStates];  // Direct memcpy of POD
  uint32 workerId;
  uint32 workerPid;
  uint32 workerEvalCounter;
  [if elite]:
    DebugSamples...
    PhysicsTrace...
```

### POD-ification Requirements

To enable direct memcpy, these structs must be POD-compatible:

**Path** (currently 36 bytes, already POD-like):
```cpp
struct PathPOD {
    float start[3];       // 12 bytes
    float orientation[3]; // 12 bytes
    float distanceFromStart; // 4 bytes
    float radiansFromStart;  // 4 bytes
    float simTimeMsec;       // 4 bytes
};  // 36 bytes, no padding needed
```

**AircraftState** (needs restructuring):
```cpp
struct AircraftStatePOD {
    float position[3];     // 12 bytes
    float velocity[3];     // 12 bytes
    float orientation[4];  // 16 bytes (quat: w,x,y,z)
    float relVel;          // 4 bytes
    float pitchCommand;    // 4 bytes
    float rollCommand;     // 4 bytes
    float throttleCommand; // 4 bytes
    float windVelocity[3]; // 12 bytes
    uint32_t simTimeMsec;  // 4 bytes (was unsigned long, wasteful)
    int32_t thisPathIndex; // 4 bytes
};  // 76 bytes, well-packed
```

**ScenarioMetadata** (currently ~40 bytes):
```cpp
struct ScenarioMetadataPOD {
    int32_t pathVariantIndex;   // 4 bytes
    int32_t windVariantIndex;   // 4 bytes
    uint32_t windSeed;          // 4 bytes
    uint64_t scenarioSequence;  // 8 bytes
    uint64_t bakeoffSequence;   // 8 bytes
    uint8_t enableDeterministicLogging; // 1 byte
    uint8_t padding[7];         // 7 bytes padding for alignment
};  // 36 bytes
```

### Implementation Helpers

```cpp
// Simple length-prefixed write
template<typename T>
void writeArray(std::vector<char>& buf, const std::vector<T>& arr) {
    uint32_t count = static_cast<uint32_t>(arr.size());
    buf.insert(buf.end(), reinterpret_cast<char*>(&count),
               reinterpret_cast<char*>(&count) + sizeof(count));
    if (count > 0) {
        const char* data = reinterpret_cast<const char*>(arr.data());
        buf.insert(buf.end(), data, data + count * sizeof(T));
    }
}

// Simple length-prefixed read
template<typename T>
std::vector<T> readArray(const char*& ptr) {
    uint32_t count;
    memcpy(&count, ptr, sizeof(count));
    ptr += sizeof(count);
    std::vector<T> arr(count);
    if (count > 0) {
        memcpy(arr.data(), ptr, count * sizeof(T));
        ptr += count * sizeof(T);
    }
    return arr;
}
```

---

## Task Batching for Scale-Up

### Current Model (1 eval = 1 RPC round-trip)
```
autoc → worker: EvalData (1 GP, N paths)
autoc ← worker: EvalResults (N path results)
```

At 50 sims/sec with 8 workers: 400 round-trips/sec
At 5000 sims/sec with 80 workers: 40,000 round-trips/sec (problematic)

### Proposed Model: Batched Evaluation

```
EvalBatch {
    uint32 batchId;
    uint32 pathGeneration;
    uint64 pathHash;
    PathList? pathList;  // Only on cache miss

    uint32 numEvals;
    EvalItem evals[numEvals];  // Each has gpBytecode + gpHash + scenarioList
}

EvalBatchResults {
    uint32 batchId;
    uint32 numResults;
    EvalItemResult results[numResults];  // Per-item crash/states
}
```

**Benefits:**
- Amortize TCP round-trip latency across batch
- Single path transmission per batch (caching)
- Worker can potentially parallelize within batch (GPU)
- Better utilization of network bandwidth

**Batch sizing:**
- CPU workers: 4-16 evals per batch (balance latency vs throughput)
- GPU workers: 64-256 evals per batch (maximize parallelism)

---

## Migration Plan

### Phase 1: Custom Binary Protocol (Week 1-2)

1. **Define POD structs** in new header `gp_wire_protocol.h`
2. **Implement serialize/deserialize** for EvalData, EvalResults
3. **Add path caching** to worker side
4. **Parallel deployment**: Run both boost and custom, compare results
5. **Switch over** once verified identical

### Phase 2: Remove Boost Serialization (Week 2-3)

1. **Update minisim.h** - remove boost archive includes
2. **Update autoc.cc** - use new protocol
3. **Update renderer.cc, gpextractor.cc** - if they read serialized files
4. **Update crrcsim** - use new protocol
5. **Remove boost::serialization** from CMakeLists.txt dependencies

### Phase 3: Task Batching (Week 3-4)

1. **Define batch protocol** in `gp_wire_protocol.h`
2. **Update worker** to handle batches
3. **Update autoc** task distribution to batch evals
4. **Benchmark** batch sizes for optimal throughput

### Phase 4: Optional FlatBuffers (Future)

If custom protocol becomes maintenance burden or cross-language needed:
1. Define `.fbs` schema files
2. Generate C++ bindings
3. Migrate incrementally (similar to Phase 1 parallel approach)

---

## Verification Checklist

- [ ] Identical fitness values between boost and custom protocol
- [ ] Path caching reduces request size by >90%
- [ ] No memory leaks in new serialization code
- [ ] Renderer can still read saved evaluation data
- [ ] gpextractor can still read saved evaluation data
- [ ] crrcsim works with new protocol
- [ ] Cross-architecture (x86/ARM64) verified

---

## UPDATE: Profiling Results Impact (Jan 2026)

**Key finding:** Boost serialization is only ~1.5% of runtime (see SCALEUP.md profiling results).

This changes the priority of this work:
- **Performance motivation: LOW** - serialization is not the bottleneck
- **Portability motivation: HIGH** - Boost version incompatibility remains painful
- **Scalability motivation: MEDIUM** - At 10,000+ sims/sec, even 1.5% matters

**Recommendation:** Defer full replacement until after quick wins (gear disable, SIMD).
Focus on **path/scenario caching** first as it provides bigger wins with less risk.

See SCALEUP.md "Scenario Caching & EvalTask Redesign" section for the caching architecture
that applies regardless of whether we use Boost or custom serialization.
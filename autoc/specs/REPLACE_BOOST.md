# Replace boost with a platform independent faster/bespoke serialization strategy
Right now we are using boost to communicate between the main autoc program and the simulators; ~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc and ~/GP/autoc/minisim.cc.  This library is relatively expensive and not compatible between ubuntu22/24 versions. We should replace this with a bespoke serialization strategy that is faster and platform independent.

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
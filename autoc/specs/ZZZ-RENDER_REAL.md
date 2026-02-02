# Flight Visualization Renderer

## Overview

VTK-based 3D renderer for visualizing:
- GP training runs (path vs actual flight trajectories)
- Blackbox flight logs from XIAO-GP hardware
- Test span navigation with highlighted segments

### Key Files

| File | Purpose |
|------|---------|
| `renderer.h` | VTK renderer class definition |
| `renderer.cc` | Main visualization implementation |

### Modes

| Mode | Description | Activation |
|------|-------------|------------|
| GP Training | Multiple arenas (up to 9 scenarios), path (green) vs actual (red) comparison, generation/fitness display | Default with S3 key |
| Xiao-Only | Single flight from hardware log, test span navigation, no generation navigation | `-x <logfile>` without S3 key |
| Decode | For analyzing bytecode execution | `-d` flag |

### Keyboard Controls

| Key | Function | Notes |
|-----|----------|-------|
| `n` | Next generation | GP mode only |
| `N` | Newest generation | GP mode only |
| `p` | Previous generation | GP mode only |
| `P` | Oldest generation | GP mode only |
| `t` | Next test span | All modes |
| `r` | Previous test span | All modes |
| `a` | All flight view | Shows absolute positioning with first arm at origin |
| `space` | Toggle playback | Animates through flight path |
| `f` | Toggle focus mode | Centers camera on current arena |
| `Left/Right/Up/Down` | Focus navigation | Moves between arenas (GP mode) |

### Animation System

```cpp
bool isPlaybackActive;
bool isPlaybackPaused;
std::chrono::steady_clock::time_point animationStartTime;
gp_scalar animationSpeed = 1.0f;
```

Playback animates through the flight path showing:
- Progressive path reveal
- Stopwatch overlay
- Control stick/throttle HUD

### Test Span Structure

```cpp
struct TestSpan {
  size_t startIndex;
  size_t endIndex;
  unsigned long startTime;
  unsigned long endTime;
  gp_vec3 origin;  // Test origin for xiao mode
  int pathIndex;   // Path index from GP State (0-5)
};
```

Test spans are extracted from `autoc=Y` segments in the xiao log.

---

## Debug Tips

### Check VTK Events
Add logging in `OnChar()`:
```cpp
std::cerr << "Key pressed: " << key << std::endl;
```

### Check Renderer Pointer
In CustomInteractorStyle:
```cpp
if (!renderer_) {
  std::cerr << "ERROR: renderer_ is null!" << std::endl;
  return;
}
```

### FPE Debug
Enable FPE trapping:
```cpp
#include <fenv.h>
feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
```

Then run under gdb to catch the exact location.

---

## Implementation Details: Xiao-Only Mode

This section describes the implementation for enabling `renderer.cc` to visualize real flight data from xiao-gp logs without needing external simulation data.

### Data Sources

**Xiao Flight Log** (e.g., `flight25a-hb1-flight_log_2026-02-01T21-43-53.txt`):
```
#00000328 0079654 0080327 i GP Control: Switch enabled - origin NED=[172.68, -107.99, -94.83] - program=autoc-...-gen9800.dmp
#00000329 0079655 0080327 i GP State: pos_raw=[172.68,-107.99,-94.83] pos=[0.00,0.00,0.00] vel=[-16.55,6.64,-1.97] quat=[-0.063,-0.095,0.403,-0.908] armed=Y fs=N servo=Y autoc=Y rabbit=Y path=2
#00000330 0079657 0080327 i GP Input: idx=1 rabbit=[-1.0,0.0,0.0] vec=[-1.0,0.0,0.0] alpha[0]=0.478 beta[0]=0.010 dtheta[0]=0.125 dphi[0]=1.004 dhome[0]=25.000 relvel=17.94
#00000331 0079663 0080327 i GP Output: rc=[1750,1406,1249]
```

**Blackbox CSV** (from INAV blackbox_decode):
- Provides high-frequency flight controller data (~100Hz)
- Key columns: `time (us)`, `navPos[0,1,2]`, `navVel[0,1,2]`, `quaternion[0-3]`, `rcData[0-3]`, `servo[0,1]`

---

## Phase 1: Fix Current `-x` Overlay ✅

### Changes Made

1. **Updated regex patterns** (lines 1790 and 2067):
   - Changed `test origin NED` → `origin NED` to match current log format

2. **Added control disable handling** (line 2067+):
   - Added `controlDisableRe` to detect span endings
   - Reset `inControlSpan = false` when disabled
   - Properly close any open autoc span on disable

3. **Fixed Z offset for xiao mode** (line 1700+):
   - Applied `SIM_INITIAL_ALTITUDE` (-25m) to virtual positions
   - Aligns xiao flight data with simulation path convention

---

## Phase 2: Standalone `-x` Mode (Xiao-Only Rendering) ✅

### Architecture

The renderer skips S3 operations when `-x` is provided without a simulation key. Key components:

1. **Mode Flag**: `inXiaoOnlyMode` bypasses `evalResults` entirely
2. **Rabbit Points**: Captured from GP Input lines for goal path reconstruction
3. **Rendering Path**: `renderFullScene()` handles xiao-only rendering
4. **Keyboard Adaptations**: Generation navigation (n/p/N/P) disabled

### Data Flow

```
Xiao Log File
    │
    ├─► parseXiaoData()
    │       ├─► GP State: pos, pos_raw, vel, quat → AircraftState (actual path)
    │       ├─► GP Input: rabbit → rabbitPoints (goal path)
    │       └─► GP Input: vec → craftToTargetVectors (arrows)
    │
    ├─► extractXiaoTestSpans()
    │       └─► autoc=Y/N transitions → TestSpan indices
    │
    └─► renderFullScene()
            ├─► rabbitPoints → paths (blue goal line)
            ├─► blackboxAircraftStates → blackboxTapes (yellow actual)
            └─► craftToTargetVectors → xiaoVecArrows (error arrows)
```

---

## Phase 3: Bug Fixes ✅

| Issue | Cause | Fix |
|-------|-------|-----|
| FPE on playback (SPACE) | Division by zero in `renderingOffset()` when `evalResults.pathList.empty()` | Early return with `vec3(0,0,0)` |
| 'f' key not working | `toggleFocusMode()` returned early on empty pathList | Added xiao-only mode handling |
| Arrow keys not working | `focusMove*()` functions returned early on empty pathList | Added xiao-only mode handling |
| Purple vec arrows cluttering view | Not needed in `-x` mode (blue error bars suffice) | Added `&& !inXiaoOnlyMode` guard |

---

## Phase 4: All-Flight Absolute Positioning ✅

### Implementation

1. **`showAllFlight()` modifications**:
   - In xiao-only mode, translates all raw positions so first arming point is at origin
   - Rebuilds `blackboxAircraftStates` with translated positions

2. **`renderFullScene()` modifications**:
   - When `showingFullFlight` is true:
     - Renders dimmed full flight tape (25% opacity)
     - Renders highlighted test span segments (full opacity)
     - Renders each span's path (rabbit points) at its absolute offset
     - Path offset = aircraft position at `autoc=Y` edge
     - Removes `SIM_INITIAL_ALTITUDE` from rabbit points in all-flight mode

3. **`parseXiaoData()` modifications**:
   - Now parses ALL GP State lines, not just ones during control spans
   - Removed `if (!inSpan) continue;` that was skipping GP State parsing
   - Added `gpStateLineToStateIndex` mapping for correct span indexing
   - Full flight path now continuous from takeoff to landing

4. **`updateBlackboxForCurrentTest()` modifications**:
   - Timestamps normalized to "time since arm" for proper playback

### Coordinate Handling

| Data | Storage | Notes |
|------|---------|-------|
| Raw positions (`pos_raw`) | `fullBlackboxAircraftStates` | Continuous through entire flight |
| Virtual positions (`pos`) | `xiaoVirtualPositions` | Reset at each arming |
| Rabbit points | Stored with `SIM_INITIAL_ALTITUDE` offset | Offset removed when rendering in all-flight mode |

---

## Testing Checklist

### Phase 1 ✅
- [x] Build renderer.cc with updated regex
- [x] Run with `-x flight25a-hb1-flight_log_2026-02-01T21-43-53.txt`
- [x] Verify test spans detected (5 spans found)
- [x] Verify span navigation works (t/r keys)
- [x] Verify Z offset applied (-25m)
- [x] Handle control disable properly

### Phase 2 ✅
- [x] Add `inXiaoOnlyMode` flag
- [x] Capture rabbit points in `parseXiaoData()`
- [x] Render goal path in `renderFullScene()`
- [x] Skip S3 operations in xiao-only mode
- [x] Adapt keyboard controls (disable n/p/N/P)
- [x] Update text display for xiao mode
- [x] Add ground plane for xiao-only mode

### Phase 3 ✅
- [x] FPE on playback fixed
- [x] Focus mode working
- [x] Arrow keys working
- [x] Purple arrows suppressed in -x mode

### Phase 4 ✅
- [x] All-flight shows continuous path from takeoff
- [x] First arm point at origin
- [x] Highlighted test spans
- [x] Paths rendered at absolute offsets
- [x] Timestamps normalized for playback
- [x] Progressive path animation in all-flight mode (paths appear when span starts, animate through duration)
- [x] Fixed timestamps to use absolute flight time (relative to arm, not span start) for continuous timeline

---

## Future Improvements

1. Add scrubbing (click timeline to jump)
2. Add speed control (+/- keys)
3. Export animation to video
4. IMU drift visualization
5. GPS quality overlay

---

## File Locations Reference

| File | Purpose |
|------|---------|
| `~/GP/autoc/renderer.cc` | Main visualization code |
| `~/GP/autoc/renderer.h` | Renderer class declaration |
| `~/xiao-gp/src/msplink.cpp` | Log format source (see `logPrint()` calls) |
| `~/xiao-gp/specs/POSTFLIGHT.md` | Log analysis workflow |
| `~/xiao-gp/specs/FLIGHT_LOG.md` | Log format documentation |
| `~/GP/autoc/specs/COORDINATE_CONVENTIONS.md` | NED coordinate system details |

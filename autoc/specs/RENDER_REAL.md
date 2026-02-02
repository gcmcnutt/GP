# Rendering Real Flight Data with renderer.cc

This document describes the plan for enabling `~/GP/autoc/renderer.cc` to fully visualize real flight data from xiao-gp logs, both for debugging the current `-x` overlay mode (Phase 1) and for a future mode that renders goal paths and flight data entirely from xiao logs without needing external playback data (Phase 2).

## Current State Analysis

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

### Current Issues with `-x` Mode (Phase 1)

The renderer.cc `-x` option uses `parseXiaoData()` to load xiao logs. **Test spans are not being detected** because of a regex mismatch:

#### Problem: Regex Pattern Outdated

**Current regex in renderer.cc** (line 1790):
```cpp
std::regex controlEnableRe(R"(GP Control: Switch enabled - test origin NED=\[([-0-9\.]+),\s*([-0-9\.]+),\s*([-0-9\.]+)\])");
```

**Actual log format** (from `xiao-gp/src/msplink.cpp:429`):
```cpp
logPrint(INFO, "GP Control: Switch enabled - origin NED=[%.2f, %.2f, %.2f] - program=%s", ...)
```

Key difference: The log format changed from `test origin NED=` to `origin NED=` and added `- program=...` suffix.

#### Additional Format Notes

Looking at the flight log, the GP State format also evolved:
- **Old format** (expected): `GP State:.*autoc=(Y|N)` with autoc flag at end
- **New format** (actual): `GP State: pos_raw=[...] pos=[...] vel=[...] quat=[...] armed=Y fs=N servo=Y autoc=Y rabbit=Y path=2`

The autoc flag regex `R"(autoc=(Y|N))"` should still work since it's a simple substring match.

---

## Phase 1: Fix Current `-x` Overlay

### Goal
Get the current `-x` mode working again to parse test spans and render xiao flight data as an overlay.

### Changes Required

#### 1. Update Control Enable Regex (renderer.cc:1790 and 2067)

**Before:**
```cpp
std::regex controlEnableRe(R"(GP Control: Switch enabled - test origin NED=\[([-0-9\.]+),\s*([-0-9\.]+),\s*([-0-9\.]+)\])");
```

**After:**
```cpp
// Updated regex to match current xiao-gp log format (changed from "test origin NED" to "origin NED")
std::regex controlEnableRe(R"(GP Control: Switch enabled - origin NED=\[([-0-9\.]+),\s*([-0-9\.]+),\s*([-0-9\.]+)\])");
```

This pattern will still work because:
- The `[...]` capture groups extract the NED origin coordinates
- The trailing `- program=...` is ignored (not matched, but not required)

#### 2. Verify Test Data

After the fix, run:
```bash
cd ~/GP/autoc
./renderer -x /mnt/c/Users/gcmcn/OneDrive/Documents/Navigation/HB1-orange-configs/flight25a-hb1-flight_log_2026-02-01T21-43-53.txt
```

Expected output should show detected test spans:
```
Found N test spans
Test 1: indices 329-500 (171 states) origin=[172.68, -107.99, -94.83] vecs=85
```

---

## Phase 1: Fix Current `-x` Overlay - COMPLETED

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

## Phase 2: Standalone `-x` Mode (Xiao-Only Rendering)

### Goal
Enable renderer.cc to render flight visualization **entirely from xiao log** without:
- S3 simulation data fetch
- Blackbox CSV decode
- External path files

### Current Architecture Issues

The renderer is deeply tied to `evalResults` from S3:

```cpp
// These all depend on evalResults.pathList:
int maxArenas = evalResults.pathList.size();           // Line 345
std::vector<vec3> p = pathToVector(evalResults.pathList[i]);  // Line 349
renderingOffset() uses evalResults.pathList.size()     // Line 97
```

The main flow (lines 1400-1464):
1. S3 lookup for latest run
2. Fetch generation file from S3
3. Deserialize into `evalResults`
4. `updateGenerationDisplay()` renders from `evalResults`

### Implementation Strategy

#### Option A: Synthetic `evalResults` (Minimal Changes)
Create a fake `evalResults` with one arena containing xiao-derived path.

**Pros**: Reuses existing rendering pipeline
**Cons**: Hackish, path reconstruction complex

#### Option B: Dedicated Xiao Rendering Path (Recommended)
Add `inXiaoOnlyMode` flag and bypass evalResults entirely.

**Pros**: Clean separation, no fake data
**Cons**: More code, some duplication

### Recommended Approach: Option B

#### 1. New Mode Flag and Validation

```cpp
// Global state
bool inXiaoOnlyMode = false;  // -x without simulation data

// In main(), after option parsing:
if (!xiaoLogFile.empty()) {
    if (!computedKeyName.empty()) {
        std::cerr << "Error: -x mode doesn't use -k (sim key)" << std::endl;
        return 1;
    }
    inXiaoOnlyMode = true;
}
```

#### 2. Skip S3 Operations in Xiao-Only Mode

```cpp
// Around line 1400:
if (!inXiaoOnlyMode) {
    // Existing S3 lookup code...
    if (computedKeyName.empty()) {
        // S3 folder listing...
    }
    // S3 generation fetch...
}
```

#### 3. Store Rabbit Points for Path Reconstruction

Add to `SpanData` struct:
```cpp
struct SpanData {
    vec3 origin;
    std::vector<TimestampedVec> vecs;
    std::vector<vec3> rabbitPoints;  // NEW: Goal path points
    size_t startStateIdx;
    size_t endStateIdx;
};
```

In `parseXiaoData()`, capture rabbit coordinates:
```cpp
// Parse GP Input (vec, relvel, AND rabbit)
if (std::regex_search(line, matches, inputRe)) {
    int idx = std::stoi(matches[1].str());

    // Capture rabbit position (already in virtual/origin-relative coords)
    scalar rabbit_x = std::stof(matches[2].str());
    scalar rabbit_y = std::stof(matches[3].str());
    scalar rabbit_z = std::stof(matches[4].str());
    vec3 rabbitPos(rabbit_x, rabbit_y, rabbit_z);

    // Apply same Z offset as actual positions
    rabbitPos[2] = rabbitPos[2] + SIM_INITIAL_ALTITUDE;
    currentSpanData.rabbitPoints.push_back(rabbitPos);

    // Existing vec handling...
}
```

#### 4. New Rendering Function for Xiao-Only Mode

```cpp
void Renderer::renderXiaoOnlyScene() {
    if (!inXiaoOnlyMode || testSpans.empty()) return;

    // Clear existing data
    paths->RemoveAllInputs();
    actuals->RemoveAllInputs();
    blackboxTapes->RemoveAllInputs();
    xiaoVecArrows->RemoveAllInputs();

    vec3 offset(0.0f, 0.0f, 0.0f);  // Single arena at origin

    // Get current span
    const TestSpan& span = testSpans[currentTestIndex];

    // Render goal path (blue) from rabbit points
    if (!xiaoSpanData[currentTestIndex].rabbitPoints.empty()) {
        paths->AddInputData(createPointSet(offset,
            xiaoSpanData[currentTestIndex].rabbitPoints));
    }

    // Render actual flight path (yellow tape)
    if (!blackboxAircraftStates.empty()) {
        std::vector<vec3> actualPoints = stateToVector(blackboxAircraftStates);
        blackboxTapes->AddInputData(createTapeSet(offset, actualPoints,
            stateToOrientation(blackboxAircraftStates)));
    }

    // Render vec arrows
    // ... existing vec arrow code ...

    // Create ground plane
    // ... single arena ground plane ...

    paths->Update();
    actuals->Update();
    blackboxTapes->Update();
    xiaoVecArrows->Update();

    renderWindow->Render();
}
```

#### 5. Adapt Keyboard Controls

Most controls work as-is. Only generation navigation (n/p/N/P) needs to be disabled.

| Key | Normal Mode | Xiao-Only Mode |
|-----|-------------|----------------|
| t/r | Next/prev test span | ✅ Same - works as-is |
| a | All flight view | ✅ Same - shows all spans |
| SPACE | Playback animation | ✅ Same - playback within span |
| f | Focus mode toggle | ✅ Same - single arena focus |
| Arrow keys | Move between arenas | ⚠️ N/A (single arena) |
| n/p | Next/prev generation | ❌ **Disable** - no generations |
| N/P | Newest/oldest gen | ❌ **Disable** - no generations |

In `CustomInteractorStyle::OnChar()` (renderer.h:229-240), add mode check:
```cpp
if (key == "n") {
    if (!inXiaoOnlyMode) {  // Skip generation nav in xiao-only mode
        this->InvokeEvent(NextModelEvent, nullptr);
    }
}
// Same for "N", "p", "P"
```

Note: Span detection uses `autoc=Y→N` transitions (not control enable/disable). The `controlDisableRe` just resets `inControlSpan` to stop counting GP State lines between control regions.

#### 6. Text Display Adaptation

```cpp
void Renderer::updateTextDisplay(...) {
    if (inXiaoOnlyMode) {
        // Show: "Xiao Flight - Span N of M"
        // Hide: Generation number, fitness
        // Show: Span duration, state count
    } else {
        // Existing display logic
    }
}
```

#### 7. Modified Main Flow

The existing xiao loading code (lines 1375-1396) already works. Main changes:
1. Skip S3 lookup when `-x` is provided alone
2. Call `renderXiaoOnlyScene()` instead of `updateGenerationDisplay()`
3. Print xiao-specific controls

```cpp
int main(int argc, char* argv[]) {
    // ... option parsing ...

    // Load xiao log data if specified (EXISTING CODE - lines 1375-1396)
    if (!xiaoLogFile.empty()) {
        if (!loadXiaoData()) return 1;
        renderer.inXiaoMode = true;
        extractXiaoTestSpans();
        // ... existing span setup ...
    }

    // NEW: Skip S3 operations in xiao-only mode
    if (!inXiaoOnlyMode) {
        // Existing S3 lookup code (lines 1400-1457)
        // ...
    }

    renderer.initialize();

    // NEW: Diverge rendering path
    if (inXiaoOnlyMode) {
        renderer.renderXiaoOnlyScene();

        // Print xiao-specific controls
        std::cout << "\nXiao Flight Mode Controls:" << std::endl;
        std::cout << "  t/r - Next/previous test span" << std::endl;
        std::cout << "  a - Show all flight" << std::endl;
        std::cout << "  SPACE - Playback animation" << std::endl;
        std::cout << "  f - Focus mode" << std::endl;
        std::cout << "  q - Quit" << std::endl;
    } else {
        // Existing: updateGenerationDisplay(), print sim controls
    }

    renderer.renderWindowInteractor->Start();
    return 0;
}
```

### Data Flow Summary

```
Xiao Log File
    │
    ├─► parseXiaoData()
    │       ├─► GP State: pos, pos_raw, vel, quat → AircraftState (actual path)
    │       ├─► GP Input: rabbit → rabbitPoints (goal path)  [NEW]
    │       └─► GP Input: vec → craftToTargetVectors (arrows)
    │
    ├─► extractXiaoTestSpans()
    │       └─► autoc=Y/N transitions → TestSpan indices
    │
    └─► renderXiaoOnlyScene()  [NEW]
            ├─► rabbitPoints → paths (blue goal line)
            ├─► blackboxAircraftStates → blackboxTapes (yellow actual)
            └─► craftToTargetVectors → xiaoVecArrows (error arrows)
```

### Files to Modify

| File | Changes |
|------|---------|
| `renderer.cc` | Add `inXiaoOnlyMode`, capture rabbit points, new render function, keyboard adaptations |
| `renderer.h` | Add `renderXiaoOnlyScene()` declaration |

### Benefits

1. **Self-contained**: Single xiao log file contains everything needed
2. **Clean separation**: Xiao mode doesn't pollute simulation rendering code
3. **Consistent visualization**: Same rendering primitives (paths, tapes, arrows)
4. **Easy debugging**: See exactly what GP saw vs what it commanded

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
- [x] Render goal path in `renderFullScene()` (reused existing function)
- [x] Skip S3 operations in xiao-only mode
- [x] Adapt keyboard controls (disable n/p/N/P)
- [x] Update text display for xiao mode (hide gen/fitness)
- [x] Add ground plane for xiao-only mode
- [ ] Test with flight25a log (5 spans)
- [ ] Verify goal path renders correctly
- [ ] Verify actual path aligns with goal path

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

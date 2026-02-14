# TEMPORAL STATE FOR GP NODES

## Overview

The GP tree is currently stateless - re-evaluated fresh each tick with no memory of prior evaluations. This limits the ability to do true integral control or trajectory prediction. This spec adds temporal terminals that enable:

1. **Error history** - What was GETDPHI/GETDTHETA N ticks ago?
2. **Error rate** - Rate of change of path errors (the D in PID)

## Related Documentation

- [TODO](../TODO) lines 111-141 - Original feature description
- [LAYERED_CONTROLLER.md](LAYERED_CONTROLLER.md) - Strategy/safety layer context

## New Terminals

| Terminal | Args | Description | Range |
|----------|------|-------------|-------|
| `GETDPHI_PREV` | 1 | GETDPHI(0) from N ticks ago (N=0..10) | [-π, π] |
| `GETDTHETA_PREV` | 1 | GETDTHETA(0) from N ticks ago (N=0..10) | [-π, π] |
| `GETDPHI_RATE` | 0 | Rate of change of GETDPHI(0) | rad/sec |
| `GETDTHETA_RATE` | 0 | Rate of change of GETDTHETA(0) | rad/sec |

All values are `gp_scalar` (float). History indices are converted via `static_cast<int>(arg)`.

### Semantic Clarification

**GETDPHI(n)** vs **GETDPHI_PREV(n)** - different `n` meanings:

| Call | `n` means | Returns |
|------|-----------|---------|
| `GETDPHI(0)` | lookahead=0 (current rabbit pos) | Current phi error to rabbit |
| `GETDPHI(2)` | lookahead=2 steps ahead | Phi error to future path point |
| `GETDPHI(-1)` | lookahead=-1 (behind rabbit) | Phi error to past path point |
| `GETDPHI_PREV(0)` | history=0 (this tick) | Most recent buffered GETDPHI(0) |
| `GETDPHI_PREV(1)` | history=1 tick ago | GETDPHI(0) from ~100ms ago |
| `GETDPHI_PREV(5)` | history=5 ticks ago | GETDPHI(0) from ~500ms ago |

The ring buffer stores **GETDPHI(0) computed each tick** - the instantaneous error to the rabbit's current position at that moment in time.

**GETDPHI_RATE**: Returns `(GETDPHI_PREV(0) - GETDPHI_PREV(1)) / dt`
- Zero-arg terminal (simpler than parameterized)
- Units: radians per second
- Positive = error increasing (getting worse)
- Negative = error decreasing (correcting)
- Clamped to [-10, 10] rad/sec to prevent extreme values

## Design Decisions

### Ring Buffer Size: 10 samples (1 second at 10Hz)

Rationale:
- At 100ms MSP/eval interval, 10 samples = 1 second of history
- Sufficient for I-term accumulation and D-term computation
- Memory cost: ~80 bytes per buffer (10 samples × 2 floats × 4 bytes)
- Embedded-friendly (XIAO BLE has ~256KB RAM)

### What Gets Buffered

Only the **instantaneous path errors at lookahead=0** are buffered:
- `dPhi_history[10]` - GETDPHI(0) values
- `dTheta_history[10]` - GETDTHETA(0) values
- `timestamp_history[10]` - Capture times in msec

Not buffered (compute on demand):
- GETDPHI at other lookaheads (would require N×M buffer)
- GETDTARGET (distance doesn't need temporal tracking)
- Position/velocity (already available from INAV)

### Recording Timing

History is recorded **before** GP evaluation each tick:
1. Compute GETDPHI(0), GETDTHETA(0) for current state
2. Push to ring buffer
3. Evaluate GP tree (which can access history)

This ensures the GP sees "fresh" previous values, not 1-cycle-stale data.

### Backwards Compatibility

- Old GPs without temporal terminals continue to work unchanged
- History buffer initialized to zeros on startup
- First 10 ticks will have incomplete history (graceful degradation)

## Implementation Checklist

### Phase 1: Core Infrastructure

#### 1.1 Enum Definition
**File**: `autoc.h` (after line 97)
```cpp
enum Operators {
  // ... existing ...
  GETDPHI, GETDTHETA, GETDTARGET, GETDHOME, GETVEL,
  GETDPHI_PREV, GETDTHETA_PREV,     // NEW: history lookback
  GETDPHI_RATE, GETDTHETA_RATE,     // NEW: error derivatives
  GETPITCH, GETROLL, GETTHROTTLE,
  // ... rest unchanged ...
};
```

#### 1.2 Node Definition Table
**File**: `autoc-eval.cc` (in allNodes[] array, around line 109)
```cpp
static const NodeDef allNodes[] = {
  // ... existing ...
  {GETDPHI, "GETDPHI", 1},
  {GETDPHI_PREV, "GETDPHI_PREV", 1},     // NEW: 1 arg (history index)
  {GETDPHI_RATE, "GETDPHI_RATE", 0},     // NEW: 0 args (nullary terminal)
  {GETDTHETA, "GETDTHETA", 1},
  {GETDTHETA_PREV, "GETDTHETA_PREV", 1}, // NEW: 1 arg (history index)
  {GETDTHETA_RATE, "GETDTHETA_RATE", 0}, // NEW: 0 args (nullary terminal)
  // ... rest ...
};
```

#### 1.3 History Buffer in AircraftState
**File**: `aircraft_state.h` (add to AircraftState class)

Note: All values use `gp_scalar` (float) for consistency with GP evaluation.

```cpp
class AircraftState {
public:
  // ... existing members ...

  // Temporal history for GP nodes
  static constexpr int HISTORY_SIZE = 10;  // 1 sec at 10Hz

private:
  gp_scalar dPhiHistory_[HISTORY_SIZE] = {0};    // gp_scalar = float
  gp_scalar dThetaHistory_[HISTORY_SIZE] = {0};
  unsigned long timeHistory_[HISTORY_SIZE] = {0}; // msec timestamps
  int historyIndex_ = 0;        // Next write position (ring buffer)
  int historyCount_ = 0;        // Valid samples (0 to HISTORY_SIZE)

public:
  // Record current errors to history (call before GP eval each tick)
  void recordErrorHistory(gp_scalar dPhi, gp_scalar dTheta, unsigned long timeMs) {
    dPhiHistory_[historyIndex_] = dPhi;
    dThetaHistory_[historyIndex_] = dTheta;
    timeHistory_[historyIndex_] = timeMs;
    historyIndex_ = (historyIndex_ + 1) % HISTORY_SIZE;
    if (historyCount_ < HISTORY_SIZE) historyCount_++;
  }

  // Get historical dPhi (n=0 is most recent, n=1 is one tick ago, etc.)
  // Returns 0.0f if history not available. Uses CLAMP_DEF for portability.
  gp_scalar getHistoricalDPhi(int n) const {
    if (historyCount_ == 0) return static_cast<gp_scalar>(0.0f);
    n = CLAMP_DEF(n, 0, historyCount_ - 1);
    int idx = (historyIndex_ - 1 - n + HISTORY_SIZE) % HISTORY_SIZE;
    return dPhiHistory_[idx];
  }

  gp_scalar getHistoricalDTheta(int n) const {
    if (historyCount_ == 0) return static_cast<gp_scalar>(0.0f);
    n = CLAMP_DEF(n, 0, historyCount_ - 1);
    int idx = (historyIndex_ - 1 - n + HISTORY_SIZE) % HISTORY_SIZE;
    return dThetaHistory_[idx];
  }

  unsigned long getHistoricalTime(int n) const {
    if (historyCount_ == 0) return 0;
    n = CLAMP_DEF(n, 0, historyCount_ - 1);
    int idx = (historyIndex_ - 1 - n + HISTORY_SIZE) % HISTORY_SIZE;
    return timeHistory_[idx];
  }

  int getHistoryCount() const { return historyCount_; }

  void clearHistory() {
    historyIndex_ = 0;
    historyCount_ = 0;
  }
};
```

**Memory footprint**: `10 * (4 + 4 + 4) + 8 = 128 bytes` per AircraftState instance.

### Phase 2: Evaluation Functions

#### 2.1 Forward Declarations
**File**: `gp_evaluator_portable.h` (after line 47)
```cpp
// Temporal terminals - PREV takes history index, RATE is nullary
gp_scalar executeGetDPhiPrev(AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDThetaPrev(AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDPhiRate(AircraftState& aircraftState);
gp_scalar executeGetDThetaRate(AircraftState& aircraftState);
```

#### 2.2 Main Switch Statement
**File**: `gp_evaluator_portable.cc` (in evaluateGPOperator switch, around line 154)
```cpp
// Temporal - PREV takes 1 arg (history index as gp_scalar, cast to int)
case GETDPHI_PREV:
    result = executeGetDPhiPrev(aircraftState, args ? args[0] : contextArg);
    break;
case GETDTHETA_PREV:
    result = executeGetDThetaPrev(aircraftState, args ? args[0] : contextArg);
    break;

// Temporal - RATE is nullary (no args)
case GETDPHI_RATE:
    result = executeGetDPhiRate(aircraftState);
    break;
case GETDTHETA_RATE:
    result = executeGetDThetaRate(aircraftState);
    break;
```

#### 2.3 Implementation Functions
**File**: `gp_evaluator_portable.cc` (after executeGetDTarget, around line 303)
```cpp
gp_scalar executeGetDPhiPrev(AircraftState& aircraftState, gp_scalar arg) {
    // arg is gp_scalar, convert to int history index
    int n = CLAMP_DEF(static_cast<int>(arg), 0, AircraftState::HISTORY_SIZE - 1);
    return aircraftState.getHistoricalDPhi(n);
}

gp_scalar executeGetDThetaPrev(AircraftState& aircraftState, gp_scalar arg) {
    int n = CLAMP_DEF(static_cast<int>(arg), 0, AircraftState::HISTORY_SIZE - 1);
    return aircraftState.getHistoricalDTheta(n);
}

gp_scalar executeGetDPhiRate(AircraftState& aircraftState) {
    // Rate = (current - previous) / dt
    gp_scalar current = aircraftState.getHistoricalDPhi(0);   // This tick
    gp_scalar previous = aircraftState.getHistoricalDPhi(1);  // Last tick

    // Time delta from history timestamps
    unsigned long t0 = aircraftState.getHistoricalTime(0);
    unsigned long t1 = aircraftState.getHistoricalTime(1);
    gp_scalar dt = (t0 > t1) ? static_cast<gp_scalar>(t0 - t1) / 1000.0f : 0.1f;

    if (dt < 0.001f) dt = 0.1f;  // Prevent divide by zero, default 100ms

    gp_scalar rate = (current - previous) / dt;
    return CLAMP_DEF(rate, -10.0f, 10.0f);  // Clamp to reasonable range (rad/sec)
}

gp_scalar executeGetDThetaRate(AircraftState& aircraftState) {
    gp_scalar current = aircraftState.getHistoricalDTheta(0);
    gp_scalar previous = aircraftState.getHistoricalDTheta(1);

    unsigned long t0 = aircraftState.getHistoricalTime(0);
    unsigned long t1 = aircraftState.getHistoricalTime(1);
    gp_scalar dt = (t0 > t1) ? static_cast<gp_scalar>(t0 - t1) / 1000.0f : 0.1f;

    if (dt < 0.001f) dt = 0.1f;

    gp_scalar rate = (current - previous) / dt;
    return CLAMP_DEF(rate, -10.0f, 10.0f);
}
```

### Phase 3: Bytecode Support

#### 3.1 Bytecode Generation
**File**: `gpextractor.cc` (in generateBytecode switch, around line 165)
```cpp
// PREV terminals: unary (1 arg = history index)
case GETDPHI_PREV:
case GETDTHETA_PREV:
    if (gene->containerSize() >= 1) {
        generateBytecode(gene->NthMyChild(0), program);
    }
    break;

// RATE terminals: nullary (0 args) - no child to generate
case GETDPHI_RATE:
case GETDTHETA_RATE:
    // Nothing to generate for children
    break;
```

#### 3.2 Stack Depth Analysis
**File**: `bytecode2cpp.cc` (in analyzeStackDepth, around line 45)
```cpp
// Unary operations (pop 1, push 1 = net 0)
case GETDPHI_PREV:
case GETDTHETA_PREV:
    // net 0
    break;

// Nullary terminals (push 1 = net +1)
case GETDPHI_RATE:
case GETDTHETA_RATE:
    currentStack += 1;
    break;
```

#### 3.3 Operator Name Mapping
**File**: `bytecode2cpp.cc` (in getOperatorName, around line 100)
```cpp
case GETDPHI_PREV: return "GETDPHI_PREV";
case GETDTHETA_PREV: return "GETDTHETA_PREV";
case GETDPHI_RATE: return "GETDPHI_RATE";
case GETDTHETA_RATE: return "GETDTHETA_RATE";
```

#### 3.4 Instruction Generation
**File**: `bytecode2cpp.cc` (in generateInstruction, around line 160)
```cpp
// PREV: unary (pop 1, push 1)
case GETDPHI_PREV:
case GETDTHETA_PREV:
    code << "    // " << opName << "\n";
    code << "    {\n";
    code << "        gp_scalar args[1] = {stack[sp-1]};\n";
    code << "        sp -= 1;\n";
    code << "        stack[sp++] = evaluateGPOperator(" << opcodeInt
         << ", pathProvider, aircraftState, args, 1, arg);\n";
    code << "    }\n";
    break;

// RATE: nullary (push 1)
case GETDPHI_RATE:
case GETDTHETA_RATE:
    code << "    // " << opName << "\n";
    code << "    stack[sp++] = evaluateGPOperator(" << opcodeInt
         << ", pathProvider, aircraftState, nullptr, 0, arg);\n";
    break;
```

#### 3.5 Bytecode Execution
**File**: `gp_evaluator_portable.cc` (in evaluateBytecodePortable switch, around line 345)
```cpp
// PREV: unary (pop 1, push 1)
case GETDPHI_PREV:
case GETDTHETA_PREV: {
    if (stack_ptr < 1) return 0.0f;
    gp_scalar args[1] = {stack[stack_ptr-1]};
    stack_ptr -= 1;
    stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider,
                                             aircraftState, args, 1, contextArg);
    break;
}

// RATE: nullary (push 1)
case GETDPHI_RATE:
case GETDTHETA_RATE:
    stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider,
                                             aircraftState, nullptr, 0, contextArg);
    break;
```

### Phase 4: Integration Points

#### 4.1 CRRCSim History Capture
**File**: `~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc/inputdev_autoc.cpp`

Before GP evaluation (around line 967), add:
```cpp
// Capture current errors to history before GP evaluation
{
    SinglePathProvider tempProvider(currentPath);
    gp_scalar dPhi = executeGetDPhi(tempProvider, aircraftState, 0.0f);
    gp_scalar dTheta = executeGetDTheta(tempProvider, aircraftState, 0.0f);
    aircraftState.recordErrorHistory(dPhi, dTheta, simTimeMsec);
}

// Now evaluate GP
interpreter->evaluate(aircraftState, path, 0.0);
```

#### 4.2 xiao-gp History Capture
**File**: `~/xiao-gp/src/msplink.cpp`

In `mspUpdateGPControl()` (around line 220), add before GP call:
```cpp
// Capture current errors to history before GP evaluation
{
    gp_scalar dPhi = executeGetDPhi(pathProvider, aircraft_state, 0.0f);
    gp_scalar dTheta = executeGetDTheta(pathProvider, aircraft_state, 0.0f);
    aircraft_state.recordErrorHistory(dPhi, dTheta, aircraft_state.getSimTimeMsec());
}

// Evaluate GP
generatedGPProgram(pathProvider, aircraft_state, 0.0f);
```

#### 4.3 Minisim History Capture
**File**: `~/GP/autoc/minisim.h` or evaluation loop in `autoc.cc`

Similar pattern - capture errors before each GP evaluation tick.

### Phase 5: Configuration

#### 5.1 Training Node Mask
**File**: `autoc.ini` (TrainingNodes line)
```ini
TrainingNodes = SETPITCH,SETROLL,SETTHROTTLE,GETDPHI,GETDTHETA,GETDTARGET,GETDPHI_PREV,GETDTHETA_PREV,GETDPHI_RATE,GETDTHETA_RATE,...
```

## Data Available from INAV

The temporal data we're adding is computable from what's already available:

| Data | INAV Source | Used For |
|------|-------------|----------|
| Position | MSP2_INAV_LOCAL_STATE pos[3] | GETDPHI/GETDTHETA calculation |
| Velocity | MSP2_INAV_LOCAL_STATE vel[3] | Could compute rate geometrically |
| Quaternion | MSP2_INAV_LOCAL_STATE q[4] | Body frame transformations |
| Timestamp | MSP2_INAV_LOCAL_STATE timestamp_us | dt calculation |
| Gyro rates | MSP_RAW_IMU (if enabled) | Alternative rate source |

The ring buffer approach captures computed path errors rather than raw IMU data, which is more directly useful for GP control.

## Testing Strategy

1. **Unit tests**: Verify ring buffer wrap-around, edge cases
2. **Replay tests**: Ensure history is maintained correctly across evaluation sequence
3. **Evolution test**: Run GP training with new terminals enabled, verify discovery
4. **Regression test**: Old GPs without temporal nodes still work correctly

## Future Extensions

If this works well, consider adding:

- `GETRELVEL_PREV(n)` / `GETRELVEL_RATE` - Velocity history and rate of change.
  **Revisit once variable target velocity is implemented.** With a non-constant
  velocity rabbit, velocity history becomes important for energy management
  and catching up/slowing down.
- `GETDTARGET_PREV(n)` - Distance history (less useful)
- `GETPOS_PREV(n)` - Position history (enables velocity change detection)
- `GETDPHI_INTEGRAL(n)` - Sum of last N GETDPHI values (I-term directly)
- `STORE(slot, value)` / `RECALL(slot)` - General purpose registers (complex)

## Files Changed Summary

| File | Changes |
|------|---------|
| `autoc.h` | Add 4 enum values |
| `autoc-eval.cc` | Add 4 entries to allNodes[] |
| `aircraft_state.h` | Add history buffer + methods |
| `gp_evaluator_portable.h` | Add 4 forward declarations |
| `gp_evaluator_portable.cc` | Add 4 case statements + 4 implementation functions |
| `gpextractor.cc` | Add 4 cases to bytecode generation |
| `bytecode2cpp.cc` | Add cases to 3 switch statements |
| `autoc.ini` | Add to TrainingNodes |
| `inputdev_autoc.cpp` (crrcsim) | Add history capture call |
| `msplink.cpp` (xiao-gp) | Add history capture call |

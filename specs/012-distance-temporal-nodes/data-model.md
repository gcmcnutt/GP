# Data Model: Distance Temporal Sensor Nodes

**Feature**: 012-distance-temporal-nodes
**Date**: 2026-03-12

## Entities

### Opcode Enum (extended)

Three new entries appended before `_END` marker in the Operators enum:

| Opcode | Name | Args | Returns | Units |
|--------|------|------|---------|-------|
| GETDIST | `"GETDIST"` | 0 (nullary) | Distance to rabbit | meters |
| GETDIST_PREV | `"GETDIST_PREV"` | 1 (unary: history index n) | Buffered distance | meters |
| GETDIST_RATE | `"GETDIST_RATE"` | 0 (nullary) | Rate of distance change | m/s |

Must appear in both `autoc.h` and `gp_evaluator_portable.h` enums in identical order.

### Distance History Buffer (new field in AircraftState)

| Field | Type | Description |
|-------|------|-------------|
| `distHistory_[HISTORY_SIZE]` | `gp_scalar[10]` | Ring buffer of raw distance values (meters) |

Shares with existing temporal infrastructure:
- `historyIndex_` — shared write pointer
- `historyCount_` — shared valid sample count
- `timeHistory_[HISTORY_SIZE]` — shared timestamp array (milliseconds)

### Buffer Lifecycle

```
New Path Start
  └─→ clearHistory()
       ├─ historyIndex_ = 0
       ├─ historyCount_ = 0
       ├─ zero dPhiHistory_[], dThetaHistory_[], distHistory_[]
       └─ zero timeHistory_[]

Each Simulation Tick (before GP evaluation)
  └─→ recordErrorHistory(dPhi, dTheta, distance, timeMs)
       ├─ dPhiHistory_[historyIndex_] = dPhi
       ├─ dThetaHistory_[historyIndex_] = dTheta
       ├─ distHistory_[historyIndex_] = distance      ← NEW
       ├─ timeHistory_[historyIndex_] = timeMs
       ├─ historyIndex_ = (historyIndex_ + 1) % HISTORY_SIZE
       └─ historyCount_ = min(historyCount_ + 1, HISTORY_SIZE)

GETDIST Query
  └─→ Returns current distance (computed fresh, not from buffer)

GETDIST_PREV(n) Query
  └─→ n = clamp(int(arg), 0, historyCount_ - 1)
       idx = (historyIndex_ - 1 - n + HISTORY_SIZE) % HISTORY_SIZE
       return distHistory_[idx]

GETDIST_RATE Query
  └─→ if historyCount_ < 2: return 0.0
       dist0 = distHistory_[(historyIndex_ - 1 + HISTORY_SIZE) % HISTORY_SIZE]
       dist1 = distHistory_[(historyIndex_ - 2 + HISTORY_SIZE) % HISTORY_SIZE]
       t0 = timeHistory_[(historyIndex_ - 1 + HISTORY_SIZE) % HISTORY_SIZE]
       t1 = timeHistory_[(historyIndex_ - 2 + HISTORY_SIZE) % HISTORY_SIZE]
       dt = (t0 > t1) ? (t0 - t1) / 1000.0 : 0.1
       if dt < 0.001: dt = 0.1
       return clamp((dist0 - dist1) / dt, -10.0, 10.0)
```

### Bytecode Instruction Format

Uses existing `GPBytecode` struct (9 bytes):

| Node | opcode | argc | constant |
|------|--------|------|----------|
| GETDIST | GETDIST | 0 | 0.0 |
| GETDIST_PREV | GETDIST_PREV | 1 | 0.0 |
| GETDIST_RATE | GETDIST_RATE | 0 | 0.0 |

Stack effects in bytecode interpreter:
- GETDIST: push 1 (nullary terminal)
- GETDIST_PREV: pop 1, push 1 (net 0 — unary operation)
- GETDIST_RATE: push 1 (nullary terminal)

### Node Registration Table Entry

| Opcode | Name String | Arg Count |
|--------|-------------|-----------|
| GETDIST | `"GETDIST"` | 0 |
| GETDIST_PREV | `"GETDIST_PREV"` | 1 |
| GETDIST_RATE | `"GETDIST_RATE"` | 0 |

### GETDTARGET Deprecation

GETDTARGET remains in:
- Opcode enum (unchanged position — backward compat)
- allNodes[] registration table (still parseable)
- All evaluator switch cases (still evaluates correctly)
- gpextractor / bytecode2cpp (still generates correctly)

GETDTARGET removed from:
- autoc.ini TrainingNodes (with deprecation comment)
- autoc-eval.ini TrainingNodes (with deprecation comment)

## Validation Rules

- History index n: cast to `int`, clamp to `[0, historyCount_ - 1]`
- Rate dt: must be >= 0.001s, default 0.1s
- Rate result: clamp to `[-10.0, 10.0]` m/s
- Distance: always >= 0 (Euclidean norm)
- Buffer reset: all values zero, count zero, index zero

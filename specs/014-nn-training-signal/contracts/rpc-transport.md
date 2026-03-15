# Contract: RPC Transport (autoc ↔ minisim)

## Overview

Binary length-prefixed messages over TCP. Single format for both train and eval modes.

## Wire Format

```
┌──────────────┬──────────────────────┐
│ uint32_t len │ payload (len bytes)  │
└──────────────┴──────────────────────┘
```

All multi-byte values are host byte order (same-machine communication only).

## Messages

### Request: EvalRequest (autoc → minisim)

```
magic: "EVAL"  (4 bytes)
version: uint32_t (1)
genome_count: uint32_t
for each genome:
  weight_count: uint32_t (531)
  weights: float[weight_count]
scenario_count: uint32_t
for each scenario:
  wind_speed: double
  wind_direction: double
  turbulence: double
  path_seed: uint64_t
  ... (scenario params TBD)
```

### Response: EvalResponse (minisim → autoc)

```
magic: "RSLT"  (4 bytes)
version: uint32_t (1)
genome_count: uint32_t
for each genome:
  scenario_count: uint32_t
  for each scenario:
    fitness: double
    crashed: uint8_t (0/1)
    timestep_count: uint32_t (0 if streaming disabled)
    for each timestep:
      time_ms: uint32_t
      distance: float
      dphi: float
      dtheta: float
      pitch_cmd: float
      roll_cmd: float
      throttle_cmd: float
```

## Constraints

- No Boost serialization — plain binary structs
- No versioned format negotiation — single version, both sides must match
- No GP tree or bytecode payloads — NN weights only
- Per-timestep data is optional (timestep_count=0 to disable)

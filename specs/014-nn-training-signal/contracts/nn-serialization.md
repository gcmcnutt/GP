# Contract: NN Serialization (NN01 format)

## Overview

Binary file format for NN genome persistence. Used for S3 archives, checkpoint files, and nn2cpp export pipeline.

## Format (existing, preserved)

```
magic: "NN01" (4 bytes)
layer_count: uint32_t
for each layer:
  size: uint32_t
weight_count: uint32_t
for each weight:
  value: float
fitness: double
generation: uint32_t
```

## S3 Archive Layout

```
nn-{timestamp}/
├── gen000.dmp    # Best genome of generation 0
├── gen001.dmp    # Best genome of generation 1
├── ...
└── genNNN.dmp
```

## Export Pipeline

```
S3 archive → nnextractor → nn_weights.dat → nn2cpp → nn_program_generated.cpp → xiao-gp
```

The sensor-in/control-out contract:
- Input: 22 normalized sensor values (see nn_topology.h)
- Output: 3 control commands (pitch, roll, throttle) via tanh → [-1, 1]

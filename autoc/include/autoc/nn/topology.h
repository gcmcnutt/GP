#ifndef NN_TOPOLOGY_H
#define NN_TOPOLOGY_H

// Compile-time NN topology constants (T110)
//
// Single source of truth for neural network dimensions.
// The topology is {22, 16, 8, 3}: 22 sensor inputs, two hidden layers
// (16 and 8 neurons), 3 control outputs (pitch, roll, throttle).
//
// Input layout (22 sensors):
//  0- 3: dPhi history    [now, -0.1s, -0.3s, -0.9s]  /π
//  4- 7: dTheta history  [now, -0.1s, -0.3s, -0.9s]  /π
//  8-11: dist history    [now, -0.1s, -0.3s, -0.9s]  /NORM_DIST
// 12-15: quaternion (w, x, y, z)                      [-1,1]
//    16: velocity                                     /NORM_VEL
//    17: alpha (angle of attack)                      /π
//    18: beta (sideslip)                              /π
// 19-21: pitchCmd, rollCmd, throttleCmd               [-1,1]
//
// Output layout (3 controls):
//    0: pitch command    [-1, 1] via tanh
//    1: roll command     [-1, 1] via tanh
//    2: throttle command [-1, 1] via tanh

#include <cstddef>

// Layer sizes
constexpr int NN_INPUT_COUNT = 22;
constexpr int NN_HIDDEN1_SIZE = 16;
constexpr int NN_HIDDEN2_SIZE = 8;
constexpr int NN_OUTPUT_COUNT = 3;

// Full topology as a compile-time array
constexpr int NN_NUM_LAYERS = 4;
constexpr int NN_TOPOLOGY[NN_NUM_LAYERS] = {
    NN_INPUT_COUNT, NN_HIDDEN1_SIZE, NN_HIDDEN2_SIZE, NN_OUTPUT_COUNT
};

// Total weight+bias count: (22*16+16) + (16*8+8) + (8*3+3) = 531
constexpr int NN_WEIGHT_COUNT =
    (NN_INPUT_COUNT * NN_HIDDEN1_SIZE + NN_HIDDEN1_SIZE) +
    (NN_HIDDEN1_SIZE * NN_HIDDEN2_SIZE + NN_HIDDEN2_SIZE) +
    (NN_HIDDEN2_SIZE * NN_OUTPUT_COUNT + NN_OUTPUT_COUNT);

// Topology as comma-separated string (for config logging/validation)
constexpr const char* NN_TOPOLOGY_STRING = "22,16,8,3";

#endif

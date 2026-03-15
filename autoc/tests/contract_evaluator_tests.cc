// Contract test: NN evaluator (sensor-in/control-out)
// Defines the surviving behavior: 22 normalized inputs → 3 control outputs in [-1,1]
// This contract MUST hold through GP removal, Boost removal, and source reorg.

#include <gtest/gtest.h>
#include "../nn_evaluator_portable.h"
#include "../nn_topology.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

// ============================================================
// Contract: NN topology is 22→16→8→3 with 531 weights
// ============================================================

TEST(ContractEvaluator, TopologyConstants) {
    EXPECT_EQ(NN_INPUT_COUNT, 22);
    EXPECT_EQ(NN_HIDDEN1_SIZE, 16);
    EXPECT_EQ(NN_HIDDEN2_SIZE, 8);
    EXPECT_EQ(NN_OUTPUT_COUNT, 3);
    EXPECT_EQ(NN_NUM_LAYERS, 4);
    EXPECT_EQ(NN_WEIGHT_COUNT, 531);
}

TEST(ContractEvaluator, WeightCountMatchesTopology) {
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    EXPECT_EQ(nn_weight_count(topology), NN_WEIGHT_COUNT);
}

// ============================================================
// Contract: Forward pass maps 22 inputs → 3 outputs in [-1,1]
// ============================================================

TEST(ContractEvaluator, ForwardPassDimensions) {
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    std::vector<float> weights(NN_WEIGHT_COUNT, 0.1f);
    float inputs[NN_INPUT_COUNT];
    float outputs[NN_OUTPUT_COUNT];

    // Fill inputs with normalized sensor-like values
    for (int i = 0; i < NN_INPUT_COUNT; i++) {
        inputs[i] = static_cast<float>(i - 11) / 11.0f;  // range [-1, 1]
    }

    nn_forward(weights.data(), topology, inputs, outputs);

    // All 3 outputs must be in [-1, 1] (tanh guarantee)
    for (int i = 0; i < NN_OUTPUT_COUNT; i++) {
        EXPECT_GE(outputs[i], -1.0f) << "Output " << i << " below -1";
        EXPECT_LE(outputs[i], 1.0f) << "Output " << i << " above 1";
        EXPECT_TRUE(std::isfinite(outputs[i])) << "Output " << i << " not finite";
    }
}

TEST(ContractEvaluator, ForwardPassOutputRange) {
    // Stress test: large weights should still produce outputs in [-1, 1]
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    std::vector<float> weights(NN_WEIGHT_COUNT);
    for (int i = 0; i < NN_WEIGHT_COUNT; i++) {
        weights[i] = static_cast<float>((i * 7 + 13) % 100 - 50) / 5.0f;
    }

    float inputs[NN_INPUT_COUNT];
    for (int i = 0; i < NN_INPUT_COUNT; i++) {
        inputs[i] = (i % 2 == 0) ? 1.0f : -1.0f;
    }
    float outputs[NN_OUTPUT_COUNT];

    nn_forward(weights.data(), topology, inputs, outputs);

    for (int i = 0; i < NN_OUTPUT_COUNT; i++) {
        EXPECT_GE(outputs[i], -1.0f);
        EXPECT_LE(outputs[i], 1.0f);
        EXPECT_TRUE(std::isfinite(outputs[i]));
    }
}

// ============================================================
// Contract: Forward pass is deterministic
// ============================================================

TEST(ContractEvaluator, ForwardPassDeterministic) {
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    std::vector<float> weights(NN_WEIGHT_COUNT);
    for (int i = 0; i < NN_WEIGHT_COUNT; i++) {
        weights[i] = static_cast<float>((i * 3 + 7) % 50 - 25) / 25.0f;
    }

    float inputs[NN_INPUT_COUNT];
    for (int i = 0; i < NN_INPUT_COUNT; i++) {
        inputs[i] = static_cast<float>(i) / static_cast<float>(NN_INPUT_COUNT);
    }

    float outputs1[NN_OUTPUT_COUNT], outputs2[NN_OUTPUT_COUNT];
    nn_forward(weights.data(), topology, inputs, outputs1);
    nn_forward(weights.data(), topology, inputs, outputs2);

    for (int i = 0; i < NN_OUTPUT_COUNT; i++) {
        EXPECT_EQ(outputs1[i], outputs2[i]) << "Non-deterministic at output " << i;
    }
}

// ============================================================
// Contract: Zero weights → zero outputs
// ============================================================

TEST(ContractEvaluator, ZeroWeightsZeroOutputs) {
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    std::vector<float> weights(NN_WEIGHT_COUNT, 0.0f);
    float inputs[NN_INPUT_COUNT];
    for (int i = 0; i < NN_INPUT_COUNT; i++) inputs[i] = 1.0f;
    float outputs[NN_OUTPUT_COUNT];

    nn_forward(weights.data(), topology, inputs, outputs);

    for (int i = 0; i < NN_OUTPUT_COUNT; i++) {
        EXPECT_NEAR(outputs[i], 0.0f, 1e-3);
    }
}

// ============================================================
// Contract: Xavier init produces valid genome
// ============================================================

TEST(ContractEvaluator, XavierInitProducesValidGenome) {
    NNGenome genome;
    genome.topology = std::vector<int>(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    nn_xavier_init(genome);

    EXPECT_EQ(static_cast<int>(genome.weights.size()), NN_WEIGHT_COUNT);

    // All weights must be finite
    for (float w : genome.weights) {
        EXPECT_TRUE(std::isfinite(w));
    }

    // Forward pass on Xavier-initialized weights must produce valid outputs
    float inputs[NN_INPUT_COUNT];
    for (int i = 0; i < NN_INPUT_COUNT; i++) inputs[i] = 0.5f;
    float outputs[NN_OUTPUT_COUNT];

    nn_forward(genome.weights.data(), genome.topology, inputs, outputs);

    for (int i = 0; i < NN_OUTPUT_COUNT; i++) {
        EXPECT_GE(outputs[i], -1.0f);
        EXPECT_LE(outputs[i], 1.0f);
        EXPECT_TRUE(std::isfinite(outputs[i]));
    }
}

// ============================================================
// Contract: fast_tanh matches std::tanh within tolerance
// ============================================================

TEST(ContractEvaluator, FastTanhAccuracy) {
    for (float x = -5.0f; x <= 5.0f; x += 0.05f) {
        gp_scalar result = testFastTanh(static_cast<gp_scalar>(x));
        gp_scalar expected = static_cast<gp_scalar>(std::tanh(x));
        EXPECT_NEAR(result, expected, 0.02)
            << "fast_tanh(" << x << ") = " << result << " vs std::tanh = " << expected;
    }
}

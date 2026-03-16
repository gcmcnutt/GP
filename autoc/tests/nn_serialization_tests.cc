#include <gtest/gtest.h>
#include "autoc/nn/serialization.h"
#include "autoc/nn/topology.h"
#include <cstring>

// ============================================================
// T043: Format detection tests
// ============================================================

TEST(NNFormatDetection, DetectsNNMagic) {
    uint8_t nn_data[] = {'N', 'N', '0', '1', 0, 0, 0, 0};
    EXPECT_TRUE(nn_detect_format(nn_data, sizeof(nn_data)));
}

TEST(NNFormatDetection, RejectsGPData) {
    uint8_t gp_data[] = {0x16, 0x00, 0x00, 0x00};
    EXPECT_FALSE(nn_detect_format(gp_data, sizeof(gp_data)));
}

TEST(NNFormatDetection, RejectsBytecodeData) {
    uint8_t bc_data[] = {0x16, 0x00, 0x04, 0x08, 0x04, 0x08};
    EXPECT_FALSE(nn_detect_format(bc_data, sizeof(bc_data)));
}

TEST(NNFormatDetection, RejectsTooShort) {
    uint8_t short_data[] = {'N', 'N', '0'};
    EXPECT_FALSE(nn_detect_format(short_data, sizeof(short_data)));
    EXPECT_FALSE(nn_detect_format(nullptr, 0));
}

TEST(NNFormatDetection, RejectsWrongMagic) {
    uint8_t wrong[] = {'N', 'N', '0', '2'};
    EXPECT_FALSE(nn_detect_format(wrong, sizeof(wrong)));
}

// ============================================================
// T042: Magic number test — serialized data starts with "NN01"
// ============================================================

TEST(NNSerialization, SerializedDataStartsWithMagic) {
    NNGenome genome;
    genome.topology = {2, 3, 1};
    genome.weights.resize(nn_weight_count(genome.topology), 0.5f);
    genome.fitness = 42.0;
    genome.generation = 7;
    genome.mutation_sigma = 0.1f;

    std::vector<uint8_t> buf;
    ASSERT_TRUE(nn_serialize(genome, buf));
    ASSERT_GE(buf.size(), 4u);
    EXPECT_EQ(buf[0], 'N');
    EXPECT_EQ(buf[1], 'N');
    EXPECT_EQ(buf[2], '0');
    EXPECT_EQ(buf[3], '1');
}

// ============================================================
// T041: Round-trip test
// ============================================================

TEST(NNSerialization, RoundTripPreservesAllFields) {
    NNGenome original;
    original.topology = std::vector<int>(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    original.weights.resize(nn_weight_count(original.topology));
    for (size_t i = 0; i < original.weights.size(); i++) {
        original.weights[i] = static_cast<float>(i) / 100.0f - 2.0f;
    }
    original.fitness = 123.456;
    original.generation = 42;
    original.mutation_sigma = 0.05f;

    std::vector<uint8_t> buf;
    ASSERT_TRUE(nn_serialize(original, buf));

    NNGenome restored;
    ASSERT_TRUE(nn_deserialize(buf.data(), buf.size(), restored));

    // Verify all fields match exactly
    EXPECT_EQ(restored.topology, original.topology);
    EXPECT_EQ(restored.weights.size(), original.weights.size());
    for (size_t i = 0; i < original.weights.size(); i++) {
        EXPECT_EQ(restored.weights[i], original.weights[i])
            << "Weight mismatch at index " << i;
    }
    EXPECT_DOUBLE_EQ(restored.fitness, original.fitness);
    EXPECT_EQ(restored.generation, original.generation);
    EXPECT_FLOAT_EQ(restored.mutation_sigma, original.mutation_sigma);
}

TEST(NNSerialization, RoundTripSmallNetwork) {
    NNGenome original;
    original.topology = {1, 1};
    original.weights = {0.5f, -0.3f};
    original.fitness = -999.0;
    original.generation = 0;
    original.mutation_sigma = 1.0f;

    std::vector<uint8_t> buf;
    ASSERT_TRUE(nn_serialize(original, buf));

    NNGenome restored;
    ASSERT_TRUE(nn_deserialize(buf.data(), buf.size(), restored));

    EXPECT_EQ(restored.topology, original.topology);
    EXPECT_EQ(restored.weights, original.weights);
    EXPECT_DOUBLE_EQ(restored.fitness, original.fitness);
}

TEST(NNSerialization, DetectFormatOnSerializedData) {
    NNGenome genome;
    genome.topology = {2, 2};
    genome.weights.resize(nn_weight_count(genome.topology), 0.0f);

    std::vector<uint8_t> buf;
    ASSERT_TRUE(nn_serialize(genome, buf));
    EXPECT_TRUE(nn_detect_format(buf.data(), buf.size()));
}

// ============================================================
// T044: Corrupt data test
// ============================================================

TEST(NNSerialization, RejectsTruncatedData) {
    NNGenome genome;
    genome.topology = std::vector<int>(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    genome.weights.resize(nn_weight_count(genome.topology), 1.0f);

    std::vector<uint8_t> buf;
    ASSERT_TRUE(nn_serialize(genome, buf));

    // Truncate to half
    NNGenome restored;
    EXPECT_FALSE(nn_deserialize(buf.data(), buf.size() / 2, restored));
}

TEST(NNSerialization, RejectsGarbageData) {
    uint8_t garbage[] = {0xFF, 0xFE, 0xFD, 0xFC, 0xFB, 0xFA};
    NNGenome restored;
    EXPECT_FALSE(nn_deserialize(garbage, sizeof(garbage), restored));
}

TEST(NNSerialization, RejectsEmptyData) {
    NNGenome restored;
    EXPECT_FALSE(nn_deserialize(nullptr, 0, restored));
}

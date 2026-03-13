#ifndef NN_SERIALIZATION_H
#define NN_SERIALIZATION_H

#include <vector>
#include <cstdint>
#include <string>
#include "nn_evaluator_portable.h"

// Magic bytes for NN binary format
constexpr char NN_MAGIC[4] = {'N', 'N', '0', '1'};
constexpr uint32_t NN_FORMAT_VERSION = 1;

// Serialize NNGenome to binary format
bool nn_serialize(const NNGenome& genome, std::vector<uint8_t>& out);

// Deserialize NNGenome from binary format
bool nn_deserialize(const uint8_t* data, size_t size, NNGenome& genome);

// Check if data starts with NN magic bytes
bool nn_detect_format(const uint8_t* data, size_t size);

#endif

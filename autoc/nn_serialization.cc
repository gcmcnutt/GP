#include "nn_serialization.h"
#include <cstring>

// Binary format layout:
// [0..3]   char[4]   magic "NN01"
// [4..7]   uint32_t  format version
// [8..11]  uint32_t  num_layers
// [12..]   uint32_t[num_layers]  topology
//          uint32_t  num_weights
//          float[num_weights]    weights
//          double    fitness
//          uint32_t  generation
//          float     mutation_sigma

namespace {

template<typename T>
void write_val(std::vector<uint8_t>& buf, const T& val) {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(&val);
    buf.insert(buf.end(), p, p + sizeof(T));
}

template<typename T>
bool read_val(const uint8_t*& ptr, size_t& remaining, T& val) {
    if (remaining < sizeof(T)) return false;
    std::memcpy(&val, ptr, sizeof(T));
    ptr += sizeof(T);
    remaining -= sizeof(T);
    return true;
}

} // namespace

bool nn_serialize(const NNGenome& genome, std::vector<uint8_t>& out) {
    out.clear();

    // Magic
    out.push_back('N');
    out.push_back('N');
    out.push_back('0');
    out.push_back('1');

    // Version
    write_val(out, NN_FORMAT_VERSION);

    // Topology
    uint32_t num_layers = static_cast<uint32_t>(genome.topology.size());
    write_val(out, num_layers);
    for (int layer_size : genome.topology) {
        uint32_t ls = static_cast<uint32_t>(layer_size);
        write_val(out, ls);
    }

    // Weights
    uint32_t num_weights = static_cast<uint32_t>(genome.weights.size());
    write_val(out, num_weights);
    for (float w : genome.weights) {
        write_val(out, w);
    }

    // Metadata
    write_val(out, genome.fitness);
    write_val(out, genome.generation);
    write_val(out, genome.mutation_sigma);

    return true;
}

bool nn_deserialize(const uint8_t* data, size_t size, NNGenome& genome) {
    if (!data || size < 4) return false;

    const uint8_t* ptr = data;
    size_t remaining = size;

    // Check magic
    if (ptr[0] != 'N' || ptr[1] != 'N' || ptr[2] != '0' || ptr[3] != '1') {
        return false;
    }
    ptr += 4;
    remaining -= 4;

    // Version
    uint32_t version;
    if (!read_val(ptr, remaining, version)) return false;
    if (version != NN_FORMAT_VERSION) return false;

    // Topology
    uint32_t num_layers;
    if (!read_val(ptr, remaining, num_layers)) return false;
    if (num_layers == 0 || num_layers > 100) return false;  // sanity check

    genome.topology.resize(num_layers);
    for (uint32_t i = 0; i < num_layers; i++) {
        uint32_t ls;
        if (!read_val(ptr, remaining, ls)) return false;
        genome.topology[i] = static_cast<int>(ls);
    }

    // Weights
    uint32_t num_weights;
    if (!read_val(ptr, remaining, num_weights)) return false;
    if (remaining < num_weights * sizeof(float)) return false;  // check before allocating

    genome.weights.resize(num_weights);
    for (uint32_t i = 0; i < num_weights; i++) {
        if (!read_val(ptr, remaining, genome.weights[i])) return false;
    }

    // Metadata
    if (!read_val(ptr, remaining, genome.fitness)) return false;
    if (!read_val(ptr, remaining, genome.generation)) return false;
    if (!read_val(ptr, remaining, genome.mutation_sigma)) return false;

    return true;
}

bool nn_detect_format(const uint8_t* data, size_t size) {
    if (!data || size < 4) return false;
    return data[0] == 'N' && data[1] == 'N' && data[2] == '0' && data[3] == '1';
}

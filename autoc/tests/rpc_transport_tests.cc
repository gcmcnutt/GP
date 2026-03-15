// Contract test: RPC transport
// Defines the surviving behavior: serialize EvalRequest/EvalResponse, round-trip,
// verify fields match. This contract defines the NEW wire format from
// contracts/rpc-transport.md (post-Boost binary protocol).
//
// IMPORTANT: Wire format is LITTLE-ENDIAN for cross-platform portability
// between x86 (little-endian) and ARM32/ARM64 (configurable, usually LE).
// All multi-byte values are explicitly serialized in little-endian order.
// IEEE 754 float/double representation is assumed (standard on all targets).
//
// NOTE: These tests define the TARGET format. They will initially test the
// serialization contract in isolation (struct→bytes→struct round-trip).
// The actual Boost→binary migration happens in Phase 3.

#include <gtest/gtest.h>
#include <vector>
#include <cstring>
#include <cstdint>

// ============================================================
// Portable little-endian serialization helpers
// Wire format: all integers and IEEE 754 floats in little-endian byte order
// ============================================================

namespace rpc_contract {

// Write uint32_t as 4 bytes, little-endian
inline void write_le_u32(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back(static_cast<uint8_t>(v & 0xFF));
    buf.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
    buf.push_back(static_cast<uint8_t>((v >> 16) & 0xFF));
    buf.push_back(static_cast<uint8_t>((v >> 24) & 0xFF));
}

// Write uint64_t as 8 bytes, little-endian
inline void write_le_u64(std::vector<uint8_t>& buf, uint64_t v) {
    for (int i = 0; i < 8; i++) {
        buf.push_back(static_cast<uint8_t>((v >> (i * 8)) & 0xFF));
    }
}

// Write float as 4 bytes, little-endian (IEEE 754)
inline void write_le_float(std::vector<uint8_t>& buf, float v) {
    uint32_t bits;
    memcpy(&bits, &v, 4);
    write_le_u32(buf, bits);
}

// Write double as 8 bytes, little-endian (IEEE 754)
inline void write_le_double(std::vector<uint8_t>& buf, double v) {
    uint64_t bits;
    memcpy(&bits, &v, 8);
    write_le_u64(buf, bits);
}

// Read uint32_t from little-endian bytes
inline uint32_t read_le_u32(const uint8_t* data, size_t& pos) {
    uint32_t v = static_cast<uint32_t>(data[pos])
               | (static_cast<uint32_t>(data[pos + 1]) << 8)
               | (static_cast<uint32_t>(data[pos + 2]) << 16)
               | (static_cast<uint32_t>(data[pos + 3]) << 24);
    pos += 4;
    return v;
}

// Read uint64_t from little-endian bytes
inline uint64_t read_le_u64(const uint8_t* data, size_t& pos) {
    uint64_t v = 0;
    for (int i = 0; i < 8; i++) {
        v |= static_cast<uint64_t>(data[pos + i]) << (i * 8);
    }
    pos += 8;
    return v;
}

// Read float from little-endian bytes (IEEE 754)
inline float read_le_float(const uint8_t* data, size_t& pos) {
    uint32_t bits = read_le_u32(data, pos);
    // Undo the pos advance from read_le_u32, it already advanced
    float v;
    memcpy(&v, &bits, 4);
    return v;
}

// Read double from little-endian bytes (IEEE 754)
inline double read_le_double(const uint8_t* data, size_t& pos) {
    uint64_t bits = read_le_u64(data, pos);
    double v;
    memcpy(&v, &bits, 8);
    return v;
}

// ============================================================
// Wire format structs (from contracts/rpc-transport.md)
// These will eventually live in autoc/rpc_protocol.h
// ============================================================

// Request: autoc → minisim
struct EvalRequest {
    static constexpr char MAGIC[4] = {'E', 'V', 'A', 'L'};
    static constexpr uint32_t VERSION = 1;

    uint32_t genome_count = 0;
    struct GenomeData {
        uint32_t weight_count = 0;
        std::vector<float> weights;
    };
    std::vector<GenomeData> genomes;

    uint32_t scenario_count = 0;
    struct ScenarioData {
        double wind_speed = 0;
        double wind_direction = 0;
        double turbulence = 0;
        uint64_t path_seed = 0;
    };
    std::vector<ScenarioData> scenarios;

    // Serialize to binary buffer (little-endian wire format)
    std::vector<uint8_t> serialize() const {
        std::vector<uint8_t> buf;

        // Magic (4 bytes, not endian-sensitive)
        buf.insert(buf.end(), MAGIC, MAGIC + 4);
        write_le_u32(buf, VERSION);
        write_le_u32(buf, genome_count);
        for (const auto& g : genomes) {
            write_le_u32(buf, g.weight_count);
            for (float w : g.weights) {
                write_le_float(buf, w);
            }
        }
        write_le_u32(buf, scenario_count);
        for (const auto& s : scenarios) {
            write_le_double(buf, s.wind_speed);
            write_le_double(buf, s.wind_direction);
            write_le_double(buf, s.turbulence);
            write_le_u64(buf, s.path_seed);
        }
        return buf;
    }

    // Deserialize from binary buffer (little-endian wire format)
    static EvalRequest deserialize(const uint8_t* data, size_t len) {
        EvalRequest req;
        size_t pos = 0;

        // Verify magic
        if (len < 4 || memcmp(data + pos, MAGIC, 4) != 0) {
            throw std::runtime_error("Invalid magic bytes in EvalRequest");
        }
        pos += 4;

        uint32_t version = read_le_u32(data, pos);
        if (version != VERSION) {
            throw std::runtime_error("Unsupported EvalRequest version");
        }

        req.genome_count = read_le_u32(data, pos);
        req.genomes.resize(req.genome_count);
        for (auto& g : req.genomes) {
            g.weight_count = read_le_u32(data, pos);
            g.weights.resize(g.weight_count);
            for (uint32_t i = 0; i < g.weight_count; i++) {
                g.weights[i] = read_le_float(data, pos);
            }
        }

        req.scenario_count = read_le_u32(data, pos);
        req.scenarios.resize(req.scenario_count);
        for (auto& s : req.scenarios) {
            s.wind_speed = read_le_double(data, pos);
            s.wind_direction = read_le_double(data, pos);
            s.turbulence = read_le_double(data, pos);
            s.path_seed = read_le_u64(data, pos);
        }

        return req;
    }
};

constexpr char EvalRequest::MAGIC[4];

// Response: minisim → autoc
struct EvalResponse {
    static constexpr char MAGIC[4] = {'R', 'S', 'L', 'T'};
    static constexpr uint32_t VERSION = 1;

    uint32_t genome_count = 0;
    struct GenomeResult {
        uint32_t scenario_count = 0;
        struct ScenarioResult {
            double fitness = 0;
            uint8_t crashed = 0;
            uint32_t timestep_count = 0;
            struct TimestepData {
                uint32_t time_ms = 0;
                float distance = 0;
                float dphi = 0;
                float dtheta = 0;
                float pitch_cmd = 0;
                float roll_cmd = 0;
                float throttle_cmd = 0;
            };
            std::vector<TimestepData> timesteps;
        };
        std::vector<ScenarioResult> scenarios;
    };
    std::vector<GenomeResult> genomes;

    std::vector<uint8_t> serialize() const {
        std::vector<uint8_t> buf;

        buf.insert(buf.end(), MAGIC, MAGIC + 4);
        write_le_u32(buf, VERSION);
        write_le_u32(buf, genome_count);
        for (const auto& g : genomes) {
            write_le_u32(buf, g.scenario_count);
            for (const auto& s : g.scenarios) {
                write_le_double(buf, s.fitness);
                buf.push_back(s.crashed);
                write_le_u32(buf, s.timestep_count);
                for (const auto& t : s.timesteps) {
                    write_le_u32(buf, t.time_ms);
                    write_le_float(buf, t.distance);
                    write_le_float(buf, t.dphi);
                    write_le_float(buf, t.dtheta);
                    write_le_float(buf, t.pitch_cmd);
                    write_le_float(buf, t.roll_cmd);
                    write_le_float(buf, t.throttle_cmd);
                }
            }
        }
        return buf;
    }

    static EvalResponse deserialize(const uint8_t* data, size_t len) {
        EvalResponse resp;
        size_t pos = 0;

        if (len < 4 || memcmp(data + pos, MAGIC, 4) != 0) {
            throw std::runtime_error("Invalid magic bytes in EvalResponse");
        }
        pos += 4;

        uint32_t version = read_le_u32(data, pos);
        if (version != VERSION) {
            throw std::runtime_error("Unsupported EvalResponse version");
        }

        resp.genome_count = read_le_u32(data, pos);
        resp.genomes.resize(resp.genome_count);
        for (auto& g : resp.genomes) {
            g.scenario_count = read_le_u32(data, pos);
            g.scenarios.resize(g.scenario_count);
            for (auto& s : g.scenarios) {
                s.fitness = read_le_double(data, pos);
                s.crashed = data[pos++];
                s.timestep_count = read_le_u32(data, pos);
                s.timesteps.resize(s.timestep_count);
                for (auto& t : s.timesteps) {
                    t.time_ms = read_le_u32(data, pos);
                    t.distance = read_le_float(data, pos);
                    t.dphi = read_le_float(data, pos);
                    t.dtheta = read_le_float(data, pos);
                    t.pitch_cmd = read_le_float(data, pos);
                    t.roll_cmd = read_le_float(data, pos);
                    t.throttle_cmd = read_le_float(data, pos);
                }
            }
        }
        return resp;
    }
};

constexpr char EvalResponse::MAGIC[4];

} // namespace rpc_contract

// ============================================================
// Contract: Wire format is explicitly little-endian
// ============================================================

TEST(ContractRPC, WireFormatIsLittleEndian) {
    // Verify that uint32_t 0x01020304 serializes as bytes [04, 03, 02, 01]
    std::vector<uint8_t> buf;
    rpc_contract::write_le_u32(buf, 0x01020304);
    ASSERT_EQ(buf.size(), 4u);
    EXPECT_EQ(buf[0], 0x04);  // LSB first
    EXPECT_EQ(buf[1], 0x03);
    EXPECT_EQ(buf[2], 0x02);
    EXPECT_EQ(buf[3], 0x01);  // MSB last
}

TEST(ContractRPC, WireFormatU64LittleEndian) {
    std::vector<uint8_t> buf;
    rpc_contract::write_le_u64(buf, 0x0102030405060708ULL);
    ASSERT_EQ(buf.size(), 8u);
    EXPECT_EQ(buf[0], 0x08);  // LSB first
    EXPECT_EQ(buf[7], 0x01);  // MSB last
}

TEST(ContractRPC, WireFormatFloatIEEE754) {
    // float 1.0f = 0x3F800000 in IEEE 754
    // Little-endian bytes: [00, 00, 80, 3F]
    std::vector<uint8_t> buf;
    rpc_contract::write_le_float(buf, 1.0f);
    ASSERT_EQ(buf.size(), 4u);
    EXPECT_EQ(buf[0], 0x00);
    EXPECT_EQ(buf[1], 0x00);
    EXPECT_EQ(buf[2], 0x80);
    EXPECT_EQ(buf[3], 0x3F);
}

TEST(ContractRPC, WireFormatDoubleIEEE754) {
    // double 1.0 = 0x3FF0000000000000 in IEEE 754
    // Little-endian bytes: [00,00,00,00,00,00,F0,3F]
    std::vector<uint8_t> buf;
    rpc_contract::write_le_double(buf, 1.0);
    ASSERT_EQ(buf.size(), 8u);
    EXPECT_EQ(buf[0], 0x00);
    EXPECT_EQ(buf[1], 0x00);
    EXPECT_EQ(buf[2], 0x00);
    EXPECT_EQ(buf[3], 0x00);
    EXPECT_EQ(buf[4], 0x00);
    EXPECT_EQ(buf[5], 0x00);
    EXPECT_EQ(buf[6], 0xF0);
    EXPECT_EQ(buf[7], 0x3F);
}

TEST(ContractRPC, LittleEndianRoundTrip) {
    // Verify read/write round-trips for all types
    std::vector<uint8_t> buf;
    rpc_contract::write_le_u32(buf, 0xDEADBEEF);
    rpc_contract::write_le_u64(buf, 0xCAFEBABE12345678ULL);
    rpc_contract::write_le_float(buf, -3.14f);
    rpc_contract::write_le_double(buf, 2.718281828);

    size_t pos = 0;
    EXPECT_EQ(rpc_contract::read_le_u32(buf.data(), pos), 0xDEADBEEF);
    EXPECT_EQ(rpc_contract::read_le_u64(buf.data(), pos), 0xCAFEBABE12345678ULL);
    EXPECT_FLOAT_EQ(rpc_contract::read_le_float(buf.data(), pos), -3.14f);
    EXPECT_DOUBLE_EQ(rpc_contract::read_le_double(buf.data(), pos), 2.718281828);
}

// ============================================================
// Contract: EvalRequest round-trip serialization
// ============================================================

TEST(ContractRPC, EvalRequestRoundTrip) {
    rpc_contract::EvalRequest req;
    req.genome_count = 2;
    req.genomes.resize(2);

    // Genome 0: 531 weights
    req.genomes[0].weight_count = 531;
    req.genomes[0].weights.resize(531);
    for (int i = 0; i < 531; i++) {
        req.genomes[0].weights[i] = static_cast<float>(i) * 0.001f;
    }

    // Genome 1: 531 weights
    req.genomes[1].weight_count = 531;
    req.genomes[1].weights.resize(531);
    for (int i = 0; i < 531; i++) {
        req.genomes[1].weights[i] = -static_cast<float>(i) * 0.002f;
    }

    req.scenario_count = 3;
    req.scenarios.resize(3);
    req.scenarios[0] = {5.0, 180.0, 0.1, 12345};
    req.scenarios[1] = {10.0, 270.0, 0.3, 67890};
    req.scenarios[2] = {0.0, 0.0, 0.0, 99999};

    // Serialize
    auto buf = req.serialize();
    EXPECT_GT(buf.size(), 0u);

    // Deserialize
    auto req2 = rpc_contract::EvalRequest::deserialize(buf.data(), buf.size());

    // Verify
    EXPECT_EQ(req2.genome_count, 2u);
    EXPECT_EQ(req2.genomes[0].weight_count, 531u);
    EXPECT_EQ(req2.genomes[1].weight_count, 531u);

    for (int i = 0; i < 531; i++) {
        EXPECT_EQ(req2.genomes[0].weights[i], req.genomes[0].weights[i]);
        EXPECT_EQ(req2.genomes[1].weights[i], req.genomes[1].weights[i]);
    }

    EXPECT_EQ(req2.scenario_count, 3u);
    EXPECT_DOUBLE_EQ(req2.scenarios[0].wind_speed, 5.0);
    EXPECT_DOUBLE_EQ(req2.scenarios[1].wind_direction, 270.0);
    EXPECT_DOUBLE_EQ(req2.scenarios[2].turbulence, 0.0);
    EXPECT_EQ(req2.scenarios[0].path_seed, 12345u);
}

// ============================================================
// Contract: EvalResponse round-trip (no per-timestep data)
// ============================================================

TEST(ContractRPC, EvalResponseRoundTripNoTimesteps) {
    rpc_contract::EvalResponse resp;
    resp.genome_count = 1;
    resp.genomes.resize(1);
    resp.genomes[0].scenario_count = 2;
    resp.genomes[0].scenarios.resize(2);

    resp.genomes[0].scenarios[0].fitness = 1234.56;
    resp.genomes[0].scenarios[0].crashed = 0;
    resp.genomes[0].scenarios[0].timestep_count = 0;

    resp.genomes[0].scenarios[1].fitness = 9999.99;
    resp.genomes[0].scenarios[1].crashed = 1;
    resp.genomes[0].scenarios[1].timestep_count = 0;

    auto buf = resp.serialize();
    auto resp2 = rpc_contract::EvalResponse::deserialize(buf.data(), buf.size());

    EXPECT_EQ(resp2.genome_count, 1u);
    EXPECT_EQ(resp2.genomes[0].scenario_count, 2u);
    EXPECT_DOUBLE_EQ(resp2.genomes[0].scenarios[0].fitness, 1234.56);
    EXPECT_EQ(resp2.genomes[0].scenarios[0].crashed, 0);
    EXPECT_EQ(resp2.genomes[0].scenarios[0].timestep_count, 0u);
    EXPECT_DOUBLE_EQ(resp2.genomes[0].scenarios[1].fitness, 9999.99);
    EXPECT_EQ(resp2.genomes[0].scenarios[1].crashed, 1);
}

// ============================================================
// Contract: EvalResponse round-trip (with per-timestep data)
// ============================================================

TEST(ContractRPC, EvalResponseRoundTripWithTimesteps) {
    rpc_contract::EvalResponse resp;
    resp.genome_count = 1;
    resp.genomes.resize(1);
    resp.genomes[0].scenario_count = 1;
    resp.genomes[0].scenarios.resize(1);

    auto& scenario = resp.genomes[0].scenarios[0];
    scenario.fitness = 500.0;
    scenario.crashed = 0;
    scenario.timestep_count = 3;
    scenario.timesteps.resize(3);

    scenario.timesteps[0] = {100, 15.0f, 0.1f, -0.2f, 0.5f, -0.3f, 0.8f};
    scenario.timesteps[1] = {200, 12.0f, 0.05f, -0.1f, 0.4f, -0.2f, 0.7f};
    scenario.timesteps[2] = {300, 8.0f, 0.02f, -0.05f, 0.3f, -0.1f, 0.6f};

    auto buf = resp.serialize();
    auto resp2 = rpc_contract::EvalResponse::deserialize(buf.data(), buf.size());

    auto& s = resp2.genomes[0].scenarios[0];
    EXPECT_DOUBLE_EQ(s.fitness, 500.0);
    EXPECT_EQ(s.timestep_count, 3u);
    EXPECT_EQ(s.timesteps[0].time_ms, 100u);
    EXPECT_FLOAT_EQ(s.timesteps[0].distance, 15.0f);
    EXPECT_FLOAT_EQ(s.timesteps[1].dphi, 0.05f);
    EXPECT_FLOAT_EQ(s.timesteps[2].pitch_cmd, 0.3f);
    EXPECT_FLOAT_EQ(s.timesteps[2].throttle_cmd, 0.6f);
}

// ============================================================
// Contract: Magic byte validation
// ============================================================

TEST(ContractRPC, EvalRequestInvalidMagic) {
    std::vector<uint8_t> bad_data = {'B', 'A', 'D', '!', 0, 0, 0, 1};
    EXPECT_THROW(
        rpc_contract::EvalRequest::deserialize(bad_data.data(), bad_data.size()),
        std::runtime_error
    );
}

TEST(ContractRPC, EvalResponseInvalidMagic) {
    std::vector<uint8_t> bad_data = {'B', 'A', 'D', '!', 0, 0, 0, 1};
    EXPECT_THROW(
        rpc_contract::EvalResponse::deserialize(bad_data.data(), bad_data.size()),
        std::runtime_error
    );
}

// ============================================================
// Contract: Empty genomes/scenarios round-trip
// ============================================================

TEST(ContractRPC, EmptyRequest) {
    rpc_contract::EvalRequest req;
    req.genome_count = 0;
    req.scenario_count = 0;

    auto buf = req.serialize();
    auto req2 = rpc_contract::EvalRequest::deserialize(buf.data(), buf.size());

    EXPECT_EQ(req2.genome_count, 0u);
    EXPECT_EQ(req2.scenario_count, 0u);
}

// ============================================================
// Contract: Negative float values survive serialization
// (Verifies sign bit handling in IEEE 754 LE encoding)
// ============================================================

TEST(ContractRPC, NegativeValuesRoundTrip) {
    rpc_contract::EvalResponse resp;
    resp.genome_count = 1;
    resp.genomes.resize(1);
    resp.genomes[0].scenario_count = 1;
    resp.genomes[0].scenarios.resize(1);

    auto& s = resp.genomes[0].scenarios[0];
    s.fitness = -12345.6789;
    s.crashed = 0;
    s.timestep_count = 1;
    s.timesteps.resize(1);
    s.timesteps[0] = {0, -1.0f, -0.5f, -0.25f, -0.75f, -0.125f, -1.0f};

    auto buf = resp.serialize();
    auto resp2 = rpc_contract::EvalResponse::deserialize(buf.data(), buf.size());

    EXPECT_DOUBLE_EQ(resp2.genomes[0].scenarios[0].fitness, -12345.6789);
    EXPECT_FLOAT_EQ(resp2.genomes[0].scenarios[0].timesteps[0].distance, -1.0f);
    EXPECT_FLOAT_EQ(resp2.genomes[0].scenarios[0].timesteps[0].dphi, -0.5f);
    EXPECT_FLOAT_EQ(resp2.genomes[0].scenarios[0].timesteps[0].roll_cmd, -0.125f);
}

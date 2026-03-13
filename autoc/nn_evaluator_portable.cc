#include "nn_evaluator_portable.h"
#include <cmath>
#include <array>
#include <algorithm>

#ifndef GP_TEST
#include "gp_evaluator_portable.h"
#endif

#if defined(GP_BUILD) && !defined(GP_TEST)
#include "gp_math_utils.h"
#endif

// ============================================================
// T036: Weight count
// ============================================================

int nn_weight_count(const std::vector<int>& topology) {
    int count = 0;
    for (size_t i = 0; i + 1 < topology.size(); i++) {
        count += topology[i] * topology[i + 1] + topology[i + 1];
    }
    return count;
}

// ============================================================
// T035: fast_tanh LUT — 512 entries, domain [-5, 5]
// ============================================================

namespace {
constexpr int TANH_LUT_SIZE = 512;
constexpr gp_scalar TANH_DOMAIN = static_cast<gp_scalar>(5.0f);
std::array<gp_scalar, TANH_LUT_SIZE + 1> TANH_LUT{};
bool TANH_LUT_INIT = false;

inline void initTanhLut() {
    if (TANH_LUT_INIT) return;
    for (int i = 0; i <= TANH_LUT_SIZE; ++i) {
        gp_scalar x = -TANH_DOMAIN + static_cast<gp_scalar>(2.0f) * TANH_DOMAIN
                       * static_cast<gp_scalar>(i) / static_cast<gp_scalar>(TANH_LUT_SIZE);
        TANH_LUT[i] = static_cast<gp_scalar>(std::tanh(x));
    }
    TANH_LUT_INIT = true;
}
} // namespace

gp_scalar fast_tanh(gp_scalar x) {
    initTanhLut();

    // Clamp to domain
    if (x <= -TANH_DOMAIN) return static_cast<gp_scalar>(-1.0f);
    if (x >= TANH_DOMAIN) return static_cast<gp_scalar>(1.0f);

    // Map x from [-5, 5] to [0, TANH_LUT_SIZE]
    gp_scalar scaled = (x + TANH_DOMAIN) / (static_cast<gp_scalar>(2.0f) * TANH_DOMAIN)
                       * static_cast<gp_scalar>(TANH_LUT_SIZE);
    int i0 = static_cast<int>(scaled);
    if (i0 < 0) i0 = 0;
    if (i0 >= TANH_LUT_SIZE) i0 = TANH_LUT_SIZE - 1;
    gp_scalar frac = scaled - static_cast<gp_scalar>(i0);
    return TANH_LUT[i0] + (TANH_LUT[i0 + 1] - TANH_LUT[i0]) * frac;
}

// ============================================================
// T034: Forward pass — feedforward with tanh activation
// Weight layout: row-major, layer-sequential [W1, B1, W2, B2, ...]
// For layer l: W has shape [out_size x in_size], B has shape [out_size]
// output[j] = tanh(sum_i(W[j*in_size + i] * input[i]) + B[j])
// ============================================================

void nn_forward(const float* weights, const std::vector<int>& topology,
                const float* inputs, float* outputs) {
    if (topology.size() < 2) return;

    // Double-buffer: alternate between buf_a and buf_b for layer I/O
    int max_layer = *std::max_element(topology.begin(), topology.end());
    // Use stack allocation for typical sizes, heap for unusual
    float buf_a_stack[64], buf_b_stack[64];
    float* buf_a = (max_layer <= 64) ? buf_a_stack : new float[max_layer];
    float* buf_b = (max_layer <= 64) ? buf_b_stack : new float[max_layer];

    // Copy inputs to buf_a
    for (int i = 0; i < topology[0]; i++) {
        buf_a[i] = inputs[i];
    }

    const float* w_ptr = weights;
    float* current_in = buf_a;
    float* current_out = buf_b;

    for (size_t layer = 0; layer + 1 < topology.size(); layer++) {
        int in_size = topology[layer];
        int out_size = topology[layer + 1];

        // Compute: current_out[j] = tanh(sum(W[j*in_size+i] * current_in[i]) + B[j])
        const float* W = w_ptr;
        const float* B = w_ptr + in_size * out_size;

        for (int j = 0; j < out_size; j++) {
            float sum = B[j];
            for (int i = 0; i < in_size; i++) {
                sum += W[j * in_size + i] * current_in[i];
            }
            current_out[j] = static_cast<float>(fast_tanh(static_cast<gp_scalar>(sum)));
        }

        w_ptr += in_size * out_size + out_size;

        // Swap buffers
        float* tmp = current_in;
        current_in = current_out;
        current_out = tmp;
    }

    // Copy final output (current_in points to last layer's output after swap)
    int out_size = topology.back();
    for (int i = 0; i < out_size; i++) {
        outputs[i] = current_in[i];
    }

    // Clean up heap allocations if used
    if (max_layer > 64) {
        delete[] buf_a;
        delete[] buf_b;
    }
}

// ============================================================
// T037: Xavier/Glorot initialization using GPrand()
// ============================================================

// Local Box-Muller Gaussian for builds without GPrand (test/embedded)
#if !defined(GP_BUILD) || defined(GP_TEST)
namespace {
static long nn_rand_state = 42;
inline long nn_local_rand() {
    nn_rand_state = (nn_rand_state * 1103515245L + 12345L) & 0x7FFFFFFF;
    return nn_rand_state;
}
inline double nn_local_gaussian(double sigma) {
    constexpr double RNMX = 2147483646.0;
    double u1 = (nn_local_rand() / RNMX) * 0.999 + 0.001;
    double u2 = nn_local_rand() / RNMX;
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    return z * sigma;
}
} // namespace
#endif

void nn_xavier_init(NNGenome& genome) {
    int total = nn_weight_count(genome.topology);
    genome.weights.resize(total);

    int idx = 0;
    for (size_t layer = 0; layer + 1 < genome.topology.size(); layer++) {
        int fan_in = genome.topology[layer];
        int fan_out = genome.topology[layer + 1];
        double stddev = std::sqrt(1.0 / fan_in);

        // Initialize weights with Gaussian(0, stddev)
        int num_weights = fan_in * fan_out;
        for (int i = 0; i < num_weights; i++) {
#if defined(GP_BUILD) && !defined(GP_TEST)
            genome.weights[idx++] = static_cast<float>(GPrandGaussian(stddev));
#else
            genome.weights[idx++] = static_cast<float>(nn_local_gaussian(stddev));
#endif
        }

        // Initialize biases to zero
        for (int i = 0; i < fan_out; i++) {
            genome.weights[idx++] = 0.0f;
        }
    }
}

// ============================================================
// T038: Gather 14 sensor inputs with normalization
// ============================================================

void nn_gather_inputs(PathProvider& pathProvider, AircraftState& aircraftState,
                      float* inputs) {
#ifdef GP_TEST
    // In test builds, zero-fill (sensor functions not available)
    for (int i = 0; i < 14; i++) inputs[i] = 0.0f;
#else
    // Desktop (GP_BUILD) and embedded builds both have executeGet* available
    // 0-1: Navigation angles
    inputs[0] = static_cast<float>(executeGetDPhi(pathProvider, aircraftState, 0.0f) / NORM_ANGLE);
    inputs[1] = static_cast<float>(executeGetDTheta(pathProvider, aircraftState, 0.0f) / NORM_ANGLE);

    // 2: Distance to target
    inputs[2] = static_cast<float>(executeGetDist(pathProvider, aircraftState) / NORM_DIST);

    // 3-5: Rate sensors
    inputs[3] = static_cast<float>(executeGetDistRate(aircraftState) / NORM_RATE);
    inputs[4] = static_cast<float>(executeGetDPhiRate(aircraftState) / NORM_RATE);
    inputs[5] = static_cast<float>(executeGetDThetaRate(aircraftState) / NORM_RATE);

    // 6-7: Attitude (roll/pitch from quaternion Euler extraction)
    {
        gp_vec3 euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
        inputs[6] = static_cast<float>(euler[2] / NORM_ANGLE);  // roll
        inputs[7] = static_cast<float>(euler[1] / NORM_ANGLE);  // pitch
    }

    // 8: Velocity
    inputs[8] = static_cast<float>(aircraftState.getRelVel() / NORM_VEL);

    // 9-10: Alpha/Beta — compute inline (same as evaluateGPOperator GETALPHA/GETBETA cases)
    {
        gp_vec3 velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
        gp_scalar alpha = std::atan2(-velocity_body.z(), velocity_body.x());
        gp_scalar beta = std::atan2(velocity_body.y(), velocity_body.x());
        inputs[9] = static_cast<float>(alpha / NORM_ANGLE);
        inputs[10] = static_cast<float>(beta / NORM_ANGLE);
    }

    // 11-13: Current control commands (already [-1, 1], no normalization)
    inputs[11] = static_cast<float>(aircraftState.getRollCommand());
    inputs[12] = static_cast<float>(aircraftState.getPitchCommand());
    inputs[13] = static_cast<float>(aircraftState.getThrottleCommand());
#endif
}

// ============================================================
// T039: NNControllerBackend
// ============================================================

#ifdef GP_BUILD
NNControllerBackend::NNControllerBackend(const NNGenome& genome)
    : genome_(genome) {}

void NNControllerBackend::evaluate(AircraftState& aircraftState, PathProvider& pathProvider) {
    float inputs[14];
    nn_gather_inputs(pathProvider, aircraftState, inputs);

    float outputs[3];
    nn_forward(genome_.weights.data(), genome_.topology, inputs, outputs);

    // Set control commands: pitch, roll, throttle (already in [-1, 1] via tanh)
    aircraftState.setPitchCommand(static_cast<gp_scalar>(outputs[0]));
    aircraftState.setRollCommand(static_cast<gp_scalar>(outputs[1]));
    aircraftState.setThrottleCommand(static_cast<gp_scalar>(outputs[2]));
}
#endif

// ============================================================
// Test helpers
// ============================================================

#ifdef GP_TEST
gp_scalar testFastTanh(gp_scalar x) {
    return fast_tanh(x);
}
#endif

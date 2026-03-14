#ifndef NN_EVALUATOR_PORTABLE_H
#define NN_EVALUATOR_PORTABLE_H

#include <vector>
#include <cstdint>
#include "gp_types.h"
#include "aircraft_state.h"

// Input normalization constants (applied before forward pass)
// See research.md R9 for rationale
constexpr gp_scalar NORM_ANGLE = static_cast<gp_scalar>(M_PI);  // angles: dPhi, dTheta, roll, pitch, alpha, beta
constexpr gp_scalar NORM_DIST  = static_cast<gp_scalar>(50.0f); // distance (meters)
constexpr gp_scalar NORM_VEL   = static_cast<gp_scalar>(16.0f); // velocity (nominal rabbit speed, m/s)
constexpr gp_scalar NORM_RATE  = static_cast<gp_scalar>(10.0f); // rate sensors (already clamped [-10, 10])

// Neural network genome — the fundamental unit of neuroevolution
struct NNGenome {
    std::vector<float> weights;       // All weights + biases, layer-major order
    std::vector<int> topology;        // Layer sizes, e.g., {22, 16, 8, 3}
    double fitness;                   // Aggregated fitness from evaluation
    uint32_t generation;              // Generation when created
    float mutation_sigma;             // Per-individual mutation step size (self-adaptive)

    NNGenome() : fitness(0.0), generation(0), mutation_sigma(0.1f) {}
};

// Compute total weight+bias count from topology
int nn_weight_count(const std::vector<int>& topology);

// Forward pass: feedforward with tanh activation
// weights layout: row-major, layer-sequential [W1, B1, W2, B2, ...]
void nn_forward(const float* weights, const std::vector<int>& topology,
                const float* inputs, float* outputs);

// Fast tanh via 512-entry LUT with linear interpolation
gp_scalar fast_tanh(gp_scalar x);

// Xavier/Glorot weight initialization
void nn_xavier_init(NNGenome& genome);

// Gather 22 sensor inputs from aircraft state, apply normalization
void nn_gather_inputs(PathProvider& pathProvider, AircraftState& aircraftState,
                      float* inputs);

#ifdef GP_BUILD
#include "eval_backend.h"

// NN controller backend — plugs into unified eval pipeline
class NNControllerBackend : public ControllerBackend {
public:
    explicit NNControllerBackend(const NNGenome& genome);

    void evaluate(AircraftState& aircraftState, PathProvider& pathProvider) override;
    const char* getName() const override { return "NeuralNet"; }

private:
    const NNGenome& genome_;
};
#endif

// Expose LUT functions for testing
#ifdef GP_TEST
gp_scalar testFastTanh(gp_scalar x);
#endif

#endif

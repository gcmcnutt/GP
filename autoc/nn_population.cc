#include "nn_population.h"
#include <cmath>
#include <algorithm>
#include <limits>

// Use local RNG for test/embedded builds, GPrand for full builds
#if defined(GP_BUILD) && !defined(GP_TEST)
#include "gp_math_utils.h"
#else
namespace {
static long nn_pop_rand_state = 54321;
inline long nn_pop_rand() {
    nn_pop_rand_state = (nn_pop_rand_state * 1103515245L + 12345L) & 0x7FFFFFFF;
    return nn_pop_rand_state;
}
inline double nn_pop_rand_double() {
    return nn_pop_rand() / 2147483646.0;
}
inline double nn_pop_rand_gaussian(double sigma) {
    double u1 = nn_pop_rand_double() * 0.999 + 0.001;
    double u2 = nn_pop_rand_double();
    double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    return z * sigma;
}
inline int nn_pop_rand_int(int max_exclusive) {
    return static_cast<int>(nn_pop_rand() % max_exclusive);
}
} // namespace
#endif

// Portable wrappers
namespace {
inline double randDouble() {
#if defined(GP_BUILD) && !defined(GP_TEST)
    return GPrandDouble();
#else
    return nn_pop_rand_double();
#endif
}

inline double randGaussian(double sigma) {
#if defined(GP_BUILD) && !defined(GP_TEST)
    return GPrandGaussian(sigma);
#else
    return nn_pop_rand_gaussian(sigma);
#endif
}

inline int randInt(int max_exclusive) {
#if defined(GP_BUILD) && !defined(GP_TEST)
    return static_cast<int>(GPrand() % max_exclusive);
#else
    return nn_pop_rand_int(max_exclusive);
#endif
}
} // namespace

// ============================================================
// T064: Population initialization
// ============================================================

void nn_init_population(NNPopulation& pop, const std::vector<int>& topology, int size) {
    pop.topology = topology;
    pop.population_size = size;
    pop.generation = 0;
    pop.best_fitness = std::numeric_limits<double>::max();
    pop.best_index = 0;

    pop.individuals.resize(size);
    for (int i = 0; i < size; i++) {
        pop.individuals[i].topology = topology;
        pop.individuals[i].generation = 0;
        pop.individuals[i].mutation_sigma = 0.1f;
        pop.individuals[i].fitness = 0.0;
        nn_xavier_init(pop.individuals[i]);
    }
}

// ============================================================
// T061: Arithmetic crossover (BLX-alpha blend)
// ============================================================

NNGenome nn_arithmetic_crossover(const NNGenome& parent1, const NNGenome& parent2, float alpha) {
    NNGenome child;
    child.topology = parent1.topology;
    child.generation = parent1.generation;
    child.mutation_sigma = (parent1.mutation_sigma + parent2.mutation_sigma) * 0.5f;
    child.fitness = 0.0;

    size_t n = parent1.weights.size();
    child.weights.resize(n);

    if (alpha < 0) {
        // Random alpha per weight (BLX-α with α=-1 means uniform random blend)
        for (size_t i = 0; i < n; i++) {
            float a = static_cast<float>(randDouble());
            child.weights[i] = (1.0f - a) * parent1.weights[i] + a * parent2.weights[i];
        }
    } else {
        // Fixed alpha blend: child = (1-alpha)*p1 + alpha*p2
        for (size_t i = 0; i < n; i++) {
            child.weights[i] = (1.0f - alpha) * parent1.weights[i] + alpha * parent2.weights[i];
        }
    }

    return child;
}

// ============================================================
// T062: Gaussian mutation with self-adaptive sigma
// ============================================================

void nn_gaussian_mutation(NNGenome& genome) {
    // Self-adapt sigma: sigma' = sigma * exp(tau * N(0,1))
    // tau = 1 / sqrt(2 * n) per Schwefel (1981)
    double n = static_cast<double>(genome.weights.size());
    double tau = 1.0 / std::sqrt(2.0 * n);
    double sigma_perturbation = std::exp(tau * randGaussian(1.0));
    genome.mutation_sigma = static_cast<float>(genome.mutation_sigma * sigma_perturbation);

    // Clamp sigma to reasonable range
    constexpr float SIGMA_MIN = 1e-6f;
    constexpr float SIGMA_MAX = 5.0f;
    if (genome.mutation_sigma < SIGMA_MIN) genome.mutation_sigma = SIGMA_MIN;
    if (genome.mutation_sigma > SIGMA_MAX) genome.mutation_sigma = SIGMA_MAX;

    // Mutate weights
    for (size_t i = 0; i < genome.weights.size(); i++) {
        genome.weights[i] += static_cast<float>(randGaussian(genome.mutation_sigma));
    }
}

// ============================================================
// T063: Tournament selection (minimization — lower fitness is better)
// ============================================================

int nn_tournament_select(const NNPopulation& pop, int tournament_size) {
    int best_idx = randInt(pop.population_size);
    double best_fitness = pop.individuals[best_idx].fitness;

    for (int i = 1; i < tournament_size; i++) {
        int candidate = randInt(pop.population_size);
        if (pop.individuals[candidate].fitness < best_fitness) {
            best_idx = candidate;
            best_fitness = pop.individuals[candidate].fitness;
        }
    }

    return best_idx;
}

// ============================================================
// T065: Evolve one generation
// ============================================================

void nn_evolve_generation(NNPopulation& pop, const NNEvolveParams& params) {
    int n = pop.population_size;
    if (n < 2) return;

    // Find elites (best fitness = lowest)
    // Sort indices by fitness to pick top elitism_count
    std::vector<int> sorted_idx(n);
    for (int i = 0; i < n; i++) sorted_idx[i] = i;
    std::sort(sorted_idx.begin(), sorted_idx.end(), [&](int a, int b) {
        return pop.individuals[a].fitness < pop.individuals[b].fitness;
    });

    int elitism_count = std::max(0, std::min(params.elitism_count, n - 1));

    // Create next generation
    std::vector<NNGenome> next_gen;
    next_gen.reserve(n);

    // Elitism: keep top individuals unchanged (GP & NN: AddBestToNewPopulation)
    for (int i = 0; i < elitism_count; i++) {
        next_gen.push_back(pop.individuals[sorted_idx[i]]);
    }

    // Probability buckets (GP & NN: CrossoverProbability, CreationProbability, SwapMutationProbability)
    // crossover_prob: crossover + mutation
    // mutation_prob: mutation only (no crossover)
    // creation_prob: fresh random individual
    // remainder: crossover only (no mutation)
    double total = params.crossover_prob + params.creation_prob + params.mutation_prob;
    double cross_thresh = params.crossover_prob / total;
    double create_thresh = cross_thresh + params.creation_prob / total;
    double mutate_thresh = create_thresh + params.mutation_prob / total;

    int tournament_size = std::max(2, params.tournament_size); // GP & NN: TournamentSize

    // Fill remaining slots
    for (int i = elitism_count; i < n; i++) {
        double r = randDouble();

        if (r < cross_thresh) {
            // Crossover + mutation
            int p1 = nn_tournament_select(pop, tournament_size);
            int p2 = nn_tournament_select(pop, tournament_size);
            NNGenome child = nn_arithmetic_crossover(
                pop.individuals[p1], pop.individuals[p2],
                params.crossover_alpha);  // NN only: NNCrossoverAlpha
            child.generation = pop.generation + 1;
            nn_gaussian_mutation(child);
            next_gen.push_back(std::move(child));
        } else if (r < create_thresh) {
            // Fresh random individual (CreationProbability)
            NNGenome fresh;
            fresh.topology = pop.topology;
            fresh.generation = pop.generation + 1;
            fresh.mutation_sigma = pop.individuals[0].mutation_sigma;
            fresh.fitness = 0.0;
            nn_xavier_init(fresh);
            next_gen.push_back(std::move(fresh));
        } else if (r < mutate_thresh) {
            // Mutation only, no crossover (SwapMutationProbability)
            int p = nn_tournament_select(pop, tournament_size);
            NNGenome child = pop.individuals[p]; // copy
            child.generation = pop.generation + 1;
            nn_gaussian_mutation(child);
            next_gen.push_back(std::move(child));
        } else {
            // Crossover only, no mutation (remainder)
            int p1 = nn_tournament_select(pop, tournament_size);
            int p2 = nn_tournament_select(pop, tournament_size);
            NNGenome child = nn_arithmetic_crossover(
                pop.individuals[p1], pop.individuals[p2],
                params.crossover_alpha);
            child.generation = pop.generation + 1;
            next_gen.push_back(std::move(child));
        }
    }

    pop.individuals = std::move(next_gen);
    pop.generation++;
    pop.best_fitness = pop.individuals[0].fitness;
    pop.best_index = 0;
}

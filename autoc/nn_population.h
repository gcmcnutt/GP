#ifndef NN_POPULATION_H
#define NN_POPULATION_H

#include <vector>
#include <cstdint>
#include "nn_evaluator_portable.h"

// Collection of NNGenome individuals with shared topology
struct NNPopulation {
    std::vector<NNGenome> individuals;
    std::vector<int> topology;        // Shared topology (all same shape)
    int population_size;
    uint32_t generation;
    double best_fitness;
    int best_index;

    NNPopulation() : population_size(0), generation(0), best_fitness(0.0), best_index(0) {}
};

// Initialize population with Xavier-initialized individuals
void nn_init_population(NNPopulation& pop, const std::vector<int>& topology, int size);

// Arithmetic crossover (BLX-alpha blend)
NNGenome nn_arithmetic_crossover(const NNGenome& parent1, const NNGenome& parent2, float alpha);

// Gaussian mutation with self-adaptive sigma
void nn_gaussian_mutation(NNGenome& genome);

// Tournament selection
int nn_tournament_select(const NNPopulation& pop, int tournament_size);

// Evolve one generation
void nn_evolve_generation(NNPopulation& pop);

#endif

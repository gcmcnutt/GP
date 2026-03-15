// Contract test: NN evolution loop
// Defines the surviving behavior: init → evaluate → select → reproduce → verify
// This contract MUST hold through GP removal and optimizer changes.

#include <gtest/gtest.h>
#include "../nn_population.h"
#include "../nn_topology.h"
#include <cmath>
#include <algorithm>
#include <limits>

// ============================================================
// Contract: Population initializes with correct structure
// ============================================================

TEST(ContractEvolution, PopulationInit) {
    NNPopulation pop;
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    nn_init_population(pop, topology, 20);

    EXPECT_EQ(pop.population_size, 20);
    EXPECT_EQ(static_cast<int>(pop.individuals.size()), 20);
    EXPECT_EQ(pop.generation, 0u);
    EXPECT_EQ(pop.topology, topology);

    for (const auto& ind : pop.individuals) {
        EXPECT_EQ(static_cast<int>(ind.weights.size()), NN_WEIGHT_COUNT);
        EXPECT_EQ(ind.topology, topology);
        EXPECT_EQ(ind.generation, 0u);
        EXPECT_GT(ind.mutation_sigma, 0.0f);
    }
}

// ============================================================
// Contract: Single generation evolution cycle
// init → assign fitness → evolve → verify structure preserved
// ============================================================

TEST(ContractEvolution, SingleGenerationCycle) {
    NNPopulation pop;
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    nn_init_population(pop, topology, 20);

    // Simulate fitness evaluation: assign random-ish fitness values
    // Lower fitness = better (minimization)
    for (int i = 0; i < pop.population_size; i++) {
        pop.individuals[i].fitness = static_cast<double>(i * 1000 + 500);
    }

    // Track best fitness before evolution
    double best_before = std::numeric_limits<double>::max();
    for (const auto& ind : pop.individuals) {
        if (ind.fitness < best_before) best_before = ind.fitness;
    }

    // Evolution parameters matching typical config
    NNEvolveParams params;
    params.tournament_size = 3;
    params.crossover_prob = 90.0;
    params.creation_prob = 5.0;
    params.mutation_prob = 5.0;
    params.crossover_alpha = -1.0f;
    params.elitism_count = 1;

    nn_evolve_generation(pop, params);

    // Verify structure after evolution
    EXPECT_EQ(pop.population_size, 20);
    EXPECT_EQ(static_cast<int>(pop.individuals.size()), 20);
    EXPECT_EQ(pop.generation, 1u);

    // All individuals should have valid weights
    for (const auto& ind : pop.individuals) {
        EXPECT_EQ(static_cast<int>(ind.weights.size()), NN_WEIGHT_COUNT);
        for (float w : ind.weights) {
            EXPECT_TRUE(std::isfinite(w)) << "Non-finite weight after evolution";
        }
        EXPECT_GT(ind.mutation_sigma, 0.0f);
    }
}

// ============================================================
// Contract: Elitism preserves best individual
// ============================================================

TEST(ContractEvolution, ElitismPreservesBest) {
    NNPopulation pop;
    std::vector<int> topology(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
    nn_init_population(pop, topology, 20);

    // Set individual 5 as clearly the best
    for (int i = 0; i < pop.population_size; i++) {
        pop.individuals[i].fitness = 1000.0;
    }
    pop.individuals[5].fitness = 1.0;  // Much better
    std::vector<float> best_weights = pop.individuals[5].weights;

    NNEvolveParams params;
    params.tournament_size = 3;
    params.crossover_prob = 90.0;
    params.creation_prob = 5.0;
    params.mutation_prob = 5.0;
    params.crossover_alpha = -1.0f;
    params.elitism_count = 1;

    nn_evolve_generation(pop, params);

    // The best individual's weights should be preserved (elitism)
    bool found_elite = false;
    for (const auto& ind : pop.individuals) {
        if (ind.weights == best_weights) {
            found_elite = true;
            break;
        }
    }
    EXPECT_TRUE(found_elite) << "Elite individual not preserved after evolution";
}

// ============================================================
// Contract: Multiple generations produce improving fitness
// (with a simple synthetic fitness function)
// ============================================================

TEST(ContractEvolution, MultiGenerationImprovement) {
    NNPopulation pop;
    // Use small topology for speed
    std::vector<int> topology = {4, 3, 2};
    nn_init_population(pop, topology, 30);

    NNEvolveParams params;
    params.tournament_size = 3;
    params.crossover_prob = 80.0;
    params.creation_prob = 5.0;
    params.mutation_prob = 15.0;
    params.crossover_alpha = -1.0f;
    params.elitism_count = 2;

    // Simple fitness: minimize sum of squared weights (trivial optimization)
    auto evaluate = [](NNPopulation& p) {
        for (auto& ind : p.individuals) {
            double sum_sq = 0;
            for (float w : ind.weights) {
                sum_sq += w * w;
            }
            ind.fitness = sum_sq;
        }
    };

    double initial_best = std::numeric_limits<double>::max();
    evaluate(pop);
    for (const auto& ind : pop.individuals) {
        if (ind.fitness < initial_best) initial_best = ind.fitness;
    }

    // Run 50 generations
    for (int gen = 0; gen < 50; gen++) {
        nn_evolve_generation(pop, params);
        evaluate(pop);
    }

    double final_best = std::numeric_limits<double>::max();
    for (const auto& ind : pop.individuals) {
        if (ind.fitness < final_best) final_best = ind.fitness;
    }

    // Fitness should improve (decrease) over 50 generations on this trivial problem
    EXPECT_LT(final_best, initial_best)
        << "Fitness did not improve over 50 generations: "
        << final_best << " >= " << initial_best;
}

// ============================================================
// Contract: Generation counter increments correctly
// ============================================================

TEST(ContractEvolution, GenerationCounterIncrements) {
    NNPopulation pop;
    std::vector<int> topology = {4, 2};
    nn_init_population(pop, topology, 10);

    NNEvolveParams params;
    params.tournament_size = 2;
    params.crossover_prob = 80.0;
    params.creation_prob = 10.0;
    params.mutation_prob = 10.0;
    params.crossover_alpha = 0.5f;
    params.elitism_count = 1;

    // Assign dummy fitness
    for (auto& ind : pop.individuals) {
        ind.fitness = 100.0;
    }

    EXPECT_EQ(pop.generation, 0u);
    nn_evolve_generation(pop, params);
    EXPECT_EQ(pop.generation, 1u);
    nn_evolve_generation(pop, params);
    EXPECT_EQ(pop.generation, 2u);
}

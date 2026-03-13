#include <gtest/gtest.h>
#include "../nn_population.h"
#include <cmath>
#include <map>
#include <algorithm>

// ============================================================
// T058: Population initialization tests
// ============================================================

TEST(NNPopulation, DefaultConstruction) {
    NNPopulation pop;
    EXPECT_EQ(pop.population_size, 0);
    EXPECT_EQ(pop.generation, 0u);
}

TEST(NNPopulation, InitPopulationCreatesCorrectSize) {
    NNPopulation pop;
    std::vector<int> topology = {14, 16, 8, 3};
    nn_init_population(pop, topology, 50);

    EXPECT_EQ(pop.population_size, 50);
    EXPECT_EQ(static_cast<int>(pop.individuals.size()), 50);
    EXPECT_EQ(pop.topology, topology);
    EXPECT_EQ(pop.generation, 0u);
}

TEST(NNPopulation, InitPopulationCorrectWeightCount) {
    NNPopulation pop;
    std::vector<int> topology = {14, 16, 8, 3};
    nn_init_population(pop, topology, 10);

    int expected_weights = nn_weight_count(topology);
    for (const auto& ind : pop.individuals) {
        EXPECT_EQ(static_cast<int>(ind.weights.size()), expected_weights);
        EXPECT_EQ(ind.topology, topology);
    }
}

TEST(NNPopulation, InitPopulationAllFinite) {
    NNPopulation pop;
    std::vector<int> topology = {14, 16, 8, 3};
    nn_init_population(pop, topology, 20);

    for (const auto& ind : pop.individuals) {
        for (float w : ind.weights) {
            EXPECT_TRUE(std::isfinite(w)) << "Non-finite weight found";
        }
    }
}

// ============================================================
// T055: Arithmetic crossover tests
// ============================================================

TEST(NNCrossover, ArithmeticBlend) {
    NNGenome parent1, parent2;
    parent1.topology = {2, 1};
    parent2.topology = {2, 1};
    parent1.weights = {1.0f, 0.0f, 0.0f};
    parent2.weights = {0.0f, 1.0f, 1.0f};

    // alpha = 0.5 -> uniform blend: child = 0.5 * p1 + 0.5 * p2
    NNGenome child = nn_arithmetic_crossover(parent1, parent2, 0.5f);

    ASSERT_EQ(child.weights.size(), 3u);
    EXPECT_NEAR(child.weights[0], 0.5f, 1e-5);
    EXPECT_NEAR(child.weights[1], 0.5f, 1e-5);
    EXPECT_NEAR(child.weights[2], 0.5f, 1e-5);
}

TEST(NNCrossover, ArithmeticBlendAlpha0) {
    NNGenome parent1, parent2;
    parent1.topology = {2, 1};
    parent2.topology = {2, 1};
    parent1.weights = {1.0f, 2.0f, 3.0f};
    parent2.weights = {4.0f, 5.0f, 6.0f};

    // alpha = 0 -> child = parent1 (no blending)
    NNGenome child = nn_arithmetic_crossover(parent1, parent2, 0.0f);

    EXPECT_NEAR(child.weights[0], 1.0f, 1e-5);
    EXPECT_NEAR(child.weights[1], 2.0f, 1e-5);
    EXPECT_NEAR(child.weights[2], 3.0f, 1e-5);
}

TEST(NNCrossover, ArithmeticBlendAlpha1) {
    NNGenome parent1, parent2;
    parent1.topology = {2, 1};
    parent2.topology = {2, 1};
    parent1.weights = {1.0f, 2.0f, 3.0f};
    parent2.weights = {4.0f, 5.0f, 6.0f};

    // alpha = 1 -> child = parent2
    NNGenome child = nn_arithmetic_crossover(parent1, parent2, 1.0f);

    EXPECT_NEAR(child.weights[0], 4.0f, 1e-5);
    EXPECT_NEAR(child.weights[1], 5.0f, 1e-5);
    EXPECT_NEAR(child.weights[2], 6.0f, 1e-5);
}

TEST(NNCrossover, PreservesTopology) {
    NNGenome parent1, parent2;
    parent1.topology = {14, 16, 8, 3};
    parent2.topology = {14, 16, 8, 3};
    int wc = nn_weight_count(parent1.topology);
    parent1.weights.resize(wc, 0.0f);
    parent2.weights.resize(wc, 1.0f);

    NNGenome child = nn_arithmetic_crossover(parent1, parent2, 0.5f);
    EXPECT_EQ(child.topology, parent1.topology);
    EXPECT_EQ(static_cast<int>(child.weights.size()), wc);
}

// ============================================================
// T056: Gaussian mutation tests
// ============================================================

TEST(NNMutation, ChangesWeights) {
    NNGenome genome;
    genome.topology = {14, 16, 8, 3};
    genome.weights.resize(nn_weight_count(genome.topology), 0.0f);
    genome.mutation_sigma = 0.1f;

    // Save original weights
    std::vector<float> original = genome.weights;

    nn_gaussian_mutation(genome);

    // At least some weights should have changed
    bool any_changed = false;
    for (size_t i = 0; i < genome.weights.size(); i++) {
        if (genome.weights[i] != original[i]) {
            any_changed = true;
            break;
        }
    }
    EXPECT_TRUE(any_changed) << "Mutation didn't change any weights";
}

TEST(NNMutation, MeanChangeNearZero) {
    NNGenome genome;
    genome.topology = {14, 16, 8, 3};
    genome.weights.resize(nn_weight_count(genome.topology), 0.0f);
    genome.mutation_sigma = 0.1f;

    nn_gaussian_mutation(genome);

    // Mean of mutated weights (from zero) should be near zero
    double sum = 0;
    for (float w : genome.weights) {
        sum += w;
    }
    double mean = sum / genome.weights.size();
    EXPECT_NEAR(mean, 0.0, 0.05) << "Mean mutation offset too large";
}

// ============================================================
// T057: Self-adaptive sigma test
// ============================================================

TEST(NNMutation, SigmaStaysPositive) {
    NNGenome genome;
    genome.topology = {2, 1};
    genome.weights.resize(nn_weight_count(genome.topology), 0.0f);
    genome.mutation_sigma = 0.01f;  // very small sigma

    for (int i = 0; i < 100; i++) {
        nn_gaussian_mutation(genome);
        EXPECT_GT(genome.mutation_sigma, 0.0f) << "Sigma went non-positive at iteration " << i;
    }
}

TEST(NNMutation, SigmaMutates) {
    NNGenome genome;
    genome.topology = {2, 1};
    genome.weights.resize(nn_weight_count(genome.topology), 0.0f);
    float original_sigma = 0.1f;
    genome.mutation_sigma = original_sigma;

    // After many mutations, sigma should have changed
    bool sigma_changed = false;
    for (int i = 0; i < 50; i++) {
        nn_gaussian_mutation(genome);
        if (std::abs(genome.mutation_sigma - original_sigma) > 1e-6f) {
            sigma_changed = true;
            break;
        }
    }
    EXPECT_TRUE(sigma_changed) << "Sigma didn't self-adapt after 50 mutations";
}

// ============================================================
// T059: Tournament selection test
// ============================================================

TEST(NNTournamentSelection, SelectsHigherFitness) {
    NNPopulation pop;
    pop.topology = {2, 1};
    pop.population_size = 10;
    pop.individuals.resize(10);

    // Lower fitness = better (this is a minimization problem)
    // Set up individuals with known fitness values
    for (int i = 0; i < 10; i++) {
        pop.individuals[i].topology = pop.topology;
        pop.individuals[i].weights.resize(nn_weight_count(pop.topology), 0.0f);
        pop.individuals[i].fitness = static_cast<double>(i) * 100.0;  // 0, 100, 200, ...
    }

    // Over many tournaments, the best individuals should be selected more often
    std::map<int, int> selection_counts;
    int trials = 1000;
    for (int i = 0; i < trials; i++) {
        int selected = nn_tournament_select(pop, 3);
        selection_counts[selected]++;
    }

    // The best individual (fitness=0, index=0) should be selected more than average
    EXPECT_GT(selection_counts[0], trials / 10)
        << "Best individual not selected frequently enough";
}

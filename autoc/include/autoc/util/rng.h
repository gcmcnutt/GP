#ifndef RNG_H
#define RNG_H

#include <random>
#include <cmath>
#include <cstdint>

// Global PRNG — single instance for deterministic evolution.
// Replaces GPrand()/GPsrand()/GPRandomPercent() from libgp.a.
// Seed via rng_seed() at startup for reproducible runs.

namespace rng {

inline std::mt19937& engine() {
    static std::mt19937 gen(42);
    return gen;
}

inline void seed(uint64_t s) {
    engine().seed(static_cast<std::mt19937::result_type>(s));
}

// Uniform integer in [0, max_exclusive)
inline int randInt(int max_exclusive) {
    std::uniform_int_distribution<int> dist(0, max_exclusive - 1);
    return dist(engine());
}

// Uniform long (compatible with old GPrand() range)
inline long randLong() {
    return static_cast<long>(engine()() & 0x7FFFFFFF);
}

// Uniform double in [0, 1)
inline double randDouble() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(engine());
}

// Gaussian N(0, sigma)
inline double randGaussian(double sigma) {
    std::normal_distribution<double> dist(0.0, sigma);
    return dist(engine());
}

// Percentage check: returns true with probability pct/100
inline bool randomPercent(double pct) {
    return randDouble() * 100.0 < pct;
}

} // namespace rng

#endif

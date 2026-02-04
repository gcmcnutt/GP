/**
 * variation_generator.h - Deterministic scenario variation generator
 *
 * Generates entry and wind variations from a seed value and sigma parameters.
 * Called by autoc to populate ScenarioMetadata before sending to crrcsim.
 *
 * All variations are Gaussian-distributed to weight toward nominal conditions
 * while still exposing the GP to off-nominal cases.
 *
 * NOTE: This runs on desktop (64-bit) only, not on flight hardware.
 * Using double precision for accuracy.
 *
 * See specs/VARIATIONS1.md for design rationale.
 */

#ifndef VARIATION_GENERATOR_H
#define VARIATION_GENERATOR_H

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * Sigma parameters for variation generation.
 * Angles in radians, speed as fraction (e.g., 0.1 = 10%).
 */
struct VariationSigmas {
    double headingSigma;      // radians
    double rollSigma;         // radians
    double pitchSigma;        // radians
    double speedSigma;        // fraction (0.1 = ±10% at 1σ)
    double windDirectionSigma; // radians

    // Construct from degrees (for convenience when reading from config)
    static VariationSigmas fromDegrees(double headingDeg, double rollDeg, double pitchDeg,
                                        double speedFrac, double windDirDeg) {
        return VariationSigmas{
            headingDeg * M_PI / 180.0,
            rollDeg * M_PI / 180.0,
            pitchDeg * M_PI / 180.0,
            speedFrac,
            windDirDeg * M_PI / 180.0
        };
    }
};

// Derived variation values - computed from windSeed, stored in ScenarioMetadata
struct VariationOffsets {
    double entryHeadingOffset;  // radians
    double entryRollOffset;     // radians
    double entryPitchOffset;    // radians
    double entrySpeedFactor;    // multiplier (1.0 = nominal)
    double windDirectionOffset; // radians
};

/**
 * Generate all variations from a seed and sigma parameters.
 * Uses simple LCG PRNG for deterministic, reproducible results.
 *
 * @param seed       The random seed (e.g., windSeed from config)
 * @param sigmas     Gaussian sigma values for each variation dimension
 * @return           Computed variation offsets
 */
inline VariationOffsets generateVariations(unsigned int seed, const VariationSigmas& sigmas) {
    VariationOffsets v;

    // LCG PRNG - simple, portable, deterministic
    auto nextDouble = [&seed]() -> double {
        seed = seed * 1103515245 + 12345;
        return static_cast<double>((seed >> 16) & 0x7FFF) / 32768.0;
    };

    // Box-Muller transform for Gaussian sampling
    auto gaussian = [&nextDouble](double sigma) -> double {
        double u1 = nextDouble() * 0.999 + 0.001;  // avoid log(0)
        double u2 = nextDouble();
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        return z * sigma;
    };

    // Generate Gaussian-distributed offsets
    v.entryHeadingOffset = gaussian(sigmas.headingSigma);
    v.entryRollOffset = gaussian(sigmas.rollSigma);
    v.entryPitchOffset = gaussian(sigmas.pitchSigma);
    v.entrySpeedFactor = 1.0 + gaussian(sigmas.speedSigma);
    v.windDirectionOffset = gaussian(sigmas.windDirectionSigma);

    return v;
}

/**
 * Convert radians to degrees for display
 */
inline double radToDeg(double rad) {
    return rad * 180.0 / M_PI;
}

#endif // VARIATION_GENERATOR_H

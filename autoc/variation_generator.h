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
#include <vector>

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

/**
 * Configuration for variable rabbit speed.
 * Set sigma=0 for constant speed at nominal.
 */
struct RabbitSpeedConfig {
    double nominal;    // m/s - center of distribution
    double sigma;      // m/s - 1σ deviation (0 = constant speed)
    double minSpeed;   // m/s - hard floor
    double maxSpeed;   // m/s - hard ceiling
    double cycleMin;   // seconds - min variation cycle duration
    double cycleMax;   // seconds - max variation cycle duration

    static RabbitSpeedConfig defaultConfig() {
        return RabbitSpeedConfig{16.0, 0.0, 8.0, 25.0, 0.5, 5.0};
    }
};

/**
 * A single point in the rabbit speed profile.
 */
struct RabbitSpeedPoint {
    double timeSec;
    double speed;
};

/**
 * Generate a speed profile for the rabbit.
 * Uses the same LCG PRNG as generateVariations() - call after generateVariations()
 * to continue the PRNG sequence.
 *
 * @param seed              Reference to seed - continues PRNG sequence
 * @param cfg               Speed configuration parameters
 * @param totalDurationSec  Total path duration in seconds
 * @return                  Vector of time-tagged speed points
 */
inline std::vector<RabbitSpeedPoint> generateSpeedProfile(
    unsigned int& seed,
    const RabbitSpeedConfig& cfg,
    double totalDurationSec
) {
    std::vector<RabbitSpeedPoint> profile;

    // Constant speed mode - just return two points at nominal speed
    if (cfg.sigma <= 0.0) {
        profile.push_back({0.0, cfg.nominal});
        profile.push_back({totalDurationSec, cfg.nominal});
        return profile;
    }

    // LCG PRNG - continues from previous seed state
    auto nextDouble = [&seed]() -> double {
        seed = seed * 1103515245 + 12345;
        return static_cast<double>((seed >> 16) & 0x7FFF) / 32768.0;
    };

    // Box-Muller for Gaussian sampling
    auto gaussian = [&nextDouble](double mean, double sigma) -> double {
        double u1 = nextDouble() * 0.999 + 0.001;  // avoid log(0)
        double u2 = nextDouble();
        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        return mean + z * sigma;
    };

    // Clamp speed to valid range
    auto clampSpeed = [&cfg](double s) -> double {
        return (s < cfg.minSpeed) ? cfg.minSpeed :
               (s > cfg.maxSpeed) ? cfg.maxSpeed : s;
    };

    double t = 0.0;
    double currentSpeed = clampSpeed(gaussian(cfg.nominal, cfg.sigma));
    profile.push_back({t, currentSpeed});

    while (t < totalDurationSec) {
        // Random cycle duration
        double cycleDuration = cfg.cycleMin + nextDouble() * (cfg.cycleMax - cfg.cycleMin);
        double targetSpeed = clampSpeed(gaussian(cfg.nominal, cfg.sigma));

        // Cosine-eased interpolation points (~10Hz resolution)
        int numSteps = static_cast<int>(cycleDuration / 0.1);
        if (numSteps < 2) numSteps = 2;

        for (int i = 1; i <= numSteps && (t + i * cycleDuration / numSteps) <= totalDurationSec; i++) {
            double frac = static_cast<double>(i) / numSteps;
            double easedFrac = 0.5 * (1.0 - cos(M_PI * frac));  // Cosine ease-in-out
            double interpTime = t + frac * cycleDuration;
            double interpSpeed = currentSpeed + easedFrac * (targetSpeed - currentSpeed);
            profile.push_back({interpTime, interpSpeed});
        }

        t += cycleDuration;
        currentSpeed = targetSpeed;
    }

    return profile;
}

/**
 * Get speed at arbitrary time by interpolating the speed profile.
 *
 * @param profile   Speed profile from generateSpeedProfile()
 * @param timeSec   Time in seconds
 * @return          Interpolated speed in m/s
 */
inline double getSpeedAtTime(const std::vector<RabbitSpeedPoint>& profile, double timeSec) {
    if (profile.empty()) return 16.0;  // Fallback to default
    if (timeSec <= profile.front().timeSec) return profile.front().speed;
    if (timeSec >= profile.back().timeSec) return profile.back().speed;

    // Binary search for bracket
    size_t lo = 0, hi = profile.size() - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (profile[mid].timeSec <= timeSec) lo = mid;
        else hi = mid;
    }

    // Linear interpolation
    double t0 = profile[lo].timeSec, t1 = profile[hi].timeSec;
    double s0 = profile[lo].speed, s1 = profile[hi].speed;
    double frac = (timeSec - t0) / (t1 - t0);
    return s0 + frac * (s1 - s0);
}

#endif // VARIATION_GENERATOR_H

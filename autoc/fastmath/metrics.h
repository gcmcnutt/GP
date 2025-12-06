#ifndef FASTMATH_METRICS_H
#define FASTMATH_METRICS_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <ostream>

namespace fastmath {

struct FastMathMetrics {
    uint64_t gp_scalar_saturations = 0;
    uint64_t range_limit_clamps = 0;
    uint64_t zero_snap_events = 0;
};

inline FastMathMetrics& getFastMathMetrics() {
    static FastMathMetrics metrics;
    return metrics;
}

inline void resetFastMathMetrics() {
    getFastMathMetrics() = FastMathMetrics{};
}

#ifdef GP_FASTMATH_TRACE
inline void recordScalarSaturation() {
    ++getFastMathMetrics().gp_scalar_saturations;
}

inline void recordRangeClamp() {
    ++getFastMathMetrics().range_limit_clamps;
}

inline void recordZeroSnap() {
    ++getFastMathMetrics().zero_snap_events;
}
#else
inline void recordScalarSaturation() {}
inline void recordRangeClamp() {}
inline void recordZeroSnap() {}
#endif

#ifdef GP_FASTMATH_TRACE
inline void logFastMathMetrics(const char* label, std::ostream& os, double divisor = 0.0) {
    const auto& metrics = getFastMathMetrics();
    if (metrics.gp_scalar_saturations == 0 &&
        metrics.range_limit_clamps == 0 &&
        metrics.zero_snap_events == 0) {
        return;
    }
    os << "[FastMath] " << label
       << " sat=" << metrics.gp_scalar_saturations
       << " clamp=" << metrics.range_limit_clamps
       << " zero=" << metrics.zero_snap_events;
    if (divisor > 0.0) {
        os << " | per_eval sat=" << (metrics.gp_scalar_saturations / divisor)
           << " clamp=" << (metrics.range_limit_clamps / divisor)
           << " zero=" << (metrics.zero_snap_events / divisor);
    }
    os << std::endl;
    resetFastMathMetrics();
}
#else
inline void logFastMathMetrics(const char*, std::ostream&, double = 0.0) {}
#endif

struct FastMathComparisonStats {
    uint64_t samples = 0;
    uint64_t violations = 0;
    double max_pitch_error = 0.0;
    double max_roll_error = 0.0;
    double max_throttle_error = 0.0;
    double max_result_error = 0.0;
};

inline FastMathComparisonStats& getFastMathComparisonStats() {
    static FastMathComparisonStats stats;
    return stats;
}

inline void resetFastMathComparisonStats() {
    getFastMathComparisonStats() = FastMathComparisonStats{};
}

inline bool isFastMathComparisonEnabled() {
    static int enabled = []() {
        const char* env = std::getenv("GP_FASTMATH_COMPARE");
        if (!env) {
            return 0;
        }
        return (std::strcmp(env, "0") == 0 || std::strcmp(env, "false") == 0 ||
                std::strcmp(env, "FALSE") == 0)
                   ? 0
                   : 1;
    }();
    return enabled != 0;
}

inline double getFastMathComparisonTolerance() {
    static double tol = []() {
        const char* env = std::getenv("GP_FASTMATH_COMPARE_TOL");
        if (!env) {
            return 0.01;
        }
        double value = std::atof(env);
        return (value <= 0.0) ? 0.01 : value;
    }();
    return tol;
}

inline void recordFastMathComparison(double pitch_scalar, double pitch_ref,
                                     double roll_scalar, double roll_ref,
                                     double throttle_scalar, double throttle_ref,
                                     double result_scalar, double result_ref) {
    if (!isFastMathComparisonEnabled()) {
        return;
    }
    auto& stats = getFastMathComparisonStats();
    stats.samples++;
    double tol = getFastMathComparisonTolerance();
    bool violated = false;

    auto trackError = [&](double error, double& maxField) {
        if (error > maxField) {
            maxField = error;
        }
        if (error > tol) {
            violated = true;
        }
    };

    trackError(std::abs(pitch_scalar - pitch_ref), stats.max_pitch_error);
    trackError(std::abs(roll_scalar - roll_ref), stats.max_roll_error);
    trackError(std::abs(throttle_scalar - throttle_ref), stats.max_throttle_error);
    trackError(std::abs(result_scalar - result_ref), stats.max_result_error);

    if (violated) {
        stats.violations++;
    }
}

inline void logFastMathComparison(const char* label, std::ostream& os) {
    if (!isFastMathComparisonEnabled()) {
        return;
    }
    const auto& stats = getFastMathComparisonStats();
    if (stats.samples == 0) {
        return;
    }
    os << "[FastMathCompare] " << label
       << " samples=" << stats.samples
       << " violations=" << stats.violations
       << " tol=" << getFastMathComparisonTolerance()
       << " pitch_max=" << stats.max_pitch_error
       << " roll_max=" << stats.max_roll_error
       << " throttle_max=" << stats.max_throttle_error
       << " result_max=" << stats.max_result_error
       << std::endl;
    resetFastMathComparisonStats();
}

}  // namespace fastmath

#endif  // FASTMATH_METRICS_H

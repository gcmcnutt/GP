#ifndef FASTMATH_METRICS_H
#define FASTMATH_METRICS_H

#include <cstdint>
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

}  // namespace fastmath

#endif  // FASTMATH_METRICS_H

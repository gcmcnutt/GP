// sensor_math.cc — Navigation sensor functions for temporal history recording
// Extracted from gp_evaluator_portable during GP removal (014-nn-training-signal)
//
// Uses LUT-based trig (fastAtan2) matching the original GP evaluator precision.
// This matters because:
//   1. NN was trained with history values computed by fastAtan2
//   2. xiao-gp embedded platform needs LUT trig for speed
//   3. Desktop/embedded must produce identical sensor values

#include "autoc/eval/sensor_math.h"
#include <cmath>
#include <array>

// ============================================================================
// LUT-based atan2 — copied from gp_evaluator_portable.cc
// 512-entry table, interpolated, ~0.001 rad max error vs std::atan2
// ============================================================================

namespace {

constexpr gp_scalar GP_PI = static_cast<gp_scalar>(M_PI);
constexpr gp_scalar GP_HALF_PI = GP_PI * static_cast<gp_scalar>(0.5f);

constexpr int ATAN_LUT_SIZE = 512;
std::array<gp_scalar, ATAN_LUT_SIZE + 1> ATAN_LUT{};
bool ATAN_LUT_INIT = false;

inline void initAtanLut() {
    if (ATAN_LUT_INIT) return;
    for (int i = 0; i <= ATAN_LUT_SIZE; ++i) {
        gp_scalar t = static_cast<gp_scalar>(i) / static_cast<gp_scalar>(ATAN_LUT_SIZE);
        ATAN_LUT[i] = std::atan(t);
    }
    ATAN_LUT_INIT = true;
}

inline gp_scalar fastAtan(gp_scalar r) {
    initAtanLut();
    gp_scalar t = ABS_DEF(r);
    gp_scalar scaled = t * static_cast<gp_scalar>(ATAN_LUT_SIZE);
    if (scaled >= ATAN_LUT_SIZE) {
        return (r < 0.0f ? -ATAN_LUT[ATAN_LUT_SIZE] : ATAN_LUT[ATAN_LUT_SIZE]);
    }
    int i0 = static_cast<int>(scaled);
    gp_scalar frac = scaled - static_cast<gp_scalar>(i0);
    gp_scalar a0 = ATAN_LUT[i0];
    gp_scalar a1 = ATAN_LUT[i0 + 1];
    gp_scalar angle = a0 + (a1 - a0) * frac;
    return (r < 0.0f) ? -angle : angle;
}

inline gp_scalar fastAtan2(gp_scalar y, gp_scalar x) {
    if (ABS_DEF(x) < static_cast<gp_scalar>(1e-6f) && ABS_DEF(y) < static_cast<gp_scalar>(1e-6f)) {
        return 0.0f;
    }

    if (ABS_DEF(x) > ABS_DEF(y)) {
        gp_scalar angle = fastAtan(y / x);
        if (x < 0.0f) {
            angle += (y >= 0.0f ? GP_PI : -GP_PI);
        }
        return angle;
    } else {
        gp_scalar angle = fastAtan(x / y);
        gp_scalar base = (y > 0.0f) ? GP_HALF_PI : -GP_HALF_PI;
        return base - angle;
    }
}

}  // namespace

// ============================================================================
// Path interpolation
// ============================================================================

gp_vec3 getInterpolatedTargetPosition(PathProvider& pathProvider,
                                       int32_t currentTimeMsec,
                                       gp_scalar offsetSteps) {
    int pathSize = pathProvider.getPathSize();
    if (pathSize == 0) {
        return gp_vec3::Zero();
    }
    if (pathSize == 1) {
        return pathProvider.getPath(0).start;
    }

    // Handle NaN offset — return current rabbit position
    if (std::isnan(offsetSteps)) {
        int currentIdx = CLAMP_DEF(pathProvider.getCurrentIndex(), 0, pathSize - 1);
        return pathProvider.getPath(currentIdx).start;
    }

    // Clamp offset to ±MAX_OFFSET_STEPS (±1 second)
    gp_scalar clampedSteps = CLAMP_DEF(offsetSteps,
                                        static_cast<gp_scalar>(-MAX_OFFSET_STEPS),
                                        static_cast<gp_scalar>(MAX_OFFSET_STEPS));

    // Calculate goal time as integer (deterministic, no float precision issues)
    int32_t goalTimeMsec = currentTimeMsec + static_cast<int32_t>(clampedSteps * static_cast<gp_scalar>(SIM_TIME_STEP_MSEC));

    // Clamp goal time to path bounds (all int32_t comparisons — exact)
    int32_t minTime = pathProvider.getPath(0).simTimeMsec;
    int32_t maxTime = pathProvider.getPath(pathSize - 1).simTimeMsec;

    if (goalTimeMsec <= minTime) {
        return pathProvider.getPath(0).start;
    }
    if (goalTimeMsec >= maxTime) {
        return pathProvider.getPath(pathSize - 1).start;
    }

    // Binary search for bracketing waypoints (int32_t vs int32_t — deterministic)
    int lo = 0;
    int hi = pathSize - 1;
    while (lo < hi - 1) {
        int mid = (lo + hi) / 2;
        if (pathProvider.getPath(mid).simTimeMsec <= goalTimeMsec) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const Path& p0 = pathProvider.getPath(lo);
    const Path& p1 = pathProvider.getPath(lo + 1);

    // Calculate interpolation fraction
    int32_t dt = p1.simTimeMsec - p0.simTimeMsec;
    gp_scalar frac = 0.0f;
    if (dt > 0) {
        frac = static_cast<gp_scalar>(goalTimeMsec - p0.simTimeMsec) / static_cast<gp_scalar>(dt);
        frac = CLAMP_DEF(frac, 0.0f, 1.0f);
    }

    return p0.start + frac * (p1.start - p0.start);
}

// ============================================================================
// Navigation sensors — angle to target in body frame
// Uses fastAtan2 for LUT-based trig matching embedded platform
// ============================================================================

gp_scalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {
    gp_vec3 targetPos = getInterpolatedTargetPosition(
        pathProvider, static_cast<int32_t>(aircraftState.getSimTimeMsec()), arg);

    gp_vec3 craftToTarget = targetPos - aircraftState.getPosition();
    gp_vec3 target_local = aircraftState.getOrientation().inverse() * craftToTarget;

    // Project onto body YZ plane, angle from body -Z axis
    return fastAtan2(target_local.y(), -target_local.z());
}

gp_scalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {
    gp_vec3 targetPos = getInterpolatedTargetPosition(
        pathProvider, static_cast<int32_t>(aircraftState.getSimTimeMsec()), arg);

    gp_vec3 craftToTarget = targetPos - aircraftState.getPosition();
    gp_vec3 target_local = aircraftState.getOrientation().inverse() * craftToTarget;

    // Pitch error: rotation about body Y to point nose at target
    return fastAtan2(-target_local.z(), target_local.x());
}

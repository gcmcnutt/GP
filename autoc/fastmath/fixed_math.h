#ifndef FASTMATH_FIXED_MATH_H
#define FASTMATH_FIXED_MATH_H

#include <cmath>
#include "gp_scalar.h"

namespace fastmath {
namespace detail {

constexpr int kTrigTableSize = 1024;
constexpr double kTwoPi = 2.0 * M_PI;
constexpr double kInvTwoPi = 1.0 / kTwoPi;
constexpr double kResolution = static_cast<double>(kTrigTableSize) / kTwoPi;

struct TrigTables {
    GPScalar sinTable[kTrigTableSize + 1];

    constexpr TrigTables() : sinTable{} {
        // constexpr constructor can't use std::sin, so populate at runtime via helper.
    }
};

inline const TrigTables& getTrigTables() {
    static TrigTables tables;
    static bool initialized = false;
    if (!initialized) {
        for (int i = 0; i <= kTrigTableSize; ++i) {
            double angle = (static_cast<double>(i) / kTrigTableSize) * kTwoPi;
            tables.sinTable[i] = GPScalar::fromDouble(std::sin(angle));
        }
        initialized = true;
    }
    return tables;
}

inline double normalizeRadians(double radians) {
    double revolutions = radians * kInvTwoPi;
    double fractional = revolutions - std::floor(revolutions);
    if (fractional < 0.0) {
        fractional += 1.0;
    }
    return fractional * kTwoPi;
}

}  // namespace detail

inline GPScalar sinApprox(GPScalar angle) {
    double normalized = detail::normalizeRadians(angle.toDouble());
    double scaled_index = normalized * detail::kResolution;
    int base_index = static_cast<int>(scaled_index);
    if (base_index >= detail::kTrigTableSize) {
        base_index = 0;
        scaled_index -= detail::kTrigTableSize;
    }
    double frac = scaled_index - static_cast<double>(base_index);
    const auto& tables = detail::getTrigTables();
    const GPScalar& y0 = tables.sinTable[base_index];
    const GPScalar& y1 = tables.sinTable[base_index + 1];
    GPScalar delta = y1 - y0;
    GPScalar frac_scalar = GPScalar::fromDouble(frac);
    return y0 + delta * frac_scalar;
}

inline GPScalar cosApprox(GPScalar angle) {
    static const GPScalar kHalfPi = GPScalar::fromDouble(M_PI_2);
    return sinApprox(angle + kHalfPi);
}

inline GPScalar atan2Approx(GPScalar y, GPScalar x) {
    double yd = y.toDouble();
    double xd = x.toDouble();
    if (xd == 0.0 && yd == 0.0) {
        return GPScalar::zero();
    }

    double abs_y = std::abs(yd) + 1e-12;  // Prevent division by zero
    double angle = 0.0;
    if (xd >= 0.0) {
        double r = (xd - abs_y) / (xd + abs_y);
        angle = M_PI_4 - M_PI_4 * r;
    } else {
        double r = (xd + abs_y) / (abs_y - xd);
        angle = 3.0 * M_PI_4 - M_PI_4 * r;
    }
    if (yd < 0.0) {
        angle = -angle;
    }
    return GPScalar::fromDouble(angle);
}

inline GPScalar sqrtApprox(GPScalar value) {
    double v = value.toDouble();
    if (v <= 0.0) {
        return GPScalar::zero();
    }

    int exponent;
    double mantissa = std::frexp(v, &exponent);
    double approx = 0.41731 + 0.59016 * mantissa;
    approx = std::ldexp(approx, exponent / 2);
    for (int i = 0; i < 2; ++i) {
        approx = 0.5 * (approx + v / approx);
    }
    return GPScalar::fromDouble(approx);
}

}  // namespace fastmath

#endif  // FASTMATH_FIXED_MATH_H

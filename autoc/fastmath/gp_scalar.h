#ifndef FASTMATH_GP_SCALAR_H
#define FASTMATH_GP_SCALAR_H

#include <cstdint>
#include <limits>
#include <cmath>

#include "metrics.h"

namespace fastmath {

class GPScalar {
public:
    static constexpr int kFractionBits = 15;
    static constexpr int32_t kScale = 1 << kFractionBits;
    static constexpr int32_t kMinRaw = std::numeric_limits<int32_t>::min();
    static constexpr int32_t kMaxRaw = std::numeric_limits<int32_t>::max();

    constexpr GPScalar() : value_(0) {}
    explicit constexpr GPScalar(int32_t raw) : value_(raw) {}
    GPScalar(double value) : value_(convertDoubleToRaw(value)) {}

    static constexpr GPScalar fromRaw(int32_t raw) { return GPScalar(raw); }
    static constexpr GPScalar fromInt(int32_t value) {
        return GPScalar::fromRaw(static_cast<int32_t>(value) << kFractionBits);
    }

    static GPScalar fromDouble(double value) {
        return GPScalar(convertDoubleToRaw(value));
    }

    double toDouble() const {
        return static_cast<double>(value_) / static_cast<double>(kScale);
    }

    constexpr GPScalar operator-() const {
        return GPScalar::fromRaw(value_ == kMinRaw ? kMaxRaw : -value_);
    }

    GPScalar operator+(const GPScalar& rhs) const {
        return GPScalar::fromRaw(saturatingAdd(value_, rhs.value_));
    }

    GPScalar operator-(const GPScalar& rhs) const {
        return GPScalar::fromRaw(saturatingSub(value_, rhs.value_));
    }

    GPScalar operator*(const GPScalar& rhs) const {
        int64_t product = static_cast<int64_t>(value_) * static_cast<int64_t>(rhs.value_);
        int64_t bias = (product >= 0) ? (1LL << (kFractionBits - 1)) : -(1LL << (kFractionBits - 1));
        int64_t shifted = (product + bias) >> kFractionBits;
        return GPScalar::fromRaw(saturatingClamp(shifted));
    }

    GPScalar operator/(const GPScalar& rhs) const {
        if (rhs.value_ == 0) {
            return GPScalar::zero();
        }
        int64_t numerator = static_cast<int64_t>(value_) << kFractionBits;
        int64_t quotient = numerator / rhs.value_;
        return GPScalar::fromRaw(saturatingClamp(quotient));
    }

    GPScalar& operator+=(const GPScalar& rhs) {
        value_ = saturatingAdd(value_, rhs.value_);
        return *this;
    }

    GPScalar& operator-=(const GPScalar& rhs) {
        value_ = saturatingSub(value_, rhs.value_);
        return *this;
    }

    GPScalar& operator*=(const GPScalar& rhs) {
        *this = *this * rhs;
        return *this;
    }

    GPScalar& operator/=(const GPScalar& rhs) {
        *this = *this / rhs;
        return *this;
    }

    friend constexpr bool operator==(const GPScalar& lhs, const GPScalar& rhs) {
        return lhs.value_ == rhs.value_;
    }

    friend constexpr bool operator!=(const GPScalar& lhs, const GPScalar& rhs) {
        return !(lhs == rhs);
    }

    friend constexpr bool operator<(const GPScalar& lhs, const GPScalar& rhs) {
        return lhs.value_ < rhs.value_;
    }

    friend constexpr bool operator>(const GPScalar& lhs, const GPScalar& rhs) {
        return rhs < lhs;
    }

    friend constexpr bool operator<=(const GPScalar& lhs, const GPScalar& rhs) {
        return !(rhs < lhs);
    }

    friend constexpr bool operator>=(const GPScalar& lhs, const GPScalar& rhs) {
        return !(lhs < rhs);
    }

    constexpr bool isZero() const { return value_ == 0; }

    constexpr GPScalar abs() const {
        return (value_ >= 0) ? *this : GPScalar::fromRaw(value_ == kMinRaw ? kMaxRaw : -value_);
    }

    constexpr int32_t raw() const { return value_; }

    static constexpr GPScalar zero() { return GPScalar::fromRaw(0); }

    operator double() const { return toDouble(); }

private:
    static inline int32_t saturatingAdd(int32_t a, int32_t b) {
        int64_t sum = static_cast<int64_t>(a) + static_cast<int64_t>(b);
        return saturatingClamp(sum);
    }

    static inline int32_t saturatingSub(int32_t a, int32_t b) {
        int64_t diff = static_cast<int64_t>(a) - static_cast<int64_t>(b);
        return saturatingClamp(diff);
    }

    static inline int32_t saturatingClamp(int64_t value) {
        if (value > static_cast<int64_t>(kMaxRaw)) {
            recordScalarSaturation();
            return kMaxRaw;
        }
        if (value < static_cast<int64_t>(kMinRaw)) {
            recordScalarSaturation();
            return kMinRaw;
        }
        return static_cast<int32_t>(value);
    }

    static int32_t convertDoubleToRaw(double value) {
        double scaled = value * static_cast<double>(kScale);
        if (scaled > static_cast<double>(kMaxRaw)) {
            recordScalarSaturation();
            scaled = static_cast<double>(kMaxRaw);
        } else if (scaled < static_cast<double>(kMinRaw)) {
            recordScalarSaturation();
            scaled = static_cast<double>(kMinRaw);
        }
        auto rounded = static_cast<int64_t>(std::llround(scaled));
        if (rounded > kMaxRaw) {
            recordScalarSaturation();
            rounded = kMaxRaw;
        }
        if (rounded < kMinRaw) {
            recordScalarSaturation();
            rounded = kMinRaw;
        }
        return static_cast<int32_t>(rounded);
    }

    int32_t value_;
};

inline GPScalar clamp(GPScalar value, GPScalar min, GPScalar max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

inline GPScalar min(GPScalar a, GPScalar b) {
    return (a < b) ? a : b;
}

inline GPScalar max(GPScalar a, GPScalar b) {
    return (a > b) ? a : b;
}

} // namespace fastmath

#endif

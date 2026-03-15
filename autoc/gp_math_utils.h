#ifndef GP_MATH_UTILS_H
#define GP_MATH_UTILS_H

#include "rng.h"

// Backward-compatible wrappers — these delegate to rng:: namespace.
// TODO: Remove once all callers use rng:: directly.

inline double GPrandGaussian(double sigma) {
    return rng::randGaussian(sigma);
}

inline double GPrandDouble() {
    return rng::randDouble();
}

#endif

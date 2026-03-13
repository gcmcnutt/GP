#ifndef GP_MATH_UTILS_H
#define GP_MATH_UTILS_H

#include <cmath>
#include "gp.h"  // For GPrand()

// Box-Muller Gaussian using GPrand() as the uniform source.
// Matches the pattern in variation_generator.h but as a reusable function.
// Returns a sample from N(0, sigma).
inline double GPrandGaussian(double sigma) {
    // GPrand() returns long in [0, 2147483646]
    // Convert to uniform double in (0, 1)
    constexpr double RNMX = 2147483646.0;
    double u1 = (GPrand() / RNMX) * 0.999 + 0.001;  // avoid log(0)
    double u2 = GPrand() / RNMX;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z * sigma;
}

// Uniform double from GPrand() in [0, 1)
inline double GPrandDouble() {
    constexpr double RNMX = 2147483646.0;
    return GPrand() / RNMX;
}

#endif

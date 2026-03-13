#include "fitness_computer.h"
#include <cmath>
#include <algorithm>

// Extracted from autoc.cc lines 1288-1290
// Per-step combined distance + attitude penalty
double FitnessComputer::computeStepPenalty(double distance, double attitude_delta, double intercept_scale) {
    double distDelta = fabs(distance - DISTANCE_TARGET);
    double dist_penalty = pow(intercept_scale * distDelta / DISTANCE_NORM, DISTANCE_POWER);
    double att_penalty = pow(intercept_scale * attitude_delta / ATTITUDE_NORM, ATTITUDE_POWER);
    return dist_penalty + att_penalty;
}

// Extracted from autoc.cc lines 1402-1408
// Crash penalty: soft lexicographic - completion dominates
double FitnessComputer::computeCrashPenalty(double fraction_completed) {
    return (1.0 - fraction_completed) * CRASH_COMPLETION_WEIGHT;
}

// Extracted from autoc.cc lines 1389-1393
// Attitude scaling from path geometry
double FitnessComputer::computeAttitudeScale(double path_distance, double path_turn_rad) {
    double denominator = std::max(path_turn_rad, 2.0 * M_PI);
    return path_distance / denominator;
}

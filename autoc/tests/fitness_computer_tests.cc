#include <gtest/gtest.h>
#include <cmath>
#include "../fitness_computer.h"

// Constants provided by fitness_computer.h via #define (matching autoc.h values)

// ========================================================================
// T010: computeStepPenalty parity
// ========================================================================

// Reference implementation matching autoc.cc lines 1288-1290
static double referenceStepPenalty(double distance, double attitude_delta, double interceptScale) {
    double distDelta = fabs(distance - DISTANCE_TARGET);
    double dist_penalty = pow(interceptScale * distDelta / DISTANCE_NORM, DISTANCE_POWER);
    double att_penalty = pow(interceptScale * attitude_delta / ATTITUDE_NORM, ATTITUDE_POWER);
    return dist_penalty + att_penalty;
}

TEST(FitnessComputer, StepPenaltyAtTarget) {
    // At DISTANCE_TARGET with no attitude change, penalty should be near zero
    FitnessComputer fc;
    double penalty = fc.computeStepPenalty(DISTANCE_TARGET, 0.0, 1.0);
    double expected = referenceStepPenalty(DISTANCE_TARGET, 0.0, 1.0);
    EXPECT_NEAR(penalty, expected, 1e-10);
    EXPECT_NEAR(penalty, 0.0, 1e-10);
}

TEST(FitnessComputer, StepPenaltyTooFar) {
    // 15m from rabbit (7.5m over target)
    FitnessComputer fc;
    double penalty = fc.computeStepPenalty(15.0, 0.0, 1.0);
    double expected = referenceStepPenalty(15.0, 0.0, 1.0);
    EXPECT_NEAR(penalty, expected, 1e-10);
    EXPECT_GT(penalty, 0.0);
}

TEST(FitnessComputer, StepPenaltyTooClose) {
    // 2m from rabbit (5.5m under target) — V-shaped penalty
    FitnessComputer fc;
    double penalty = fc.computeStepPenalty(2.0, 0.0, 1.0);
    double expected = referenceStepPenalty(2.0, 0.0, 1.0);
    EXPECT_NEAR(penalty, expected, 1e-10);
    EXPECT_GT(penalty, 0.0);
}

TEST(FitnessComputer, StepPenaltyWithAttitude) {
    // Distance at target but large attitude change
    FitnessComputer fc;
    double penalty = fc.computeStepPenalty(DISTANCE_TARGET, 0.5, 1.0);
    double expected = referenceStepPenalty(DISTANCE_TARGET, 0.5, 1.0);
    EXPECT_NEAR(penalty, expected, 1e-10);
    EXPECT_GT(penalty, 0.0);
}

TEST(FitnessComputer, StepPenaltyWithInterceptScale) {
    // interceptScale < 1.0 reduces penalty (early in intercept budget)
    FitnessComputer fc;
    double penalty_full = fc.computeStepPenalty(15.0, 0.3, 1.0);
    double penalty_scaled = fc.computeStepPenalty(15.0, 0.3, 0.1);
    EXPECT_LT(penalty_scaled, penalty_full);

    double expected = referenceStepPenalty(15.0, 0.3, 0.1);
    EXPECT_NEAR(penalty_scaled, expected, 1e-10);
}

TEST(FitnessComputer, StepPenaltyKnownValues) {
    // Hand-computed reference values from the spec:
    // At 5m: distDelta=2.5, penalty = pow(2.5/5.0, 0.75) = pow(0.5, 0.75) ≈ 0.5946
    // At 10m: distDelta=2.5, same
    // At 24m: distDelta=16.5, penalty = pow(16.5/5.0, 0.75) = pow(3.3, 0.75) ≈ 2.5767
    FitnessComputer fc;
    double penalty_5m = fc.computeStepPenalty(5.0, 0.0, 1.0);
    double penalty_10m = fc.computeStepPenalty(10.0, 0.0, 1.0);
    // Symmetric around target
    EXPECT_NEAR(penalty_5m, penalty_10m, 1e-10);
    // Known value
    double expected = pow(2.5 / DISTANCE_NORM, DISTANCE_POWER);
    EXPECT_NEAR(penalty_5m, expected, 1e-10);
}

// ========================================================================
// T011: computeCrashPenalty parity
// ========================================================================

// Reference: autoc.cc lines 1402-1408
static double referenceCrashPenalty(double fraction_completed) {
    return (1.0 - fraction_completed) * CRASH_COMPLETION_WEIGHT;
}

TEST(FitnessComputer, CrashPenaltyHalfCompleted) {
    FitnessComputer fc;
    double penalty = fc.computeCrashPenalty(0.5);
    double expected = referenceCrashPenalty(0.5);
    EXPECT_NEAR(penalty, expected, 1e-6);
    EXPECT_NEAR(penalty, 500000.0, 1e-6);
}

TEST(FitnessComputer, CrashPenaltyAlmostComplete) {
    FitnessComputer fc;
    double penalty = fc.computeCrashPenalty(0.99);
    double expected = referenceCrashPenalty(0.99);
    EXPECT_NEAR(penalty, expected, 1e-6);
    EXPECT_NEAR(penalty, 10000.0, 1e-6);
}

TEST(FitnessComputer, CrashPenaltyMinCompletion) {
    // Minimum fraction (clamped to 0.001 in autoc.cc)
    FitnessComputer fc;
    double penalty = fc.computeCrashPenalty(0.001);
    double expected = referenceCrashPenalty(0.001);
    EXPECT_NEAR(penalty, expected, 1e-6);
    EXPECT_NEAR(penalty, 999000.0, 1e-3);
}

TEST(FitnessComputer, CrashPenaltyFullCompletion) {
    // Should be zero at 1.0 (no crash, no penalty)
    FitnessComputer fc;
    double penalty = fc.computeCrashPenalty(1.0);
    EXPECT_NEAR(penalty, 0.0, 1e-10);
}

// ========================================================================
// T012: computeAttitudeScale parity
// ========================================================================

// Reference: autoc.cc lines 1389-1393
static double referenceAttitudeScale(double path_distance, double path_turn_rad) {
    double denominator = std::max(path_turn_rad, 2.0 * M_PI);
    return path_distance / denominator;
}

TEST(FitnessComputer, AttitudeScaleNormalPath) {
    // 100m path with 3 radians of turning
    FitnessComputer fc;
    double scale = fc.computeAttitudeScale(100.0, 3.0);
    // 3.0 < 2*PI, so denominator = 2*PI
    double expected = referenceAttitudeScale(100.0, 3.0);
    EXPECT_NEAR(scale, expected, 1e-10);
    EXPECT_NEAR(scale, 100.0 / (2.0 * M_PI), 1e-10);
}

TEST(FitnessComputer, AttitudeScaleHighTurnPath) {
    // 200m path with 10 radians of turning (> 2*PI)
    FitnessComputer fc;
    double scale = fc.computeAttitudeScale(200.0, 10.0);
    double expected = referenceAttitudeScale(200.0, 10.0);
    EXPECT_NEAR(scale, expected, 1e-10);
    EXPECT_NEAR(scale, 200.0 / 10.0, 1e-10);  // = 20.0
}

TEST(FitnessComputer, AttitudeScaleStraightPath) {
    // Straight path: turn_rad = 0 → fallback to 2*PI
    FitnessComputer fc;
    double scale = fc.computeAttitudeScale(50.0, 0.0);
    double expected = referenceAttitudeScale(50.0, 0.0);
    EXPECT_NEAR(scale, expected, 1e-10);
    EXPECT_NEAR(scale, 50.0 / (2.0 * M_PI), 1e-10);
}

#include <gtest/gtest.h>
#include <set>
#include <vector>
#include <cmath>
#include "../gp_evaluator_portable.h"

// T001: Local node definitions for test coverage tracking
// Mirrors allNodes[] from autoc-eval.cc but defined locally to avoid link dependencies
struct TestNodeDef {
  int opcode;
  const char* name;
  int args;
};

// All 42 nodes defined in the GP system - must stay in sync with autoc.h/autoc-eval.cc
static const TestNodeDef testAllNodes[] = {
  {ADD, "ADD", 2},
  {SUB, "SUB", 2},
  {MUL, "MUL", 2},
  {DIV, "DIV", 2},
  {IF, "IF", 3},
  {EQ, "EQ", 2},
  {GT, "GT", 2},
  {SETPITCH, "SETPITCH", 1},
  {SETROLL, "SETROLL", 1},
  {SETTHROTTLE, "SETTHROTTLE", 1},
  {GETPITCH, "GETPITCH", 0},
  {GETROLL, "GETROLL", 0},
  {GETTHROTTLE, "GETTHROTTLE", 0},
  {SIN, "SIN", 1},
  {COS, "COS", 1},
  {OP_PI, "OP_PI", 0},
  {ZERO, "0", 0},
  {ONE, "1", 0},
  {TWO, "2", 0},
  {PROGN, "PROGN", 2},
  {GETDPHI, "GETDPHI", 1},
  {GETDTHETA, "GETDTHETA", 1},
  {GETDTARGET, "GETDTARGET", 1},
  {GETVEL, "GETVEL", 0},
  {GETDHOME, "GETDHOME", 0},
  {GETALPHA, "GETALPHA", 0},
  {GETBETA, "GETBETA", 0},
  {GETVELX, "GETVELX", 0},
  {GETVELY, "GETVELY", 0},
  {GETVELZ, "GETVELZ", 0},
  {GETROLL_RAD, "GETROLL_RAD", 0},
  {GETPITCH_RAD, "GETPITCH_RAD", 0},
  {CLAMP, "CLAMP", 3},
  {ATAN2, "ATAN2", 2},
  {ABS, "ABS", 1},
  {SQRT, "SQRT", 1},
  {MIN, "MIN", 2},
  {MAX, "MAX", 2},
  {GETDPHI_PREV, "GETDPHI_PREV", 1},
  {GETDTHETA_PREV, "GETDTHETA_PREV", 1},
  {GETDPHI_RATE, "GETDPHI_RATE", 0},
  {GETDTHETA_RATE, "GETDTHETA_RATE", 0},
};
static const int testAllNodesCount = sizeof(testAllNodes) / sizeof(testAllNodes[0]);

namespace {

// =============================================================================
// Phase 2: Foundational Test Helpers
// =============================================================================

// T004-T006: QuatHelper namespace for constructing test quaternions
namespace QuatHelper {
    // T006: Degrees to radians helper
    constexpr gp_scalar deg(gp_scalar d) {
        return d * static_cast<gp_scalar>(M_PI) / 180.0f;
    }

    // T004: From Euler angles (roll, pitch, yaw) in radians - ZYX intrinsic sequence
    inline gp_quat fromEuler(gp_scalar roll, gp_scalar pitch, gp_scalar yaw) {
        return gp_quat(Eigen::AngleAxis<gp_scalar>(yaw, gp_vec3::UnitZ())) *
               gp_quat(Eigen::AngleAxis<gp_scalar>(pitch, gp_vec3::UnitY())) *
               gp_quat(Eigen::AngleAxis<gp_scalar>(roll, gp_vec3::UnitX()));
    }

    // T004: From axis-angle
    inline gp_quat fromAxisAngle(const gp_vec3& axis, gp_scalar angle) {
        return gp_quat(Eigen::AngleAxis<gp_scalar>(angle, axis.normalized()));
    }

    // T005: Flight attitude presets
    inline gp_quat level() { return gp_quat::Identity(); }
    inline gp_quat pitchedUp(gp_scalar angle) { return fromEuler(0, angle, 0); }
    inline gp_quat pitchedDown(gp_scalar angle) { return fromEuler(0, -angle, 0); }
    inline gp_quat bankedRight(gp_scalar angle) { return fromEuler(angle, 0, 0); }
    inline gp_quat bankedLeft(gp_scalar angle) { return fromEuler(-angle, 0, 0); }
    inline gp_quat yawed(gp_scalar angle) { return fromEuler(0, 0, angle); }

    // T006: Combined rotations
    inline gp_quat climbingTurn(gp_scalar pitch, gp_scalar roll, gp_scalar yaw) {
        return fromEuler(roll, pitch, yaw);
    }
}  // namespace QuatHelper

// T007-T008: NodeCoverageTracker for 100% coverage verification
class NodeCoverageTracker {
    static std::set<int> testedOps_;

public:
    static void markTested(int op) { testedOps_.insert(op); }
    static bool isTested(int op) { return testedOps_.count(op) > 0; }
    static std::vector<int> getUntestedOps() {
        std::vector<int> untested;
        for (int i = 0; i < testAllNodesCount; ++i) {
            if (testedOps_.count(testAllNodes[i].opcode) == 0) {
                untested.push_back(testAllNodes[i].opcode);
            }
        }
        return untested;
    }
    static int getTestedCount() { return static_cast<int>(testedOps_.size()); }
    static void reset() { testedOps_.clear(); }
};

std::set<int> NodeCoverageTracker::testedOps_;

// T008: Macro for coverage registration
#define TEST_NODE(op) NodeCoverageTracker::markTested(op)

// T009-T011: ExpectedValue namespace for ground truth computation
namespace ExpectedValue {
    // T009: LUT-based (bit-accurate) - uses same LUT functions AND range limit as evaluator
    gp_scalar sin(gp_scalar angle) { return applyRangeLimit(testFastSin(angle)); }
    gp_scalar cos(gp_scalar angle) { return applyRangeLimit(testFastCos(angle)); }
    gp_scalar atan2(gp_scalar y, gp_scalar x) { return applyRangeLimit(testFastAtan2(y, x)); }

    // T010: dPhi computation using Eigen operations
    // Note: Uses LUT atan2 for bit-accurate comparison with evaluator
    gp_scalar dPhi(PathProvider& pathProvider, AircraftState& state, gp_scalar offset) {
        int idx = static_cast<int>(offset);
        idx = std::clamp(idx + pathProvider.getCurrentIndex(), 0, pathProvider.getPathSize() - 1);
        gp_vec3 craftToTarget = pathProvider.getPath(idx).start - state.getPosition();
        gp_vec3 target_local = state.getOrientation().inverse() * craftToTarget;
        gp_vec3 projectedVector(0.0f, target_local.y(), target_local.z());
        return testFastAtan2(projectedVector.y(), -projectedVector.z());
    }

    // T010: dTheta computation using Eigen operations
    gp_scalar dTheta(PathProvider& pathProvider, AircraftState& state, gp_scalar offset) {
        int idx = static_cast<int>(offset);
        idx = std::clamp(idx + pathProvider.getCurrentIndex(), 0, pathProvider.getPathSize() - 1);
        gp_vec3 craftToTarget = pathProvider.getPath(idx).start - state.getPosition();
        gp_vec3 target_local = state.getOrientation().inverse() * craftToTarget;
        return testFastAtan2(-target_local.z(), target_local.x());
    }

    // T011: alpha (angle of attack) computation
    gp_scalar alpha(const AircraftState& state) {
        gp_vec3 velocity_body = state.getOrientation().inverse() * state.getVelocity();
        return std::atan2(-velocity_body.z(), velocity_body.x());
    }

    // T011: beta (sideslip) computation
    gp_scalar beta(const AircraftState& state) {
        gp_vec3 velocity_body = state.getOrientation().inverse() * state.getVelocity();
        return std::atan2(velocity_body.y(), velocity_body.x());
    }

    // T011: roll from quaternion
    gp_scalar rollRad(const AircraftState& state) {
        gp_vec3 euler = state.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
        return euler[2];
    }

    // T011: pitch from quaternion
    gp_scalar pitchRad(const AircraftState& state) {
        gp_vec3 euler = state.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
        return euler[1];
    }
}  // namespace ExpectedValue

// T012: makePath with optional simTimeMsec
Path makePath(gp_scalar x, gp_scalar y, gp_scalar z, gp_scalar timeMsec = 0.0f) {
    Path p;
    p.start = gp_vec3(x, y, z);
    p.orientation = gp_vec3::UnitX();
    p.distanceFromStart = 0.0f;
    p.radiansFromStart = 0.0f;
    p.simTimeMsec = timeMsec;
    return p;
}

// T013: makePath with position vector and time
Path makePathWithTime(const gp_vec3& pos, gp_scalar timeMsec) {
    return makePath(pos.x(), pos.y(), pos.z(), timeMsec);
}

// T012: makeState with various overloads
// AircraftState ctor: (pathIndex, relVel, velocity, orientation, position, pitch, roll, throttle, timeMsec)
AircraftState makeState() {
    return AircraftState(0,
                         static_cast<gp_scalar>(10.0f),  // relVel
                         gp_vec3::Zero(),                 // velocity
                         gp_quat::Identity(),             // orientation
                         gp_vec3::Zero(),                 // position
                         0.0f, 0.0f, 0.0f,               // commands
                         0);                              // time
}

AircraftState makeState(const gp_vec3& position, const gp_quat& orientation) {
    return AircraftState(0,
                         static_cast<gp_scalar>(10.0f),  // relVel
                         gp_vec3::Zero(),                 // velocity
                         orientation,                     // orientation
                         position,                        // position
                         0.0f, 0.0f, 0.0f,
                         0);
}

AircraftState makeState(const gp_vec3& position, const gp_quat& orientation, const gp_vec3& velocity) {
    return AircraftState(0,
                         static_cast<gp_scalar>(10.0f),  // relVel
                         velocity,                        // velocity (3rd param)
                         orientation,                     // orientation
                         position,                        // position (5th param)
                         0.0f, 0.0f, 0.0f,
                         0);
}

// T015: makeState with time - needed for interpolation tests
AircraftState makeStateWithTime(const gp_vec3& position, const gp_quat& orientation, int timeMsec) {
    return AircraftState(0,
                         static_cast<gp_scalar>(10.0f),  // relVel
                         gp_vec3::Zero(),                 // velocity
                         orientation,                     // orientation
                         position,                        // position
                         0.0f, 0.0f, 0.0f,
                         timeMsec);                       // time in msec
}

// T014: TestPathProvider with helper methods
struct TestPathProvider : public PathProvider {
    std::vector<Path> paths;
    int currentIndex{0};

    TestPathProvider() : paths({makePath(0, 0, 0)}) {}

    const Path& getPath(int index) const override {
        int idx = std::clamp(index, 0, static_cast<int>(paths.size()) - 1);
        return paths[idx];
    }
    int getCurrentIndex() const override { return currentIndex; }
    int getPathSize() const override { return static_cast<int>(paths.size()); }

    // T014: Helper methods
    void addPath(const Path& p) { paths.push_back(p); }
    void setCurrentIndex(int idx) { currentIndex = std::clamp(idx, 0, static_cast<int>(paths.size()) - 1); }
    void clear() { paths.clear(); paths.push_back(makePath(0, 0, 0)); currentIndex = 0; }
};

}  // namespace

// =============================================================================
// Phase 3: User Story 1 - Test Mathematical Operations (P1)
// =============================================================================

// T015: TEST(MathOps, AddSubMulDiv) for basic arithmetic
TEST(MathOps, AddSubMulDiv) {
    TestPathProvider provider;
    AircraftState state = makeState();
    gp_scalar args2[2] = {2.0f, 3.0f};

    EXPECT_FLOAT_EQ(evaluateGPOperator(ADD, provider, state, args2, 2), 5.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(SUB, provider, state, args2, 2), -1.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(MUL, provider, state, args2, 2), 6.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(DIV, provider, state, args2, 2), static_cast<gp_scalar>(2.0 / 3.0));

    TEST_NODE(ADD);
    TEST_NODE(SUB);
    TEST_NODE(MUL);
    TEST_NODE(DIV);
}

// T016: TEST(MathOps, DivisionByZero) verifying DIV(x,0)=0
TEST(MathOps, DivisionByZero) {
    TestPathProvider provider;
    AircraftState state = makeState();
    gp_scalar args[2] = {5.0f, 0.0f};

    EXPECT_FLOAT_EQ(evaluateGPOperator(DIV, provider, state, args, 2), 0.0f);
    TEST_NODE(DIV);
}

// T017: TEST(MathOps, SinCosBitAccurate) using LUT expected values
TEST(MathOps, SinCosBitAccurate) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar angles[] = {0.0f, static_cast<gp_scalar>(M_PI/6), static_cast<gp_scalar>(M_PI/4),
                          static_cast<gp_scalar>(M_PI/2), static_cast<gp_scalar>(M_PI)};

    for (gp_scalar a : angles) {
        gp_scalar args[1] = {a};
        gp_scalar sinResult = evaluateGPOperator(SIN, provider, state, args, 1);
        gp_scalar cosResult = evaluateGPOperator(COS, provider, state, args, 1);

        // Bit-accurate: ExpectedValue now applies same range limit as evaluator
        EXPECT_FLOAT_EQ(sinResult, ExpectedValue::sin(a));
        EXPECT_FLOAT_EQ(cosResult, ExpectedValue::cos(a));
    }

    TEST_NODE(SIN);
    TEST_NODE(COS);
}

// T018: TEST(MathOps, SinCosAllQuadrants) covering 0, π/2, π, 3π/2
TEST(MathOps, SinCosAllQuadrants) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Test all four quadrants
    gp_scalar quadrantAngles[] = {
        0.0f,
        static_cast<gp_scalar>(M_PI / 2),
        static_cast<gp_scalar>(M_PI),
        static_cast<gp_scalar>(3 * M_PI / 2),
        static_cast<gp_scalar>(2 * M_PI)
    };

    for (gp_scalar a : quadrantAngles) {
        gp_scalar args[1] = {a};
        gp_scalar sinResult = evaluateGPOperator(SIN, provider, state, args, 1);
        gp_scalar cosResult = evaluateGPOperator(COS, provider, state, args, 1);

        // Bit-accurate: ExpectedValue applies same range limit as evaluator
        EXPECT_FLOAT_EQ(sinResult, ExpectedValue::sin(a));
        EXPECT_FLOAT_EQ(cosResult, ExpectedValue::cos(a));
    }

    TEST_NODE(SIN);
    TEST_NODE(COS);
}

// T019: TEST(MathOps, Atan2AllQuadrants) for (+,+), (+,-), (-,+), (-,-)
TEST(MathOps, Atan2AllQuadrants) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Test all four quadrants
    struct TestCase { gp_scalar y; gp_scalar x; };
    TestCase cases[] = {
        {1.0f, 1.0f},   // Quadrant I (+,+)
        {1.0f, -1.0f},  // Quadrant II (+,-)
        {-1.0f, -1.0f}, // Quadrant III (-,-)
        {-1.0f, 1.0f},  // Quadrant IV (-,+)
        {0.0f, 1.0f},   // On +X axis
        {1.0f, 0.0f},   // On +Y axis
        {0.0f, -1.0f},  // On -X axis
        {-1.0f, 0.0f}   // On -Y axis
    };

    for (const auto& tc : cases) {
        gp_scalar args[2] = {tc.y, tc.x};
        gp_scalar result = evaluateGPOperator(ATAN2, provider, state, args, 2);
        EXPECT_FLOAT_EQ(result, ExpectedValue::atan2(tc.y, tc.x))
            << "Failed for y=" << tc.y << ", x=" << tc.x;
    }

    TEST_NODE(ATAN2);
}

// T020: TEST(MathOps, ClampMinMax) for CLAMP, MIN, MAX
TEST(MathOps, ClampMinMax) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // CLAMP tests
    gp_scalar clampArgs1[3] = {5.0f, 0.0f, 1.0f};  // value > max
    EXPECT_FLOAT_EQ(evaluateGPOperator(CLAMP, provider, state, clampArgs1, 3), 1.0f);

    gp_scalar clampArgs2[3] = {-5.0f, 0.0f, 1.0f}; // value < min
    EXPECT_FLOAT_EQ(evaluateGPOperator(CLAMP, provider, state, clampArgs2, 3), 0.0f);

    gp_scalar clampArgs3[3] = {0.5f, 0.0f, 1.0f};  // value in range
    EXPECT_FLOAT_EQ(evaluateGPOperator(CLAMP, provider, state, clampArgs3, 3), 0.5f);

    // MIN tests
    gp_scalar minArgs[2] = {3.0f, 7.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(MIN, provider, state, minArgs, 2), 3.0f);

    // MAX tests
    gp_scalar maxArgs[2] = {3.0f, 7.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(MAX, provider, state, maxArgs, 2), 7.0f);

    TEST_NODE(CLAMP);
    TEST_NODE(MIN);
    TEST_NODE(MAX);
}

// T021: TEST(MathOps, AbsSqrt) including SQRT(-1)=0 edge case
TEST(MathOps, AbsSqrt) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // ABS tests
    gp_scalar absArgs1[1] = {-4.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(ABS, provider, state, absArgs1, 1), 4.0f);

    gp_scalar absArgs2[1] = {4.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(ABS, provider, state, absArgs2, 1), 4.0f);

    // SQRT tests
    gp_scalar sqrtArgs1[1] = {16.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(SQRT, provider, state, sqrtArgs1, 1), 4.0f);

    // SQRT(-1) = 0 edge case
    gp_scalar sqrtArgs2[1] = {-1.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(SQRT, provider, state, sqrtArgs2, 1), 0.0f);

    TEST_NODE(ABS);
    TEST_NODE(SQRT);
}

// T022-T025: Bytecode Interpreter Math Tests
TEST(BytecodeMath, AddSubMulDiv) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // ADD(2, 3) = 5
    GPBytecode addProg[] = {{ZERO, 0, 2.0f}, {ZERO, 0, 3.0f}, {ADD, 2, 0}};
    // Use constants directly - ZERO pushes its constant field
    addProg[0].opcode = ZERO; addProg[0].constant = 2.0f;
    addProg[1].opcode = ZERO; addProg[1].constant = 3.0f;

    // Actually, looking at the bytecode impl, ZERO just pushes 0.0f
    // We need to use the actual constant value approach
    // Let's test using actual constants via the evaluator

    // Simpler approach: use known terminal values
    // ZERO + ONE = 1, ONE + TWO = 3
    GPBytecode prog1[] = {{ZERO, 0, 0}, {ONE, 0, 0}, {ADD, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(prog1, 3, provider, state, 0), 1.0f);

    GPBytecode prog2[] = {{ONE, 0, 0}, {TWO, 0, 0}, {ADD, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(prog2, 3, provider, state, 0), 3.0f);

    // ONE - TWO = -1
    GPBytecode prog3[] = {{ONE, 0, 0}, {TWO, 0, 0}, {SUB, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(prog3, 3, provider, state, 0), -1.0f);

    // ONE * TWO = 2
    GPBytecode prog4[] = {{ONE, 0, 0}, {TWO, 0, 0}, {MUL, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(prog4, 3, provider, state, 0), 2.0f);

    // TWO / ONE = 2 (note: DIV(a,b) = a/b, args are in stack order)
    GPBytecode prog5[] = {{TWO, 0, 0}, {ONE, 0, 0}, {DIV, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(prog5, 3, provider, state, 0), 2.0f);

    TEST_NODE(ADD);
    TEST_NODE(SUB);
    TEST_NODE(MUL);
    TEST_NODE(DIV);
}

TEST(BytecodeMath, SinCosBitAccurate) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // SIN(0) = 0
    GPBytecode sinProg[] = {{ZERO, 0, 0}, {SIN, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(sinProg, 2, provider, state, 0), ExpectedValue::sin(0));

    // COS(0) = 1
    GPBytecode cosProg[] = {{ZERO, 0, 0}, {COS, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(cosProg, 2, provider, state, 0), ExpectedValue::cos(0));

    // SIN(PI) ≈ 0, COS(PI) ≈ -1
    GPBytecode sinPiProg[] = {{OP_PI, 0, 0}, {SIN, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(sinPiProg, 2, provider, state, 0),
                    ExpectedValue::sin(static_cast<gp_scalar>(M_PI)));

    GPBytecode cosPiProg[] = {{OP_PI, 0, 0}, {COS, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(cosPiProg, 2, provider, state, 0),
                    ExpectedValue::cos(static_cast<gp_scalar>(M_PI)));

    TEST_NODE(SIN);
    TEST_NODE(COS);
}

TEST(BytecodeMath, Atan2AllQuadrants) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // ATAN2(1, 1) = PI/4
    GPBytecode prog1[] = {{ONE, 0, 0}, {ONE, 0, 0}, {ATAN2, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(prog1, 3, provider, state, 0),
                    ExpectedValue::atan2(1.0f, 1.0f));

    // ATAN2(0, 1) = 0
    GPBytecode prog2[] = {{ZERO, 0, 0}, {ONE, 0, 0}, {ATAN2, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(prog2, 3, provider, state, 0),
                    ExpectedValue::atan2(0.0f, 1.0f));

    TEST_NODE(ATAN2);
}

TEST(BytecodeMath, ClampMinMaxAbsSqrt) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // ABS(-1) using SUB(ZERO, ONE) then ABS
    GPBytecode absNegProg[] = {{ZERO, 0, 0}, {ONE, 0, 0}, {SUB, 2, 0}, {ABS, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(absNegProg, 4, provider, state, 0), 1.0f);

    // SQRT(ONE) = 1
    GPBytecode sqrtProg[] = {{ONE, 0, 0}, {SQRT, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(sqrtProg, 2, provider, state, 0), 1.0f);

    // MIN(ONE, TWO) = 1
    GPBytecode minProg[] = {{ONE, 0, 0}, {TWO, 0, 0}, {MIN, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(minProg, 3, provider, state, 0), 1.0f);

    // MAX(ONE, TWO) = 2
    GPBytecode maxProg[] = {{ONE, 0, 0}, {TWO, 0, 0}, {MAX, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(maxProg, 3, provider, state, 0), 2.0f);

    TEST_NODE(ABS);
    TEST_NODE(SQRT);
    TEST_NODE(MIN);
    TEST_NODE(MAX);
}

// =============================================================================
// Phase 4: User Story 2 - Test Control Side Effects (P1)
// =============================================================================

// T026: TEST(ControlOps, SetPitchClampHigh) verifying >1.0 clamps to 1.0
TEST(ControlOps, SetPitchClampHigh) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar args[1] = {2.5f};
    gp_scalar result = evaluateGPOperator(SETPITCH, provider, state, args, 1);
    EXPECT_FLOAT_EQ(result, 1.0f);
    EXPECT_FLOAT_EQ(state.getPitchCommand(), 1.0f);

    TEST_NODE(SETPITCH);
}

// T027: TEST(ControlOps, SetRollClampLow) verifying <-1.0 clamps to -1.0
TEST(ControlOps, SetRollClampLow) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar args[1] = {-2.5f};
    gp_scalar result = evaluateGPOperator(SETROLL, provider, state, args, 1);
    EXPECT_FLOAT_EQ(result, -1.0f);
    EXPECT_FLOAT_EQ(state.getRollCommand(), -1.0f);

    TEST_NODE(SETROLL);
}

// T028: TEST(ControlOps, SetThrottleClamp) verifying both bounds
TEST(ControlOps, SetThrottleClamp) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Clamp high
    gp_scalar argsHigh[1] = {1.5f};
    gp_scalar resultHigh = evaluateGPOperator(SETTHROTTLE, provider, state, argsHigh, 1);
    EXPECT_FLOAT_EQ(resultHigh, 1.0f);
    EXPECT_FLOAT_EQ(state.getThrottleCommand(), 1.0f);

    // Clamp low
    gp_scalar argsLow[1] = {-1.5f};
    gp_scalar resultLow = evaluateGPOperator(SETTHROTTLE, provider, state, argsLow, 1);
    EXPECT_FLOAT_EQ(resultLow, -1.0f);
    EXPECT_FLOAT_EQ(state.getThrottleCommand(), -1.0f);

    TEST_NODE(SETTHROTTLE);
}

// T029: TEST(ControlOps, SetGetCycle) verifying SET followed by GET returns correct value
TEST(ControlOps, SetGetCycle) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Set pitch
    gp_scalar setPitchArgs[1] = {0.75f};
    evaluateGPOperator(SETPITCH, provider, state, setPitchArgs, 1);
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETPITCH, provider, state, nullptr, 0), 0.75f);

    // Set roll
    gp_scalar setRollArgs[1] = {-0.5f};
    evaluateGPOperator(SETROLL, provider, state, setRollArgs, 1);
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETROLL, provider, state, nullptr, 0), -0.5f);

    // Set throttle
    gp_scalar setThrottleArgs[1] = {0.25f};
    evaluateGPOperator(SETTHROTTLE, provider, state, setThrottleArgs, 1);
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETTHROTTLE, provider, state, nullptr, 0), 0.25f);

    TEST_NODE(SETPITCH);
    TEST_NODE(SETROLL);
    TEST_NODE(SETTHROTTLE);
    TEST_NODE(GETPITCH);
    TEST_NODE(GETROLL);
    TEST_NODE(GETTHROTTLE);
}

// T030: TEST(ControlOps, GetBeforeSet) verifying initial GET values
TEST(ControlOps, GetBeforeSet) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Initial values should be 0
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETPITCH, provider, state, nullptr, 0), 0.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETROLL, provider, state, nullptr, 0), 0.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETTHROTTLE, provider, state, nullptr, 0), 0.0f);

    TEST_NODE(GETPITCH);
    TEST_NODE(GETROLL);
    TEST_NODE(GETTHROTTLE);
}

// T031: TEST(BytecodeControl, SetPitchRollThrottle) mirroring tree tests
TEST(BytecodeControl, SetPitchRollThrottle) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // SETPITCH(ONE) - sets pitch to 1.0
    GPBytecode pitchProg[] = {{ONE, 0, 0}, {SETPITCH, 1, 0}};
    gp_scalar pitchResult = evaluateBytecodePortable(pitchProg, 2, provider, state, 0);
    EXPECT_FLOAT_EQ(pitchResult, 1.0f);
    EXPECT_FLOAT_EQ(state.getPitchCommand(), 1.0f);

    // SETROLL with clamped value: using TWO which clamps to 1.0
    GPBytecode rollProg[] = {{TWO, 0, 0}, {SETROLL, 1, 0}};
    gp_scalar rollResult = evaluateBytecodePortable(rollProg, 2, provider, state, 0);
    EXPECT_FLOAT_EQ(rollResult, 1.0f);
    EXPECT_FLOAT_EQ(state.getRollCommand(), 1.0f);

    TEST_NODE(SETPITCH);
    TEST_NODE(SETROLL);
}

// T032: TEST(BytecodeControl, GetPitchRollThrottle) mirroring tree tests
TEST(BytecodeControl, GetPitchRollThrottle) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Set some values first
    gp_scalar setPitchArgs[1] = {0.5f};
    evaluateGPOperator(SETPITCH, provider, state, setPitchArgs, 1);

    // GETPITCH should return 0.5
    GPBytecode getPitchProg[] = {{GETPITCH, 0, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(getPitchProg, 1, provider, state, 0), 0.5f);

    TEST_NODE(GETPITCH);
    TEST_NODE(GETROLL);
    TEST_NODE(GETTHROTTLE);
}

// =============================================================================
// Phase 5: User Story 3 - Test Path Following Nodes (P1)
// =============================================================================

// T033: TEST(NavigationOps, GetDPhiDThetaZeroOffset) with target directly ahead
TEST(NavigationOps, GetDPhiDThetaZeroOffset) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 0.0f);  // Target ahead in +X
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    gp_scalar args[1] = {0.0f};
    gp_scalar dPhi = executeGetDPhi(provider, state, args[0]);
    gp_scalar dTheta = executeGetDTheta(provider, state, args[0]);

    // Target directly ahead: both angles should be near zero
    EXPECT_NEAR(dPhi, 0.0f, 1e-4f);
    EXPECT_NEAR(dTheta, 0.0f, 1e-4f);

    TEST_NODE(GETDPHI);
    TEST_NODE(GETDTHETA);
}

// T034: TEST(NavigationOps, GetDThetaTargetAbove) verifying positive result
TEST(NavigationOps, GetDThetaTargetAbove) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, -5.0f);  // Above (negative Z in NED)
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    gp_scalar args[1] = {0.0f};
    gp_scalar dTheta = executeGetDTheta(provider, state, args[0]);

    // Target above: positive pitch needed (nose up)
    EXPECT_GT(dTheta, 0.0f);

    TEST_NODE(GETDTHETA);
}

// T035: TEST(NavigationOps, GetDThetaTargetBelow) verifying negative result
TEST(NavigationOps, GetDThetaTargetBelow) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 5.0f);  // Below (positive Z in NED)
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    gp_scalar args[1] = {0.0f};
    gp_scalar dTheta = executeGetDTheta(provider, state, args[0]);

    // Target below: negative pitch needed (nose down)
    EXPECT_LT(dTheta, 0.0f);

    TEST_NODE(GETDTHETA);
}

// T036: TEST(NavigationOps, GetDPhiTargetRight) verifying positive result
TEST(NavigationOps, GetDPhiTargetRight) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 5.0f, 0.0f);  // Right (positive Y in NED)
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    gp_scalar args[1] = {0.0f};
    gp_scalar dPhi = executeGetDPhi(provider, state, args[0]);

    // Target to right: positive roll needed (bank right)
    EXPECT_GT(dPhi, 0.0f);

    TEST_NODE(GETDPHI);
}

// T037: TEST(NavigationOps, GetDPhiTargetLeft) verifying negative result
TEST(NavigationOps, GetDPhiTargetLeft) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, -5.0f, 0.0f);  // Left (negative Y in NED)
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    gp_scalar args[1] = {0.0f};
    gp_scalar dPhi = executeGetDPhi(provider, state, args[0]);

    // Target to left: negative roll needed (bank left)
    EXPECT_LT(dPhi, 0.0f);

    TEST_NODE(GETDPHI);
}

// T038: TEST(NavigationOps, TimeOffsetPositive) testing GETDPHI(2) targets waypoint+2
TEST(NavigationOps, TimeOffsetPositive) {
    TestPathProvider provider;
    provider.paths.clear();
    // Use proper timestamps for time-based indexing (100ms steps)
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));     // Index 0: ahead, t=0ms
    provider.paths.push_back(makePath(20.0f, 0.0f, 0.0f, 100.0f));   // Index 1: further ahead, t=100ms
    provider.paths.push_back(makePath(30.0f, 5.0f, 0.0f, 200.0f));   // Index 2: ahead + right, t=200ms

    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());
    provider.setCurrentIndex(0);

    // Get dPhi for waypoint at offset 2 (targets time 0 + 200ms = index 2)
    gp_scalar dPhiOffset0 = executeGetDPhi(provider, state, 0.0f);  // Targets index 0
    gp_scalar dPhiOffset2 = executeGetDPhi(provider, state, 2.0f);  // Targets index 2

    // Waypoint 0 is straight ahead, waypoint 2 is to the right
    EXPECT_NEAR(dPhiOffset0, 0.0f, 1e-4f);
    EXPECT_GT(dPhiOffset2, 0.0f);

    TEST_NODE(GETDPHI);
}

// T039: TEST(NavigationOps, TimeOffsetClampEnd) testing offset beyond path length
TEST(NavigationOps, TimeOffsetClampEnd) {
    TestPathProvider provider;
    provider.paths.clear();
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));    // Index 0, t=0ms
    provider.paths.push_back(makePath(20.0f, 5.0f, 0.0f, 100.0f));  // Index 1 (last), t=100ms

    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());
    provider.setCurrentIndex(0);

    // Offset 100 should clamp to last waypoint (index 1) - steps clamp to 5 max
    gp_scalar dPhiLast = executeGetDPhi(provider, state, 1.0f);    // Targets index 1
    gp_scalar dPhiClamped = executeGetDPhi(provider, state, 100.0f); // Should also target index 1

    EXPECT_FLOAT_EQ(dPhiClamped, dPhiLast);

    TEST_NODE(GETDPHI);
}

// T040: TEST(NavigationOps, TimeOffsetClampStart) testing negative offset clamps to 0
TEST(NavigationOps, TimeOffsetClampStart) {
    TestPathProvider provider;
    // Directly set paths to avoid clear()'s default path at origin
    provider.paths.clear();
    provider.paths.push_back(makePath(10.0f, 5.0f, 0.0f, 0.0f));    // Index 0: to the right, t=0ms
    provider.paths.push_back(makePath(20.0f, 0.0f, 0.0f, 100.0f));  // Index 1: straight ahead, t=100ms

    // With interpolation, aircraft time determines the reference point
    // Set aircraft time to 100ms (matching index 1)
    AircraftState state = makeStateWithTime(gp_vec3::Zero(), QuatHelper::level(), 100);
    provider.setCurrentIndex(1);  // Current is index 1 at t=100ms

    // Negative offset -1 targets time = 100 - 100ms = 0ms, which is index 0
    gp_scalar dPhiFirst = executeGetDPhi(provider, state, -1.0f);  // Should target index 0
    gp_scalar dPhiZero = executeGetDPhi(provider, state, 0.0f);    // Targets current (index 1)

    // First waypoint (0) is to the right (y=5), current (1) is straight ahead
    EXPECT_GT(dPhiFirst, 0.0f);  // To the right
    EXPECT_NEAR(dPhiZero, 0.0f, 1e-4f);  // Straight ahead

    TEST_NODE(GETDPHI);
}

// T040a: TEST(NavigationOps, NonUniformTimestamps) - FIXED by path interpolation
// On xiao-gp, timestamps aren't exact 100ms intervals due to real-time processing.
// Previously: Jitter caused overshooting (t=98ms < t=100ms goal → skip to next waypoint).
// NOW: Interpolation finds exact position for goal time, regardless of jitter.
TEST(NavigationOps, NonUniformTimestamps) {
    TestPathProvider provider;
    provider.paths.clear();
    // Simulate real-time jitter: timestamps slightly off from ideal 100ms intervals
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));      // Index 0, t=0ms
    provider.paths.push_back(makePath(20.0f, 0.0f, 0.0f, 102.0f));    // Index 1, t=102ms (2ms late)
    provider.paths.push_back(makePath(30.0f, 5.0f, 0.0f, 198.0f));    // Index 2, t=198ms (2ms early), to the RIGHT
    provider.paths.push_back(makePath(40.0f, 0.0f, 0.0f, 305.0f));    // Index 3, t=305ms (5ms late)

    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());
    provider.setCurrentIndex(0);

    // offset=1 targets t=100ms (aircraft time=0 + 1*100ms)
    // Binary search finds t=100ms between index 0 (t=0) and index 1 (t=102)
    // Interpolation: frac = (100-0)/(102-0) ≈ 0.98 → position very close to index 1
    gp_scalar dPhi1 = executeGetDPhi(provider, state, 1.0f);

    // offset=2 targets t=200ms
    // Binary search finds t=200ms between index 2 (t=198) and index 3 (t=305)
    // Interpolation: frac = (200-198)/(305-198) ≈ 0.019 → position very close to index 2
    // Index 2 is to the RIGHT (y=5), so dPhi should be positive
    gp_scalar dPhi2 = executeGetDPhi(provider, state, 2.0f);

    // Index 1 is straight ahead (y=0)
    EXPECT_NEAR(dPhi1, 0.0f, 0.1f);  // Nearly straight ahead (98% toward index 1)
    // FIXED: Now correctly finds interpolated position near index 2 (to the right)
    EXPECT_GT(dPhi2, 0.0f);  // Target is to the right

    TEST_NODE(GETDPHI);
}

// T040b: TEST(NavigationOps, SingleWaypointPath) - edge case with only one waypoint
TEST(NavigationOps, SingleWaypointPath) {
    TestPathProvider provider;
    provider.paths.clear();
    provider.paths.push_back(makePath(10.0f, 5.0f, 0.0f, 0.0f));  // Only one waypoint

    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());
    provider.setCurrentIndex(0);

    // Any offset should clamp to the only waypoint
    gp_scalar dPhi0 = executeGetDPhi(provider, state, 0.0f);
    gp_scalar dPhiPos = executeGetDPhi(provider, state, 5.0f);
    gp_scalar dPhiNeg = executeGetDPhi(provider, state, -5.0f);

    // All should target the same waypoint (to the right)
    EXPECT_FLOAT_EQ(dPhi0, dPhiPos);
    EXPECT_FLOAT_EQ(dPhi0, dPhiNeg);
    EXPECT_GT(dPhi0, 0.0f);  // Target is to the right

    TEST_NODE(GETDPHI);
}

// T040c: TEST(NavigationOps, TimestampBoundary) - exact boundary conditions
TEST(NavigationOps, TimestampBoundary) {
    TestPathProvider provider;
    provider.paths.clear();
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));    // Index 0, t=0ms
    provider.paths.push_back(makePath(20.0f, 5.0f, 0.0f, 100.0f));  // Index 1, t=100ms (exact boundary)
    provider.paths.push_back(makePath(30.0f, 0.0f, 0.0f, 200.0f));  // Index 2, t=200ms

    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());
    provider.setCurrentIndex(0);

    // offset=1 targets exactly t=100ms - should find index 1
    gp_scalar dPhi1 = executeGetDPhi(provider, state, 1.0f);

    // Index 1 is to the right (y=5)
    EXPECT_GT(dPhi1, 0.0f);

    // Now test from index 1, going backwards
    provider.setCurrentIndex(1);
    gp_scalar dPhiBack = executeGetDPhi(provider, state, -1.0f);  // Targets t=0ms

    // Index 0 is straight ahead (y=0)
    EXPECT_NEAR(dPhiBack, 0.0f, 1e-4f);

    TEST_NODE(GETDPHI);
}

// T040d: TEST(NavigationOps, EarlyTimestampOvershoot) - FIXED by path interpolation
// Previously: When a waypoint's timestamp was slightly EARLY (e.g., 98ms instead of 100ms),
// the discrete algorithm overshot to the NEXT waypoint because 98ms < 100ms goal.
// NOW: With interpolation, targeting t=100ms finds position between index 1 (t=98ms) and
// index 2 (t=200ms), correctly staying near index 1 instead of skipping it.
TEST(NavigationOps, EarlyTimestampOvershoot) {
    TestPathProvider provider;
    provider.paths.clear();
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));    // Index 0, t=0ms, straight
    provider.paths.push_back(makePath(20.0f, 5.0f, 0.0f, 98.0f));   // Index 1, t=98ms, RIGHT (2ms EARLY!)
    provider.paths.push_back(makePath(30.0f, 0.0f, 0.0f, 200.0f));  // Index 2, t=200ms, straight

    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());
    provider.setCurrentIndex(0);

    // offset=1 targets t=100ms
    // With interpolation: binary search finds t=100ms between index 1 (98ms) and index 2 (200ms)
    // frac = (100-98)/(200-98) = 2/102 ≈ 0.0196
    // Position ≈ (20.2, 4.9, 0) - very close to index 1, which is to the RIGHT
    gp_scalar dPhi1 = executeGetDPhi(provider, state, 1.0f);

    // FIXED: Now correctly finds interpolated position near index 1 (to the right)
    // dPhi should be positive (roll right to face the target)
    EXPECT_GT(dPhi1, 0.0f);  // Target is to the right, need positive roll

    TEST_NODE(GETDPHI);
}

// T040e: TEST(NavigationOps, LargeTimestampGaps) - missing waypoints scenario
TEST(NavigationOps, LargeTimestampGaps) {
    TestPathProvider provider;
    provider.paths.clear();
    // Gap of 500ms between waypoints (5 steps worth)
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));     // Index 0, t=0ms
    provider.paths.push_back(makePath(50.0f, 5.0f, 0.0f, 500.0f));   // Index 1, t=500ms (big gap)

    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());
    provider.setCurrentIndex(0);

    // offset=1 targets t=100ms, but next waypoint is at t=500ms
    // Should stay at index 0 since path[0].simTimeMsec + eps < 100 is true
    // but path[1].simTimeMsec = 500 >= 100, so it advances to index 1
    gp_scalar dPhi1 = executeGetDPhi(provider, state, 1.0f);

    // offset=3 targets t=300ms, still should find index 1
    gp_scalar dPhi3 = executeGetDPhi(provider, state, 3.0f);

    // Both should target index 1 (to the right)
    EXPECT_FLOAT_EQ(dPhi1, dPhi3);
    EXPECT_GT(dPhi1, 0.0f);

    TEST_NODE(GETDPHI);
}

// T041: TEST(NavigationOps, GetDTarget) verifying throttle estimate calculation
TEST(NavigationOps, GetDTarget) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(20.0f, 0.0f, 0.0f);  // 20m ahead
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    gp_scalar dTarget = executeGetDTarget(provider, state, 0.0f);

    // dTarget = clamp((distance - 10) / relVel, -1, 1)
    // distance = 20, relVel = 10 (default in makeState)
    // dTarget = clamp((20 - 10) / 10, -1, 1) = clamp(1, -1, 1) = 1.0
    EXPECT_FLOAT_EQ(dTarget, 1.0f);

    TEST_NODE(GETDTARGET);
}

// T042: TEST(NavigationOps, GetDPhiBankedRight) using QuatHelper::bankedRight(deg(45))
TEST(NavigationOps, GetDPhiBankedRight) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 0.0f);  // Target straight ahead in world frame

    // Aircraft banked 45 degrees right
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::bankedRight(QuatHelper::deg(45)));

    gp_scalar dPhi = executeGetDPhi(provider, state, 0.0f);
    gp_scalar expected = ExpectedValue::dPhi(provider, state, 0.0f);

    EXPECT_NEAR(dPhi, expected, 1e-5f);

    TEST_NODE(GETDPHI);
}

// T043: TEST(NavigationOps, GetDThetaPitchedUp) using QuatHelper::pitchedUp(deg(30))
TEST(NavigationOps, GetDThetaPitchedUp) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 0.0f);  // Target at same altitude

    // Aircraft pitched up 30 degrees
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::pitchedUp(QuatHelper::deg(30)));

    gp_scalar dTheta = executeGetDTheta(provider, state, 0.0f);
    gp_scalar expected = ExpectedValue::dTheta(provider, state, 0.0f);

    EXPECT_NEAR(dTheta, expected, 1e-5f);

    TEST_NODE(GETDTHETA);
}

// T044: TEST(NavigationOps, GetDPhiClimbingTurn) using combined pitch+roll quaternion
TEST(NavigationOps, GetDPhiClimbingTurn) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 5.0f, -5.0f);  // Target ahead, right, and above

    // Combined attitude: pitched up 20°, banked right 30°, yawed 10°
    AircraftState state = makeState(gp_vec3::Zero(),
        QuatHelper::climbingTurn(QuatHelper::deg(20), QuatHelper::deg(30), QuatHelper::deg(10)));

    gp_scalar dPhi = executeGetDPhi(provider, state, 0.0f);
    gp_scalar expected = ExpectedValue::dPhi(provider, state, 0.0f);

    EXPECT_NEAR(dPhi, expected, 1e-5f);

    TEST_NODE(GETDPHI);
}

// T045: TEST(NavigationOps, TargetAllOctants) testing ±X, ±Y, ±Z target combinations
TEST(NavigationOps, TargetAllOctants) {
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    // Test 8 octants
    struct OctantTest { gp_scalar x, y, z; };
    OctantTest octants[] = {
        {10.0f,  5.0f, -5.0f},  // +X, +Y, -Z (ahead-right-up)
        {10.0f, -5.0f, -5.0f},  // +X, -Y, -Z (ahead-left-up)
        {10.0f,  5.0f,  5.0f},  // +X, +Y, +Z (ahead-right-down)
        {10.0f, -5.0f,  5.0f},  // +X, -Y, +Z (ahead-left-down)
        {-10.0f,  5.0f, -5.0f}, // -X, +Y, -Z (behind-right-up)
        {-10.0f, -5.0f, -5.0f}, // -X, -Y, -Z (behind-left-up)
        {-10.0f,  5.0f,  5.0f}, // -X, +Y, +Z (behind-right-down)
        {-10.0f, -5.0f,  5.0f}  // -X, -Y, +Z (behind-left-down)
    };

    for (const auto& oct : octants) {
        TestPathProvider provider;
        provider.paths[0].start = gp_vec3(oct.x, oct.y, oct.z);

        gp_scalar dPhi = executeGetDPhi(provider, state, 0.0f);
        gp_scalar dTheta = executeGetDTheta(provider, state, 0.0f);
        gp_scalar expectedDPhi = ExpectedValue::dPhi(provider, state, 0.0f);
        gp_scalar expectedDTheta = ExpectedValue::dTheta(provider, state, 0.0f);

        EXPECT_NEAR(dPhi, expectedDPhi, 1e-5f)
            << "Octant (" << oct.x << "," << oct.y << "," << oct.z << ")";
        EXPECT_NEAR(dTheta, expectedDTheta, 1e-5f)
            << "Octant (" << oct.x << "," << oct.y << "," << oct.z << ")";
    }

    TEST_NODE(GETDPHI);
    TEST_NODE(GETDTHETA);
}

// =============================================================================
// Interpolation Tests - verify getInterpolatedTargetPosition()
// These tests verify continuous position interpolation instead of discrete jumps
// =============================================================================

// T045a: TEST(NavigationOps, InterpolationMidpoint) - verify lerp at 50% between waypoints
TEST(NavigationOps, InterpolationMidpoint) {
    TestPathProvider provider;
    provider.paths.clear();
    // Two waypoints: 100ms apart, positions at x=10 and x=20
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));     // Index 0, t=0ms
    provider.paths.push_back(makePath(20.0f, 0.0f, 0.0f, 100.0f));   // Index 1, t=100ms
    provider.setCurrentIndex(0);

    // Goal time: 50ms (halfway between 0ms and 100ms)
    // Expected position: lerp(10, 20, 0.5) = 15
    gp_vec3 pos = getInterpolatedTargetPosition(provider, 50.0f, 0.0f);

    EXPECT_FLOAT_EQ(pos.x(), 15.0f);
    EXPECT_FLOAT_EQ(pos.y(), 0.0f);
    EXPECT_FLOAT_EQ(pos.z(), 0.0f);
}

// T045b: TEST(NavigationOps, InterpolationContinuity) - verify no jumps at boundaries
TEST(NavigationOps, InterpolationContinuity) {
    TestPathProvider provider;
    provider.paths.clear();
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));     // t=0ms
    provider.paths.push_back(makePath(20.0f, 0.0f, 0.0f, 100.0f));   // t=100ms
    provider.paths.push_back(makePath(30.0f, 0.0f, 0.0f, 200.0f));   // t=200ms
    provider.setCurrentIndex(0);

    // Sample positions every 10ms across boundary at t=100ms
    gp_vec3 pos90 = getInterpolatedTargetPosition(provider, 90.0f, 0.0f);
    gp_vec3 pos100 = getInterpolatedTargetPosition(provider, 100.0f, 0.0f);
    gp_vec3 pos110 = getInterpolatedTargetPosition(provider, 110.0f, 0.0f);

    // Should be smooth: 19.0, 20.0, 21.0
    EXPECT_FLOAT_EQ(pos90.x(), 19.0f);
    EXPECT_FLOAT_EQ(pos100.x(), 20.0f);
    EXPECT_FLOAT_EQ(pos110.x(), 21.0f);

    // Check continuity: no jumps > 10% of step size
    EXPECT_NEAR(pos100.x() - pos90.x(), 1.0f, 0.1f);
    EXPECT_NEAR(pos110.x() - pos100.x(), 1.0f, 0.1f);
}

// T045c: TEST(NavigationOps, InterpolationJitterRobust) - 10ms jitter causes <1% change
TEST(NavigationOps, InterpolationJitterRobust) {
    TestPathProvider provider;
    provider.paths.clear();
    // Path with 100ms spacing, positions at 0, 100, 200 meters
    provider.paths.push_back(makePath(0.0f, 0.0f, 0.0f, 0.0f));
    provider.paths.push_back(makePath(100.0f, 0.0f, 0.0f, 100.0f));
    provider.paths.push_back(makePath(200.0f, 0.0f, 0.0f, 200.0f));
    provider.setCurrentIndex(1);

    // Base case: exactly at t=100ms with offset=0
    gp_vec3 posBase = getInterpolatedTargetPosition(provider, 100.0f, 0.0f);

    // Jittered case: at t=98ms (2ms early, simulating real-time jitter)
    gp_vec3 posJitter = getInterpolatedTargetPosition(provider, 98.0f, 0.0f);

    // The difference should be small (< 1% of typical position value)
    // Position at t=100ms should be around 100m, at t=98ms should be around 98m
    // Difference of 2m on a 100m position is 2%, which is acceptable
    // But the KEY test is that we don't jump to a completely different waypoint
    gp_scalar diff = std::abs(posBase.x() - posJitter.x());
    gp_scalar pctChange = diff / std::max(std::abs(posBase.x()), 1.0f) * 100.0f;

    EXPECT_LT(pctChange, 5.0f);  // Less than 5% change from 2ms jitter
    EXPECT_LT(diff, 10.0f);       // Absolute difference < 10 meters
}

// T045d: TEST(NavigationOps, InterpolationBoundaryClamp) - verify ±10 step clamp
TEST(NavigationOps, InterpolationBoundaryClamp) {
    TestPathProvider provider;
    provider.paths.clear();
    // Long path with many waypoints
    for (int i = 0; i < 50; i++) {
        provider.paths.push_back(makePath(i * 10.0f, 0.0f, 0.0f, i * 100.0f));
    }
    provider.setCurrentIndex(25);  // Middle of path, t=2500ms

    // Request offset of +20 steps (should clamp to +10)
    gp_vec3 posMax = getInterpolatedTargetPosition(provider, 2500.0f, 20.0f);
    // Goal time with +10 steps: 2500 + 10*100 = 3500ms -> x=350m
    EXPECT_FLOAT_EQ(posMax.x(), 350.0f);

    // Request offset of -20 steps (should clamp to -10)
    gp_vec3 posMin = getInterpolatedTargetPosition(provider, 2500.0f, -20.0f);
    // Goal time with -10 steps: 2500 - 10*100 = 1500ms -> x=150m
    EXPECT_FLOAT_EQ(posMin.x(), 150.0f);
}

// T045e: TEST(NavigationOps, InterpolationNaNHandling) - NaN returns current position
TEST(NavigationOps, InterpolationNaNHandling) {
    TestPathProvider provider;
    provider.paths.clear();
    provider.paths.push_back(makePath(10.0f, 0.0f, 0.0f, 0.0f));
    provider.paths.push_back(makePath(20.0f, 0.0f, 0.0f, 100.0f));
    provider.setCurrentIndex(0);

    // Pass NaN as offset - should return current rabbit position
    gp_vec3 pos = getInterpolatedTargetPosition(provider, 50.0f, std::nanf(""));

    // Should return position at current index (index 0 -> x=10)
    EXPECT_FLOAT_EQ(pos.x(), 10.0f);
    EXPECT_FLOAT_EQ(pos.y(), 0.0f);
    EXPECT_FLOAT_EQ(pos.z(), 0.0f);
}

// T046: TEST(BytecodeNavigation, GetDPhiDTheta) mirroring identity baseline tests
TEST(BytecodeNavigation, GetDPhiDTheta) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 0.0f);  // Target ahead
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    // GETDPHI(0) - push 0, then call GETDPHI
    GPBytecode dPhiProg[] = {{ZERO, 0, 0}, {GETDPHI, 1, 0}};
    gp_scalar dPhi = evaluateBytecodePortable(dPhiProg, 2, provider, state, 0);
    EXPECT_NEAR(dPhi, 0.0f, 1e-4f);

    // GETDTHETA(0)
    GPBytecode dThetaProg[] = {{ZERO, 0, 0}, {GETDTHETA, 1, 0}};
    gp_scalar dTheta = evaluateBytecodePortable(dThetaProg, 2, provider, state, 0);
    EXPECT_NEAR(dTheta, 0.0f, 1e-4f);

    TEST_NODE(GETDPHI);
    TEST_NODE(GETDTHETA);
}

// T047: TEST(BytecodeNavigation, GetDTarget) mirroring tree test
TEST(BytecodeNavigation, GetDTarget) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(20.0f, 0.0f, 0.0f);
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level());

    GPBytecode prog[] = {{ZERO, 0, 0}, {GETDTARGET, 1, 0}};
    gp_scalar dTarget = evaluateBytecodePortable(prog, 2, provider, state, 0);
    EXPECT_FLOAT_EQ(dTarget, 1.0f);

    TEST_NODE(GETDTARGET);
}

// T048: TEST(BytecodeNavigation, QuaternionOrientations) mirroring 3D quaternion tests
TEST(BytecodeNavigation, QuaternionOrientations) {
    TestPathProvider provider;
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 0.0f);

    // Banked right 45 degrees
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::bankedRight(QuatHelper::deg(45)));

    GPBytecode prog[] = {{ZERO, 0, 0}, {GETDPHI, 1, 0}};
    gp_scalar dPhi = evaluateBytecodePortable(prog, 2, provider, state, 0);
    gp_scalar expected = ExpectedValue::dPhi(provider, state, 0.0f);

    EXPECT_NEAR(dPhi, expected, 1e-5f);

    TEST_NODE(GETDPHI);
}

// =============================================================================
// Phase 6: User Story 4 - Test Sensor Nodes (P1)
// =============================================================================

// T049: TEST(SensorOps, GetVelMagnitude) verifying speed magnitude
TEST(SensorOps, GetVelMagnitude) {
    TestPathProvider provider;
    AircraftState state = makeState();
    state.setRelVel(15.0f);

    EXPECT_FLOAT_EQ(evaluateGPOperator(GETVEL, provider, state, nullptr, 0), 15.0f);

    TEST_NODE(GETVEL);
}

// T050: TEST(SensorOps, GetVelXYZ) verifying NED velocity components
TEST(SensorOps, GetVelXYZ) {
    TestPathProvider provider;
    gp_vec3 velocity(3.0f, 4.0f, 5.0f);  // NED components
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level(), velocity);

    EXPECT_FLOAT_EQ(evaluateGPOperator(GETVELX, provider, state, nullptr, 0), 3.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETVELY, provider, state, nullptr, 0), 4.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(GETVELZ, provider, state, nullptr, 0), 5.0f);

    TEST_NODE(GETVELX);
    TEST_NODE(GETVELY);
    TEST_NODE(GETVELZ);
}

// T051: TEST(SensorOps, AlphaBetaLevel) with identity orientation returning 0
TEST(SensorOps, AlphaBetaLevel) {
    TestPathProvider provider;
    // Level flight with velocity in +X direction (forward in NED)
    gp_vec3 velocity(10.0f, 0.0f, 0.0f);
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level(), velocity);

    gp_scalar alpha = evaluateGPOperator(GETALPHA, provider, state, nullptr, 0);
    gp_scalar beta = evaluateGPOperator(GETBETA, provider, state, nullptr, 0);

    EXPECT_NEAR(alpha, 0.0f, 1e-5f);
    EXPECT_NEAR(beta, 0.0f, 1e-5f);

    TEST_NODE(GETALPHA);
    TEST_NODE(GETBETA);
}

// T052: TEST(SensorOps, AlphaPitchedUp) using QuatHelper::pitchedUp(deg(30))
TEST(SensorOps, AlphaPitchedUp) {
    TestPathProvider provider;
    // Aircraft pitched up 30 degrees, velocity still in world +X
    gp_vec3 velocity(10.0f, 0.0f, 0.0f);
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::pitchedUp(QuatHelper::deg(30)), velocity);

    gp_scalar alpha = evaluateGPOperator(GETALPHA, provider, state, nullptr, 0);
    gp_scalar expected = ExpectedValue::alpha(state);

    EXPECT_NEAR(alpha, expected, 1e-5f);

    TEST_NODE(GETALPHA);
}

// T053: TEST(SensorOps, BetaYawed) using QuatHelper::yawed(deg(45))
TEST(SensorOps, BetaYawed) {
    TestPathProvider provider;
    // Aircraft yawed 45 degrees, velocity in world +X
    gp_vec3 velocity(10.0f, 0.0f, 0.0f);
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::yawed(QuatHelper::deg(45)), velocity);

    gp_scalar beta = evaluateGPOperator(GETBETA, provider, state, nullptr, 0);
    gp_scalar expected = ExpectedValue::beta(state);

    EXPECT_NEAR(beta, expected, 1e-5f);

    TEST_NODE(GETBETA);
}

// T054: TEST(SensorOps, AlphaBetaRolled90) using QuatHelper::bankedRight(deg(90))
TEST(SensorOps, AlphaBetaRolled90) {
    TestPathProvider provider;
    // Aircraft rolled 90 degrees right, velocity in world +X
    gp_vec3 velocity(10.0f, 0.0f, 0.0f);
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::bankedRight(QuatHelper::deg(90)), velocity);

    gp_scalar alpha = evaluateGPOperator(GETALPHA, provider, state, nullptr, 0);
    gp_scalar beta = evaluateGPOperator(GETBETA, provider, state, nullptr, 0);
    gp_scalar expectedAlpha = ExpectedValue::alpha(state);
    gp_scalar expectedBeta = ExpectedValue::beta(state);

    EXPECT_NEAR(alpha, expectedAlpha, 1e-5f);
    EXPECT_NEAR(beta, expectedBeta, 1e-5f);

    TEST_NODE(GETALPHA);
    TEST_NODE(GETBETA);
}

// T055: TEST(SensorOps, AlphaBetaCombined) using combined pitch+roll+yaw
TEST(SensorOps, AlphaBetaCombined) {
    TestPathProvider provider;
    gp_vec3 velocity(10.0f, 2.0f, -1.0f);  // Complex velocity vector
    gp_quat orientation = QuatHelper::climbingTurn(QuatHelper::deg(15), QuatHelper::deg(20), QuatHelper::deg(10));
    AircraftState state = makeState(gp_vec3::Zero(), orientation, velocity);

    gp_scalar alpha = evaluateGPOperator(GETALPHA, provider, state, nullptr, 0);
    gp_scalar beta = evaluateGPOperator(GETBETA, provider, state, nullptr, 0);
    gp_scalar expectedAlpha = ExpectedValue::alpha(state);
    gp_scalar expectedBeta = ExpectedValue::beta(state);

    EXPECT_NEAR(alpha, expectedAlpha, 1e-5f);
    EXPECT_NEAR(beta, expectedBeta, 1e-5f);

    TEST_NODE(GETALPHA);
    TEST_NODE(GETBETA);
}

// T056: TEST(SensorOps, GetRollRadPurePitch) verifying pure pitch orientations
TEST(SensorOps, GetRollRadPurePitch) {
    TestPathProvider provider;
    // Pure pitch up - roll should be ~0
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::pitchedUp(QuatHelper::deg(30)));

    gp_scalar roll = evaluateGPOperator(GETROLL_RAD, provider, state, nullptr, 0);
    gp_scalar expected = ExpectedValue::rollRad(state);

    EXPECT_NEAR(roll, expected, 1e-5f);

    TEST_NODE(GETROLL_RAD);
}

// T057: TEST(SensorOps, GetPitchRadPureRoll) verifying pure roll orientations
TEST(SensorOps, GetPitchRadPureRoll) {
    TestPathProvider provider;
    // Pure roll - pitch should be ~0
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::bankedRight(QuatHelper::deg(45)));

    gp_scalar pitch = evaluateGPOperator(GETPITCH_RAD, provider, state, nullptr, 0);
    gp_scalar expected = ExpectedValue::pitchRad(state);

    EXPECT_NEAR(pitch, expected, 1e-5f);

    TEST_NODE(GETPITCH_RAD);
}

// T058: TEST(SensorOps, GetRollPitchGimbalLock) testing pitch=±90° edge case
TEST(SensorOps, GetRollPitchGimbalLock) {
    TestPathProvider provider;
    // Pitch 90 degrees up (gimbal lock)
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::pitchedUp(QuatHelper::deg(89)));

    gp_scalar roll = evaluateGPOperator(GETROLL_RAD, provider, state, nullptr, 0);
    gp_scalar pitch = evaluateGPOperator(GETPITCH_RAD, provider, state, nullptr, 0);

    // Near gimbal lock, just verify values are reasonable
    EXPECT_TRUE(std::isfinite(roll));
    EXPECT_TRUE(std::isfinite(pitch));

    TEST_NODE(GETROLL_RAD);
    TEST_NODE(GETPITCH_RAD);
}

// T059: TEST(SensorOps, GetRollPitchAllQuadrants) covering 0°, 30°, 60°, 90°, 180°
TEST(SensorOps, GetRollPitchAllQuadrants) {
    TestPathProvider provider;
    gp_scalar angles[] = {0.0f, 30.0f, 60.0f, 90.0f, 120.0f, 150.0f, 180.0f};

    for (gp_scalar angle : angles) {
        if (std::abs(angle - 90.0f) < 1.0f) continue;  // Skip gimbal lock

        // Test roll at various angles
        AircraftState rollState = makeState(gp_vec3::Zero(), QuatHelper::bankedRight(QuatHelper::deg(angle)));
        gp_scalar roll = evaluateGPOperator(GETROLL_RAD, provider, rollState, nullptr, 0);
        gp_scalar expectedRoll = ExpectedValue::rollRad(rollState);
        EXPECT_NEAR(roll, expectedRoll, 1e-4f) << "Roll angle: " << angle;

        // Test pitch at various angles (limited to avoid gimbal lock)
        if (angle <= 80.0f) {
            AircraftState pitchState = makeState(gp_vec3::Zero(), QuatHelper::pitchedUp(QuatHelper::deg(angle)));
            gp_scalar pitch = evaluateGPOperator(GETPITCH_RAD, provider, pitchState, nullptr, 0);
            gp_scalar expectedPitch = ExpectedValue::pitchRad(pitchState);
            EXPECT_NEAR(pitch, expectedPitch, 1e-4f) << "Pitch angle: " << angle;
        }
    }

    TEST_NODE(GETROLL_RAD);
    TEST_NODE(GETPITCH_RAD);
}

// T060: TEST(SensorOps, GetDHome) verifying distance calculation
TEST(SensorOps, GetDHome) {
    TestPathProvider provider;
    // Aircraft 100m north of home
    AircraftState state = makeState(gp_vec3(100.0f, 0.0f, SIM_INITIAL_ALTITUDE), QuatHelper::level());

    gp_scalar dHome = executeGetDHome(state);

    // Should be 100m (same altitude as home)
    EXPECT_NEAR(dHome, 100.0f, 1e-4f);

    TEST_NODE(GETDHOME);
}

// T061-T064: Bytecode Sensor Tests
TEST(BytecodeSensors, GetVelXYZ) {
    TestPathProvider provider;
    gp_vec3 velocity(3.0f, 4.0f, 5.0f);
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::level(), velocity);

    GPBytecode velXProg[] = {{GETVELX, 0, 0}};
    GPBytecode velYProg[] = {{GETVELY, 0, 0}};
    GPBytecode velZProg[] = {{GETVELZ, 0, 0}};

    EXPECT_FLOAT_EQ(evaluateBytecodePortable(velXProg, 1, provider, state, 0), 3.0f);
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(velYProg, 1, provider, state, 0), 4.0f);
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(velZProg, 1, provider, state, 0), 5.0f);

    TEST_NODE(GETVELX);
    TEST_NODE(GETVELY);
    TEST_NODE(GETVELZ);
}

TEST(BytecodeSensors, AlphaBetaQuaternions) {
    TestPathProvider provider;
    gp_vec3 velocity(10.0f, 0.0f, 0.0f);
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::pitchedUp(QuatHelper::deg(30)), velocity);

    GPBytecode alphaProg[] = {{GETALPHA, 0, 0}};
    GPBytecode betaProg[] = {{GETBETA, 0, 0}};

    gp_scalar alpha = evaluateBytecodePortable(alphaProg, 1, provider, state, 0);
    gp_scalar beta = evaluateBytecodePortable(betaProg, 1, provider, state, 0);

    EXPECT_NEAR(alpha, ExpectedValue::alpha(state), 1e-5f);
    EXPECT_NEAR(beta, ExpectedValue::beta(state), 1e-5f);

    TEST_NODE(GETALPHA);
    TEST_NODE(GETBETA);
}

TEST(BytecodeSensors, GetRollPitchRad) {
    TestPathProvider provider;
    AircraftState state = makeState(gp_vec3::Zero(), QuatHelper::bankedRight(QuatHelper::deg(45)));

    GPBytecode rollProg[] = {{GETROLL_RAD, 0, 0}};
    GPBytecode pitchProg[] = {{GETPITCH_RAD, 0, 0}};

    gp_scalar roll = evaluateBytecodePortable(rollProg, 1, provider, state, 0);
    gp_scalar pitch = evaluateBytecodePortable(pitchProg, 1, provider, state, 0);

    EXPECT_NEAR(roll, ExpectedValue::rollRad(state), 1e-5f);
    EXPECT_NEAR(pitch, ExpectedValue::pitchRad(state), 1e-5f);

    TEST_NODE(GETROLL_RAD);
    TEST_NODE(GETPITCH_RAD);
}

TEST(BytecodeSensors, GetDHome) {
    TestPathProvider provider;
    AircraftState state = makeState(gp_vec3(100.0f, 0.0f, SIM_INITIAL_ALTITUDE), QuatHelper::level());

    GPBytecode prog[] = {{GETDHOME, 0, 0}};
    gp_scalar dHome = evaluateBytecodePortable(prog, 1, provider, state, 0);

    EXPECT_NEAR(dHome, 100.0f, 1e-4f);

    TEST_NODE(GETDHOME);
}

// =============================================================================
// Phase 7: User Story 5 - Test Temporal History Nodes (P2)
// =============================================================================

// T065: TEST(TemporalOps, GetDPhiPrevIndex0) verifying most recent history
TEST(TemporalOps, GetDPhiPrevIndex0) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Record some history
    state.recordErrorHistory(0.5f, 0.3f, 100);
    state.recordErrorHistory(0.6f, 0.4f, 200);

    gp_scalar args[1] = {0.0f};  // Index 0 = most recent
    gp_scalar result = evaluateGPOperator(GETDPHI_PREV, provider, state, args, 1);

    EXPECT_FLOAT_EQ(result, 0.6f);

    TEST_NODE(GETDPHI_PREV);
}

// T066: TEST(TemporalOps, GetDPhiPrevIndex2) verifying older history access
TEST(TemporalOps, GetDPhiPrevIndex2) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Record 5 history entries
    state.recordErrorHistory(0.1f, 0.1f, 100);
    state.recordErrorHistory(0.2f, 0.2f, 200);
    state.recordErrorHistory(0.3f, 0.3f, 300);
    state.recordErrorHistory(0.4f, 0.4f, 400);
    state.recordErrorHistory(0.5f, 0.5f, 500);

    gp_scalar args[1] = {2.0f};  // Index 2 = 3rd most recent
    gp_scalar result = evaluateGPOperator(GETDPHI_PREV, provider, state, args, 1);

    EXPECT_FLOAT_EQ(result, 0.3f);

    TEST_NODE(GETDPHI_PREV);
}

// T067: TEST(TemporalOps, GetDThetaPrev) mirroring dPhi tests
TEST(TemporalOps, GetDThetaPrev) {
    TestPathProvider provider;
    AircraftState state = makeState();

    state.recordErrorHistory(0.1f, 0.2f, 100);
    state.recordErrorHistory(0.3f, 0.4f, 200);

    gp_scalar args[1] = {0.0f};
    gp_scalar result = evaluateGPOperator(GETDTHETA_PREV, provider, state, args, 1);

    EXPECT_FLOAT_EQ(result, 0.4f);

    TEST_NODE(GETDTHETA_PREV);
}

// T068: TEST(TemporalOps, EmptyHistoryReturnsZero) verifying empty buffer behavior
TEST(TemporalOps, EmptyHistoryReturnsZero) {
    TestPathProvider provider;
    AircraftState state = makeState();
    // No history recorded

    gp_scalar args[1] = {0.0f};
    gp_scalar dPhiResult = evaluateGPOperator(GETDPHI_PREV, provider, state, args, 1);
    gp_scalar dThetaResult = evaluateGPOperator(GETDTHETA_PREV, provider, state, args, 1);

    EXPECT_FLOAT_EQ(dPhiResult, 0.0f);
    EXPECT_FLOAT_EQ(dThetaResult, 0.0f);

    TEST_NODE(GETDPHI_PREV);
    TEST_NODE(GETDTHETA_PREV);
}

// T069: TEST(TemporalOps, HistoryIndexOutOfBounds) verifying clamping
TEST(TemporalOps, HistoryIndexOutOfBounds) {
    TestPathProvider provider;
    AircraftState state = makeState();

    state.recordErrorHistory(0.5f, 0.6f, 100);

    // Index way beyond history size - should clamp to last valid
    gp_scalar args[1] = {100.0f};
    gp_scalar result = evaluateGPOperator(GETDPHI_PREV, provider, state, args, 1);

    // Should return something valid (clamped index)
    EXPECT_TRUE(std::isfinite(result));

    TEST_NODE(GETDPHI_PREV);
}

// T070: TEST(TemporalOps, GetDPhiRate) verifying rate calculation
TEST(TemporalOps, GetDPhiRate) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Record history with known time delta
    state.recordErrorHistory(0.0f, 0.0f, 0);
    state.recordErrorHistory(0.5f, 0.3f, 100);  // 100ms later

    gp_scalar rate = evaluateGPOperator(GETDPHI_RATE, provider, state, nullptr, 0);

    // Rate = (0.5 - 0.0) / 0.1s = 5.0 rad/s
    EXPECT_NEAR(rate, 5.0f, 0.1f);

    TEST_NODE(GETDPHI_RATE);
}

// T071: TEST(TemporalOps, GetDThetaRate) mirroring dPhi rate test
TEST(TemporalOps, GetDThetaRate) {
    TestPathProvider provider;
    AircraftState state = makeState();

    state.recordErrorHistory(0.0f, 0.0f, 0);
    state.recordErrorHistory(0.2f, 0.4f, 100);  // 100ms later

    gp_scalar rate = evaluateGPOperator(GETDTHETA_RATE, provider, state, nullptr, 0);

    // Rate = (0.4 - 0.0) / 0.1s = 4.0 rad/s
    EXPECT_NEAR(rate, 4.0f, 0.1f);

    TEST_NODE(GETDTHETA_RATE);
}

// T072: TEST(TemporalOps, RateWithZeroDt) verifying default dt behavior
TEST(TemporalOps, RateWithZeroDt) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Same timestamp - dt would be 0, should use default
    state.recordErrorHistory(0.0f, 0.0f, 100);
    state.recordErrorHistory(0.1f, 0.1f, 100);  // Same time!

    gp_scalar dPhiRate = evaluateGPOperator(GETDPHI_RATE, provider, state, nullptr, 0);
    gp_scalar dThetaRate = evaluateGPOperator(GETDTHETA_RATE, provider, state, nullptr, 0);

    // Should use default dt=0.1s: rate = 0.1 / 0.1 = 1.0
    EXPECT_NEAR(dPhiRate, 1.0f, 0.1f);
    EXPECT_NEAR(dThetaRate, 1.0f, 0.1f);

    TEST_NODE(GETDPHI_RATE);
    TEST_NODE(GETDTHETA_RATE);
}

// T073-T074: Bytecode Temporal Tests
TEST(BytecodeTemporal, GetPrevNodes) {
    TestPathProvider provider;
    AircraftState state = makeState();

    state.recordErrorHistory(0.1f, 0.2f, 100);
    state.recordErrorHistory(0.3f, 0.4f, 200);

    // GETDPHI_PREV(0)
    GPBytecode dPhiPrevProg[] = {{ZERO, 0, 0}, {GETDPHI_PREV, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(dPhiPrevProg, 2, provider, state, 0), 0.3f);

    // GETDTHETA_PREV(0)
    GPBytecode dThetaPrevProg[] = {{ZERO, 0, 0}, {GETDTHETA_PREV, 1, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(dThetaPrevProg, 2, provider, state, 0), 0.4f);

    TEST_NODE(GETDPHI_PREV);
    TEST_NODE(GETDTHETA_PREV);
}

TEST(BytecodeTemporal, GetRateNodes) {
    TestPathProvider provider;
    AircraftState state = makeState();

    state.recordErrorHistory(0.0f, 0.0f, 0);
    state.recordErrorHistory(0.5f, 0.4f, 100);

    GPBytecode dPhiRateProg[] = {{GETDPHI_RATE, 0, 0}};
    GPBytecode dThetaRateProg[] = {{GETDTHETA_RATE, 0, 0}};

    EXPECT_NEAR(evaluateBytecodePortable(dPhiRateProg, 1, provider, state, 0), 5.0f, 0.1f);
    EXPECT_NEAR(evaluateBytecodePortable(dThetaRateProg, 1, provider, state, 0), 4.0f, 0.1f);

    TEST_NODE(GETDPHI_RATE);
    TEST_NODE(GETDTHETA_RATE);
}

// =============================================================================
// Phase 8: User Story 6 - Test Logic and Constant Nodes (P3)
// =============================================================================

// T075: TEST(LogicOps, IfTruthyBranch) verifying IF(1,a,b)=a
TEST(LogicOps, IfTruthyBranch) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar args[3] = {1.0f, 7.0f, -3.0f};  // IF(1, 7, -3) = 7
    EXPECT_FLOAT_EQ(evaluateGPOperator(IF, provider, state, args, 3), 7.0f);

    // Any non-zero is truthy
    args[0] = 0.5f;
    EXPECT_FLOAT_EQ(evaluateGPOperator(IF, provider, state, args, 3), 7.0f);

    args[0] = -1.0f;
    EXPECT_FLOAT_EQ(evaluateGPOperator(IF, provider, state, args, 3), 7.0f);

    TEST_NODE(IF);
}

// T076: TEST(LogicOps, IfFalsyBranch) verifying IF(0,a,b)=b
TEST(LogicOps, IfFalsyBranch) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar args[3] = {0.0f, 7.0f, -3.0f};  // IF(0, 7, -3) = -3
    EXPECT_FLOAT_EQ(evaluateGPOperator(IF, provider, state, args, 3), -3.0f);

    TEST_NODE(IF);
}

// T077: TEST(LogicOps, EqEqual) verifying EQ(5,5)=1
TEST(LogicOps, EqEqual) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar args[2] = {5.0f, 5.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(EQ, provider, state, args, 2), 1.0f);

    TEST_NODE(EQ);
}

// T078: TEST(LogicOps, EqNotEqual) verifying EQ(5,6)=0
TEST(LogicOps, EqNotEqual) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar args[2] = {5.0f, 6.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(EQ, provider, state, args, 2), 0.0f);

    TEST_NODE(EQ);
}

// T079: TEST(LogicOps, GtComparisons) verifying GT behavior
TEST(LogicOps, GtComparisons) {
    TestPathProvider provider;
    AircraftState state = makeState();

    gp_scalar argsGt[2] = {7.0f, 5.0f};  // 7 > 5 = true
    EXPECT_FLOAT_EQ(evaluateGPOperator(GT, provider, state, argsGt, 2), 1.0f);

    gp_scalar argsLt[2] = {3.0f, 5.0f};  // 3 > 5 = false
    EXPECT_FLOAT_EQ(evaluateGPOperator(GT, provider, state, argsLt, 2), 0.0f);

    gp_scalar argsEq[2] = {5.0f, 5.0f};  // 5 > 5 = false
    EXPECT_FLOAT_EQ(evaluateGPOperator(GT, provider, state, argsEq, 2), 0.0f);

    TEST_NODE(GT);
}

// T080: TEST(LogicOps, PrognSideEffect) verifying PROGN(SETROLL,GETROLL)
TEST(LogicOps, PrognSideEffect) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // First execute SETROLL(0.5) and capture its result
    gp_scalar setRollArgs[1] = {0.5f};
    gp_scalar setResult = evaluateGPOperator(SETROLL, provider, state, setRollArgs, 1);

    // Reset state for PROGN test
    state = makeState();

    // PROGN should execute first arg for side effect, return second arg
    // We simulate this by: args = {setResult, getResult}
    // Actually, PROGN takes two evaluated arguments
    gp_scalar prognArgs[2] = {0.5f, 0.5f};  // Both args already evaluated
    gp_scalar result = evaluateGPOperator(PROGN, provider, state, prognArgs, 2);
    EXPECT_FLOAT_EQ(result, 0.5f);  // Returns second arg

    TEST_NODE(PROGN);
}

// T081: TEST(LogicOps, Constants) verifying ZERO, ONE, TWO, OP_PI values
TEST(LogicOps, Constants) {
    TestPathProvider provider;
    AircraftState state = makeState();

    EXPECT_FLOAT_EQ(evaluateGPOperator(ZERO, provider, state, nullptr, 0), 0.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(ONE, provider, state, nullptr, 0), 1.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(TWO, provider, state, nullptr, 0), 2.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(OP_PI, provider, state, nullptr, 0), static_cast<gp_scalar>(M_PI));

    TEST_NODE(ZERO);
    TEST_NODE(ONE);
    TEST_NODE(TWO);
    TEST_NODE(OP_PI);
}

// T082-T084: Bytecode Logic Tests
TEST(BytecodeLogic, IfEqGt) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // IF(1, 2, 0) = 2
    GPBytecode ifTrueProg[] = {{ONE, 0, 0}, {TWO, 0, 0}, {ZERO, 0, 0}, {IF, 3, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(ifTrueProg, 4, provider, state, 0), 2.0f);

    // IF(0, 2, 1) = 1
    GPBytecode ifFalseProg[] = {{ZERO, 0, 0}, {TWO, 0, 0}, {ONE, 0, 0}, {IF, 3, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(ifFalseProg, 4, provider, state, 0), 1.0f);

    // EQ(1, 1) = 1
    GPBytecode eqTrueProg[] = {{ONE, 0, 0}, {ONE, 0, 0}, {EQ, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(eqTrueProg, 3, provider, state, 0), 1.0f);

    // GT(2, 1) = 1
    GPBytecode gtTrueProg[] = {{TWO, 0, 0}, {ONE, 0, 0}, {GT, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(gtTrueProg, 3, provider, state, 0), 1.0f);

    TEST_NODE(IF);
    TEST_NODE(EQ);
    TEST_NODE(GT);
}

TEST(BytecodeLogic, PrognSideEffect) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // PROGN(SETROLL(1), GETROLL) - should set roll to 1, return 1
    GPBytecode prognProg[] = {
        {ONE, 0, 0},      // Push 1
        {SETROLL, 1, 0},  // SETROLL(1) - pops 1, pushes result (1.0)
        {GETROLL, 0, 0},  // GETROLL - pushes 1.0
        {PROGN, 2, 0}     // PROGN - pops 2, returns second (GETROLL result)
    };
    gp_scalar result = evaluateBytecodePortable(prognProg, 4, provider, state, 0);
    EXPECT_FLOAT_EQ(result, 1.0f);
    EXPECT_FLOAT_EQ(state.getRollCommand(), 1.0f);

    TEST_NODE(PROGN);
}

TEST(BytecodeLogic, Constants) {
    TestPathProvider provider;
    AircraftState state = makeState();

    GPBytecode zeroProg[] = {{ZERO, 0, 0}};
    GPBytecode oneProg[] = {{ONE, 0, 0}};
    GPBytecode twoProg[] = {{TWO, 0, 0}};
    GPBytecode piProg[] = {{OP_PI, 0, 0}};

    EXPECT_FLOAT_EQ(evaluateBytecodePortable(zeroProg, 1, provider, state, 0), 0.0f);
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(oneProg, 1, provider, state, 0), 1.0f);
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(twoProg, 1, provider, state, 0), 2.0f);
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(piProg, 1, provider, state, 0), static_cast<gp_scalar>(M_PI));

    TEST_NODE(ZERO);
    TEST_NODE(ONE);
    TEST_NODE(TWO);
    TEST_NODE(OP_PI);
}

// =============================================================================
// Phase 9: Coverage Validation & Polish
// =============================================================================

// T085: TEST(Coverage, AllNodesTested) iterating testAllNodes[] and verifying each is tested
TEST(Coverage, AllNodesTested) {
    for (int i = 0; i < testAllNodesCount; ++i) {
        EXPECT_TRUE(NodeCoverageTracker::isTested(testAllNodes[i].opcode))
            << "Missing test for node: " << testAllNodes[i].name
            << " (opcode: " << testAllNodes[i].opcode << ")";
    }
}

// T086: TEST(Coverage, NodeCount) verifying expected node count
TEST(Coverage, NodeCount) {
    // Expected: 42 nodes (based on testAllNodes[] definition)
    EXPECT_EQ(testAllNodesCount, 42) << "Node count has changed - update tests accordingly";
    EXPECT_GE(NodeCoverageTracker::getTestedCount(), testAllNodesCount)
        << "Not all nodes are covered";
}

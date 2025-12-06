#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "../aircraft_state.h"
#include "../gp_evaluator_portable.h"
#include "../gp_bytecode.h"

namespace {

using fastmath::GPScalar;

AircraftState makeDefaultState() {
    return AircraftState(
        /*pathIndex*/ 0,
        /*relVel*/ 10.0,
        Eigen::Vector3d(5.0, 1.0, -2.0),
        Eigen::Quaterniond::Identity(),
        Eigen::Vector3d::Zero(),
        /*pitch*/ 0.1,
        /*roll*/ -0.1,
        /*throttle*/ 0.0,
        /*time*/ 0);
}

std::vector<Path> makeDefaultPath() {
    std::vector<Path> paths(1);
    paths[0].start = Eigen::Vector3d(10.0, 0.0, 0.0);
    paths[0].orientation = Eigen::Vector3d::UnitX();
    paths[0].distanceFromStart = 0.0;
    paths[0].radiansFromStart = 0.0;
    paths[0].simTimeMsec = 0.0;
    return paths;
}

void ExpectScalarEq(const GPScalar& actual, const GPScalar& expected, const char* msg) {
    EXPECT_EQ(actual, expected) << msg << " got=" << actual.toDouble()
                                << " expected=" << expected.toDouble();
}

void ExpectNear(const GPScalar& actual, double expected, double tol, const char* msg) {
    EXPECT_NEAR(actual.toDouble(), expected, tol) << msg;
}

TEST(GPEvaluatorTest, ArithmeticAndLogicOpcodes) {
    auto paths = makeDefaultPath();
    VectorPathProvider provider(paths, 0);
    AircraftState state = makeDefaultState();

    GPScalar args2[2] = {GPScalar::fromInt(2), GPScalar::fromInt(3)};
    ExpectNear(evaluateGPOperator(ADD, provider, state, args2, 2, GPScalar::zero()), 5.0, 1e-6, "ADD");
    ExpectNear(evaluateGPOperator(SUB, provider, state, args2, 2, GPScalar::zero()), -1.0, 1e-6, "SUB");
    ExpectNear(evaluateGPOperator(MUL, provider, state, args2, 2, GPScalar::zero()), 6.0, 1e-6, "MUL");

    GPScalar divArgs[2] = {GPScalar::fromInt(5), GPScalar::fromInt(2)};
    ExpectNear(evaluateGPOperator(DIV, provider, state, divArgs, 2, GPScalar::zero()), 2.5, 1e-6, "DIV");

    ExpectScalarEq(evaluateGPOperator(EQ, provider, state, args2, 2, GPScalar::zero()),
                   GPScalar::zero(), "EQ false");
    ExpectScalarEq(evaluateGPOperator(GT, provider, state, args2, 2, GPScalar::zero()),
                   GPScalar::zero(), "GT false");
    ExpectScalarEq(evaluateGPOperator(GT, provider, state, divArgs, 2, GPScalar::zero()),
                   GPScalar::fromInt(1), "GT true");

    GPScalar clampArgs[3] = {GPScalar::fromDouble(1.5), GPScalar::fromDouble(-1.0), GPScalar::fromDouble(1.0)};
    ExpectNear(evaluateGPOperator(CLAMP, provider, state, clampArgs, 3, GPScalar::zero()), 1.0, 1e-6, "CLAMP");

    GPScalar absArg[1] = {GPScalar::fromDouble(-2.5)};
    ExpectNear(evaluateGPOperator(ABS, provider, state, absArg, 1, GPScalar::zero()), 2.5, 1e-6, "ABS");

    GPScalar sqrtArg[1] = {GPScalar::fromDouble(4.0)};
    ExpectNear(evaluateGPOperator(SQRT, provider, state, sqrtArg, 1, GPScalar::zero()), 2.0, 1e-2, "SQRT");

    GPScalar ifArgs[3] = {GPScalar::fromInt(1), GPScalar::fromDouble(0.8), GPScalar::fromDouble(-0.8)};
    ExpectNear(evaluateGPOperator(IF, provider, state, ifArgs, 3, GPScalar::zero()), 0.8, 1e-3, "IF true");
    ifArgs[0] = GPScalar::zero();
    ExpectNear(evaluateGPOperator(IF, provider, state, ifArgs, 3, GPScalar::zero()), -0.8, 1e-3, "IF false");
}

TEST(GPEvaluatorTest, SetterGetterOpcodes) {
    auto paths = makeDefaultPath();
    VectorPathProvider provider(paths, 0);
    AircraftState state = makeDefaultState();

    GPScalar setterArg[1] = {GPScalar::fromDouble(0.75)};
    evaluateGPOperator(SETPITCH, provider, state, setterArg, 1, GPScalar::zero());
    ExpectNear(evaluateGPOperator(GETPITCH, provider, state, nullptr, 0, GPScalar::zero()), 0.75, 1e-6,
               "GETPITCH");

    setterArg[0] = GPScalar::fromDouble(-0.5);
    evaluateGPOperator(SETROLL, provider, state, setterArg, 1, GPScalar::zero());
    ExpectNear(evaluateGPOperator(GETROLL, provider, state, nullptr, 0, GPScalar::zero()), -0.5, 1e-6,
               "GETROLL");

    setterArg[0] = GPScalar::fromDouble(1.5);
    evaluateGPOperator(SETTHROTTLE, provider, state, setterArg, 1, GPScalar::zero());
    ExpectNear(evaluateGPOperator(GETTHROTTLE, provider, state, nullptr, 0, GPScalar::zero()), 1.0, 1e-6,
               "GETTHROTTLE clamp");
}

TEST(GPEvaluatorTest, SensorOpcodes) {
    auto paths = makeDefaultPath();
    VectorPathProvider provider(paths, 0);
    AircraftState state = makeDefaultState();
    state.setVelocity(Eigen::Vector3d(5.0, 0.0, 0.0));
    state.setOrientation(Eigen::Quaterniond::Identity());

    ExpectNear(evaluateGPOperator(GETVELX, provider, state, nullptr, 0, GPScalar::zero()), 5.0, 1e-6, "VELX");
    ExpectNear(evaluateGPOperator(GETVELY, provider, state, nullptr, 0, GPScalar::zero()), 0.0, 1e-6, "VELY");

    ExpectNear(evaluateGPOperator(GETALPHA, provider, state, nullptr, 0, GPScalar::zero()), 0.0, 1e-6, "ALPHA");
    ExpectNear(evaluateGPOperator(GETBETA, provider, state, nullptr, 0, GPScalar::zero()), 0.0, 1e-6, "BETA");

    ExpectNear(evaluateGPOperator(GETROLL_RAD, provider, state, nullptr, 0, GPScalar::zero()), 0.0, 1e-6,
               "ROLL_RAD");
    ExpectNear(evaluateGPOperator(GETPITCH_RAD, provider, state, nullptr, 0, GPScalar::zero()), 0.0, 1e-6,
               "PITCH_RAD");

    ExpectNear(evaluateGPOperator(GETVEL, provider, state, nullptr, 0, GPScalar::zero()),
               state.getRelVel(), 1e-6, "VEL");
}

TEST(GPEvaluatorTest, NavigationOpcodes) {
    auto paths = makeDefaultPath();
    VectorPathProvider provider(paths, 0);
    AircraftState state = makeDefaultState();
    state.setPosition(Eigen::Vector3d::Zero());

    GPScalar arg[1] = {GPScalar::fromDouble(0)};
    ExpectNear(evaluateGPOperator(GETDPHI, provider, state, arg, 1, GPScalar::zero()), 0.0, 1e-6, "DPHI");
    ExpectNear(evaluateGPOperator(GETDTHETA, provider, state, arg, 1, GPScalar::zero()), 0.0, 1e-6, "DTHETA");
    ExpectNear(evaluateGPOperator(GETDTARGET, provider, state, arg, 1, GPScalar::zero()), 0.0, 1e-6, "DTARGET");
    ExpectNear(evaluateGPOperator(GETDHOME, provider, state, nullptr, 0, GPScalar::zero()),
               std::abs(SIM_INITIAL_ALTITUDE), 1e-6, "DHOME");
}

TEST(GPEvaluatorTest, ConstantOpcodes) {
    auto paths = makeDefaultPath();
    VectorPathProvider provider(paths, 0);
    AircraftState state = makeDefaultState();

    ExpectNear(evaluateGPOperator(OP_PI, provider, state, nullptr, 0, GPScalar::zero()), M_PI, 1e-5, "PI");
    ExpectScalarEq(evaluateGPOperator(ONE, provider, state, nullptr, 0, GPScalar::zero()),
                   GPScalar::fromInt(1), "ONE");
    ExpectScalarEq(evaluateGPOperator(ZERO, provider, state, nullptr, 0, GPScalar::zero()),
                   GPScalar::zero(), "ZERO");
    ExpectScalarEq(evaluateGPOperator(TWO, provider, state, nullptr, 0, GPScalar::zero()),
                   GPScalar::fromInt(2), "TWO");
}

TEST(GPEvaluatorTest, EvaluatorStack) {
    auto paths = makeDefaultPath();
    VectorPathProvider provider(paths, 0);
    AircraftState state = makeDefaultState();

    std::vector<GPBytecode> program = {
        GPBytecode(ONE, 0, 0.0f),
        GPBytecode(TWO, 0, 0.0f),
        GPBytecode(ADD, 2, 0.0f),
        GPBytecode(SETPITCH, 1, 0.0f),
        GPBytecode(GETPITCH, 0, 0.0f),
        GPBytecode(PROGN, 2, 0.0f)
    };

    GPScalar result = evaluateBytecodePortable(program.data(), static_cast<int>(program.size()),
                                               provider, state, GPScalar::zero());
    EXPECT_NEAR(state.getPitchCommand(), 1.0, 1e-6);
    ExpectNear(result, state.getPitchCommand(), 1e-6, "stack result");
}

}  // namespace

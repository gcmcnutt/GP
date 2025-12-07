#include <gtest/gtest.h>
#define GP_TEST
#define GP_BUILD
#include "../gp_evaluator_portable.h"

namespace {

Path makePath(gp_scalar x, gp_scalar y, gp_scalar z) {
  Path p;
  p.start = gp_vec3(x, y, z);
  p.orientation = gp_vec3::UnitX();
  p.distanceFromStart = 0.0f;
  p.radiansFromStart = 0.0f;
  p.simTimeMsec = 0.0f;
  return p;
}

AircraftState makeState() {
  return AircraftState(0,
                       static_cast<gp_scalar>(10.0f),
                       gp_vec3::Zero(),
                       gp_quat::Identity(),
                       gp_vec3::Zero(),
                       0.0f, 0.0f, 0.0f,
                       0);
}

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
};

}  // namespace

TEST(GPOperatorTest, BasicMath) {
  TestPathProvider provider;
  AircraftState state = makeState();
  gp_scalar args2[2] = {2.0f, 3.0f};
  EXPECT_FLOAT_EQ(evaluateGPOperator(ADD, provider, state, args2, 2), 5.0f);
  EXPECT_FLOAT_EQ(evaluateGPOperator(SUB, provider, state, args2, 2), -1.0f);
  EXPECT_FLOAT_EQ(evaluateGPOperator(MUL, provider, state, args2, 2), 6.0f);
  EXPECT_FLOAT_EQ(evaluateGPOperator(DIV, provider, state, args2, 2), static_cast<gp_scalar>(2.0 / 3.0));

  gp_scalar args1[1] = {-4.0f};
  EXPECT_FLOAT_EQ(evaluateGPOperator(ABS, provider, state, args1, 1), 4.0f);
  EXPECT_FLOAT_EQ(evaluateGPOperator(SQRT, provider, state, args1, 1), 0.0f);
}

TEST(GPOperatorTest, TrigLUT) {
  TestPathProvider provider;
  AircraftState state = makeState();
  gp_scalar args1[1] = {static_cast<gp_scalar>(M_PI / 6)};
  gp_scalar sinVal = evaluateGPOperator(SIN, provider, state, args1, 1);
  gp_scalar cosVal = evaluateGPOperator(COS, provider, state, args1, 1);
  EXPECT_NEAR(sinVal, 0.5f, 0.02f);
  EXPECT_NEAR(cosVal, static_cast<gp_scalar>(std::sqrt(3.0) / 2.0), 0.02f);

  gp_scalar args2[2] = {0.5f, 0.5f};
  gp_scalar atanVal = evaluateGPOperator(ATAN2, provider, state, args2, 2);
  EXPECT_NEAR(atanVal, static_cast<gp_scalar>(M_PI / 4), 0.02f);
}

TEST(GPOperatorTest, ClampAndLogic) {
  TestPathProvider provider;
  AircraftState state = makeState();
  gp_scalar args3[3] = {5.0f, 0.0f, 1.0f};
  EXPECT_FLOAT_EQ(evaluateGPOperator(CLAMP, provider, state, args3, 3), 1.0f);

  gp_scalar argsIf[3] = {1.0f, 7.0f, -3.0f};
  EXPECT_FLOAT_EQ(evaluateGPOperator(IF, provider, state, argsIf, 3), 7.0f);
  argsIf[0] = 0.0f;
  EXPECT_FLOAT_EQ(evaluateGPOperator(IF, provider, state, argsIf, 3), -3.0f);
}

TEST(GPOperatorTest, ControlAndSensors) {
  TestPathProvider provider;
  AircraftState state = makeState();

  // Control side effects clamp to [-1,1]
  gp_scalar setArgs[1] = {2.5f};
  EXPECT_FLOAT_EQ(evaluateGPOperator(SETPITCH, provider, state, setArgs, 1), 1.0f);
  EXPECT_FLOAT_EQ(state.getPitchCommand(), 1.0f);
  setArgs[0] = -2.0f;
  EXPECT_FLOAT_EQ(evaluateGPOperator(SETROLL, provider, state, setArgs, 1), -1.0f);
  EXPECT_FLOAT_EQ(state.getRollCommand(), -1.0f);

  // Sensor reads reflect state
  state.setRelVel(9.0f);
  EXPECT_FLOAT_EQ(evaluateGPOperator(GETVEL, provider, state, nullptr, 0), 9.0f);
}

TEST(GPOperatorTest, NavigationAngles) {
  TestPathProvider provider;
  provider.paths[0].start = gp_vec3(1.0f, 0.0f, 0.0f);
  AircraftState state = makeState();
  state.setPosition(gp_vec3::Zero());
  state.setOrientation(gp_quat::Identity());

  gp_scalar args[1] = {0.0f};
  EXPECT_NEAR(executeGetDPhi(provider, state, args[0]), 0.0f, 1e-4f);
  EXPECT_NEAR(executeGetDTheta(provider, state, args[0]), 0.0f, 1e-4f);
  EXPECT_NEAR(executeGetDTarget(provider, state, args[0]), -0.9f, 1e-4f);
  EXPECT_NEAR(executeGetDHome(state), (gp_vec3(0.0f, 0.0f, SIM_INITIAL_ALTITUDE) - state.getPosition()).norm(), 1e-4f);
}

TEST(BytecodeTest, SimpleProgram) {
  TestPathProvider provider;
  AircraftState state = makeState();
  GPBytecode program[] = {
      {ZERO, 0, 0.0f},
      {ONE, 0, 0.0f},
      {ADD, 2, 0.0f},
      {TWO, 0, 0.0f},
      {MUL, 2, 0.0f}
  };
  gp_scalar result = evaluateBytecodePortable(program, 5, provider, state, 0.0f);
  EXPECT_FLOAT_EQ(result, 2.0f);
}

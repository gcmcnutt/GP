# Quickstart: GP Evaluator Regression Tests

**Date**: 2026-03-06 | **Branch**: `001-gp-eval-tests`

## Prerequisites

- GP repository cloned at `~/GP`
- CMake 3.10+ installed
- g++ with C++17 support
- Eigen3 development headers (`sudo apt-get install libeigen3-dev`)

## Build & Run Tests

### Initial Build

```bash
# From scratch - creates build directory and compiles with tests enabled
cd ~/GP/autoc && bash rebuild.sh

# Or manually:
cd ~/GP
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_AUTOC_TESTS=ON ../autoc
make
```

### Run Tests

```bash
# Option 1: Via ctest (recommended)
cd ~/GP/build && ctest --output-on-failure

# Option 2: Direct execution
cd ~/GP/build && ./autoc_tests

# Option 3: Via make (runs as part of build)
cd ~/GP && make
```

### Incremental Rebuild After Test Changes

```bash
# Fast rebuild - just recompiles changed files
cd ~/GP && make

# Tests run automatically as part of build (run_autoc_tests target)
```

## Test File Location

All tests are in a single file:
```
autoc/tests/gp_evaluator_tests.cc
```

## Writing New Tests

### Basic Pattern

```cpp
#include <gtest/gtest.h>
#include "../gp_evaluator_portable.h"

TEST(TestSuiteName, TestName) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Setup state/path as needed
    // ...

    // Evaluate operator
    gp_scalar args[N] = {...};
    gp_scalar result = evaluateGPOperator(OP_NAME, provider, state, args, N);

    // Verify
    EXPECT_FLOAT_EQ(result, expected);

    // Register coverage
    TEST_NODE(OP_NAME);
}
```

### Quaternion Test Pattern

```cpp
TEST(Navigation, GetDPhiWithBankedAircraft) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Use QuatHelper for orientation setup
    state.setOrientation(QuatHelper::bankedRight(QuatHelper::deg(45)));
    state.setPosition(gp_vec3::Zero());
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 0.0f);

    // Compute expected using ExpectedValue helper
    gp_scalar expected = ExpectedValue::dPhi(provider, state, 0.0f);

    // Evaluate and compare
    gp_scalar args[1] = {0.0f};
    gp_scalar result = executeGetDPhi(provider, state, args[0]);
    EXPECT_NEAR(result, expected, 1e-5f);

    TEST_NODE(GETDPHI);
}
```

### Bytecode Test Pattern

```cpp
TEST(BytecodeEval, MathOps) {
    TestPathProvider provider;
    AircraftState state = makeState();

    // Build bytecode program: ADD(2, 3)
    GPBytecode prog[] = {
        {CONSTANT, 0, 2.0f},  // Push 2
        {CONSTANT, 0, 3.0f},  // Push 3
        {ADD, 2, 0.0f}        // Pop 2 args, push result
    };

    gp_scalar result = evaluateBytecodePortable(prog, 3, provider, state, 0.0f);
    EXPECT_FLOAT_EQ(result, 5.0f);

    TEST_NODE(ADD);  // Count towards coverage
}
```

## Test Organization

Tests are organized by user story priority:

| Test Suite | Nodes Covered | Priority |
|------------|---------------|----------|
| `MathOps` | ADD, SUB, MUL, DIV, SIN, COS, ATAN2, CLAMP, ABS, SQRT, MIN, MAX | P1 |
| `ControlOps` | SETPITCH, SETROLL, SETTHROTTLE, GETPITCH, GETROLL, GETTHROTTLE | P1 |
| `NavigationOps` | GETDPHI, GETDTHETA, GETDTARGET, GETDHOME | P1 |
| `SensorOps` | GETVEL, GETVELX, GETVELY, GETVELZ, GETALPHA, GETBETA, GETROLL_RAD, GETPITCH_RAD | P1 |
| `TemporalOps` | GETDPHI_PREV, GETDTHETA_PREV, GETDPHI_RATE, GETDTHETA_RATE | P2 |
| `LogicOps` | IF, EQ, GT, PROGN, ZERO, ONE, TWO, OP_PI | P3 |
| `BytecodeEval` | All nodes via bytecode interpreter | P1 |
| `Coverage` | Validates 100% node coverage | P1 |

## Coverage Verification

The final test in the suite verifies 100% coverage:

```cpp
TEST(Coverage, AllNodesTested) {
    extern const NodeDef allNodes[];
    extern const int allNodesCount;

    for (int i = 0; i < allNodesCount; ++i) {
        EXPECT_TRUE(NodeCoverageTracker::isTested(allNodes[i].op))
            << "Missing test for node: " << allNodes[i].name;
    }

    // Expected: 42 nodes
    EXPECT_EQ(NodeCoverageTracker::getTestedCount(), allNodesCount);
}
```

## Debugging Failed Tests

```bash
# Run with verbose output
cd ~/GP/build && ./autoc_tests --gtest_filter="TestSuite.TestName" --gtest_break_on_failure

# Run specific test suite
cd ~/GP/build && ./autoc_tests --gtest_filter="NavigationOps.*"

# List all tests
cd ~/GP/build && ./autoc_tests --gtest_list_tests
```

## Common Issues

### "fastSin/fastCos not declared"

Ensure test file includes the evaluator header with GP_TEST defined:
```cpp
// In CMakeLists.txt, already configured:
target_compile_definitions(autoc_tests PRIVATE GP_BUILD GP_TEST)
```

### "allNodes not found"

The `allNodes[]` array is in `autoc-eval.cc`. For test access, add extern declarations:
```cpp
extern const NodeDef allNodes[];
extern const int allNodesCount;
```

### Quaternion tests produce unexpected values

Check coordinate convention: ZYX intrinsic (yaw→pitch→roll). Verify QuatHelper::fromEuler matches codebase convention in COORDINATE_CONVENTIONS.md.

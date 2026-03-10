# Research: GP Evaluator Regression Tests

**Date**: 2026-03-06 | **Branch**: `001-gp-eval-tests`

## Research Questions

### R1: Node Enumeration & Coverage Detection (FR-010, SC-005)

**Question**: How to detect when a new node is added to `allNodes[]` but lacks test coverage?

**Decision**: Runtime enumeration check at test startup

**Rationale**:
- `allNodes[]` is defined in `autoc-eval.cc` with `allNodesCount` computed at compile time
- A test-time check can iterate `allNodes[]` and verify each opcode has a corresponding test
- Compile-time static_assert is not feasible because test coverage is runtime-determined

**Implementation Pattern**:
```cpp
// In test file: register each tested node
static std::set<GPOperator> testedNodes;

#define TEST_NODE(op) testedNodes.insert(op)

// At test suite end: verify coverage
TEST(Coverage, AllNodesTested) {
    extern const NodeDef allNodes[];
    extern const int allNodesCount;

    for (int i = 0; i < allNodesCount; ++i) {
        EXPECT_TRUE(testedNodes.count(allNodes[i].op) > 0)
            << "Missing test for node: " << allNodes[i].name;
    }
}
```

**Alternatives Considered**:
1. Manual checklist in code comments - Rejected: error-prone, no enforcement
2. Compile-time template metaprogramming - Rejected: overly complex for the benefit
3. External script to parse test file - Rejected: fragile, requires tooling maintenance

---

### R2: Quaternion Helper Design (FR-007)

**Question**: Best approach for constructing test quaternions from Euler angles and axis-angle?

**Decision**: Eigen `AngleAxis` composition following ZYX aerospace convention

**Rationale**:
- Existing codebase uses `Eigen::AngleAxis<gp_scalar>` for quaternion construction (see `minisim.cc:191-193`)
- Euler angles follow ZYX (yaw-pitch-roll) intrinsic sequence per `COORDINATE_CONVENTIONS.md`
- Eigen's quaternion multiplication handles composition correctly

**Implementation Pattern**:
```cpp
namespace QuatHelper {
    // From Euler angles (roll, pitch, yaw) in radians - ZYX intrinsic sequence
    inline gp_quat fromEuler(gp_scalar roll, gp_scalar pitch, gp_scalar yaw) {
        return gp_quat(Eigen::AngleAxis<gp_scalar>(yaw, gp_vec3::UnitZ())) *
               gp_quat(Eigen::AngleAxis<gp_scalar>(pitch, gp_vec3::UnitY())) *
               gp_quat(Eigen::AngleAxis<gp_scalar>(roll, gp_vec3::UnitX()));
    }

    // From axis-angle
    inline gp_quat fromAxisAngle(const gp_vec3& axis, gp_scalar angle) {
        return gp_quat(Eigen::AngleAxis<gp_scalar>(angle, axis.normalized()));
    }

    // Common flight attitudes (convenience)
    inline gp_quat level() { return gp_quat::Identity(); }
    inline gp_quat pitchedUp(gp_scalar angle) { return fromEuler(0, angle, 0); }
    inline gp_quat pitchedDown(gp_scalar angle) { return fromEuler(0, -angle, 0); }
    inline gp_quat bankedRight(gp_scalar angle) { return fromEuler(angle, 0, 0); }
    inline gp_quat bankedLeft(gp_scalar angle) { return fromEuler(-angle, 0, 0); }
    inline gp_quat yawed(gp_scalar angle) { return fromEuler(0, 0, angle); }

    // Degrees to radians helper
    constexpr gp_scalar deg(gp_scalar d) { return d * GP_PI / 180.0f; }
}
```

**Alternatives Considered**:
1. Direct quaternion component specification - Rejected: error-prone, unintuitive
2. Rotation matrix conversion - Rejected: unnecessary intermediate step
3. External quaternion library - Rejected: Eigen already available and trusted

---

### R3: LUT-Based Expected Value Computation

**Question**: How to achieve bit-accurate SIN/COS comparisons?

**Decision**: Pre-compute expected values using the same `fastSin`/`fastCos` functions

**Rationale**:
- Per clarification session: "Pre-compute expected values using the same LUT, compare exactly"
- The LUT functions (`fastSin`, `fastCos`) in `gp_evaluator_portable.cc` are the canonical implementation
- Test computes expected value, evaluator computes result, both use same LUT → bit-accurate match

**Implementation Pattern**:
```cpp
// Expose LUT functions for testing (via GP_TEST define)
#ifdef GP_TEST
gp_scalar fastSin(gp_scalar angle);
gp_scalar fastCos(gp_scalar angle);
gp_scalar fastAtan2(gp_scalar y, gp_scalar x);
#endif

TEST(MathOps, SinBitAccurate) {
    gp_scalar angles[] = {0, GP_PI/6, GP_PI/4, GP_PI/2, GP_PI, 3*GP_PI/2, 2*GP_PI};
    for (gp_scalar a : angles) {
        gp_scalar expected = fastSin(a);  // Pre-compute using same LUT
        gp_scalar args[1] = {a};
        gp_scalar result = evaluateGPOperator(SIN, provider, state, args, 1);
        EXPECT_FLOAT_EQ(result, expected);  // Bit-accurate comparison
    }
}
```

**Alternatives Considered**:
1. Compare against `std::sin` with tolerance - Rejected: LUT has interpolation differences
2. Store golden values in test - Rejected: not self-documenting, platform-dependent
3. Skip LUT tests - Rejected: LUT is critical path, must be tested

---

### R4: Bytecode Test Pattern (FR-001 dual-mode coverage)

**Question**: How to structure tests for bytecode interpreter parity?

**Decision**: Parallel test structure - each tree evaluator test has corresponding bytecode test

**Rationale**:
- Constitution requires dual-mode parity (Principle III)
- Existing pattern in `gp_evaluator_tests.cc` shows `BytecodeTest::SimpleProgram`
- Bytecode uses same operators, just different execution path (stack-based vs recursive)

**Implementation Pattern**:
```cpp
// Tree evaluator test
TEST(TreeEval, AddSubMul) {
    gp_scalar args2[2] = {2.0f, 3.0f};
    EXPECT_FLOAT_EQ(evaluateGPOperator(ADD, provider, state, args2, 2), 5.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(SUB, provider, state, args2, 2), -1.0f);
    EXPECT_FLOAT_EQ(evaluateGPOperator(MUL, provider, state, args2, 2), 6.0f);
}

// Corresponding bytecode test
TEST(BytecodeEval, AddSubMul) {
    // ADD(2, 3) = 5
    GPBytecode addProg[] = {{CONSTANT, 0, 2.0f}, {CONSTANT, 0, 3.0f}, {ADD, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(addProg, 3, provider, state, 0), 5.0f);

    // SUB(2, 3) = -1
    GPBytecode subProg[] = {{CONSTANT, 0, 2.0f}, {CONSTANT, 0, 3.0f}, {SUB, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(subProg, 3, provider, state, 0), -1.0f);

    // MUL(2, 3) = 6
    GPBytecode mulProg[] = {{CONSTANT, 0, 2.0f}, {CONSTANT, 0, 3.0f}, {MUL, 2, 0}};
    EXPECT_FLOAT_EQ(evaluateBytecodePortable(mulProg, 3, provider, state, 0), 6.0f);
}
```

**Alternatives Considered**:
1. Single parameterized test for both modes - Rejected: different call signatures
2. Skip bytecode tests - Rejected: violates Constitution Principle III
3. Only test bytecode for subset - Rejected: incomplete coverage

---

### R5: Ground Truth for Quaternion-Dependent Nodes

**Question**: How to compute expected values for GETDPHI, GETDTHETA, GETALPHA, GETBETA?

**Decision**: Use Eigen library operations as independent ground truth

**Rationale**:
- Per clarification session: "Eigen library quaternion operations - independent computation path"
- Eigen is trusted, well-tested math library
- Expected values computed using Eigen's rotation and vector operations
- Tests verify evaluator matches Eigen's computation

**Implementation Pattern**:
```cpp
// Expected value computation using Eigen
gp_scalar computeExpectedDPhi(const gp_vec3& target, const gp_vec3& position, const gp_quat& orientation) {
    // World-frame direction to target
    gp_vec3 toTarget = (target - position).normalized();

    // Transform to body frame using inverse rotation
    gp_vec3 bodyTarget = orientation.conjugate() * toTarget;

    // Roll angle needed: atan2(y, x) in body frame XY plane
    return std::atan2(bodyTarget.y(), bodyTarget.x());
}

TEST(Navigation, GetDPhiWithBankedAircraft) {
    // Setup: aircraft at origin, banked 45° right, target ahead in world frame
    state.setPosition(gp_vec3::Zero());
    state.setOrientation(QuatHelper::bankedRight(QuatHelper::deg(45)));
    provider.paths[0].start = gp_vec3(10.0f, 0.0f, 0.0f);  // Target ahead

    // Compute expected using Eigen operations
    gp_scalar expected = computeExpectedDPhi(
        provider.paths[0].start, state.getPosition(), state.getOrientation());

    // Evaluate and compare
    gp_scalar args[1] = {0.0f};
    gp_scalar result = executeGetDPhi(provider, state, args[0]);
    EXPECT_NEAR(result, expected, 1e-5f);  // Near for FP accumulation
}
```

**Note**: `EXPECT_NEAR` used for quaternion-dependent nodes due to potential floating-point accumulation differences. Tolerance is tight (1e-5f) but not bit-exact.

**Alternatives Considered**:
1. Hand-computed rotation matrices - Rejected: error-prone, less trustworthy than Eigen
2. Reference implementation in Python - Rejected: adds external dependency
3. Compare against current implementation - Rejected: violates "clean slate" principle

---

## Summary of Decisions

| Research Question | Decision | Key Artifact |
|-------------------|----------|--------------|
| R1: Node enumeration | Runtime coverage check | `testedNodes` set + Coverage test |
| R2: Quaternion helpers | Eigen AngleAxis + ZYX convention | `QuatHelper` namespace |
| R3: LUT expected values | Pre-compute using same LUT | Expose `fastSin`/`fastCos` for test |
| R4: Bytecode tests | Parallel test structure | `BytecodeEval` test suite |
| R5: Quaternion ground truth | Eigen library operations | `computeExpected*` helpers |

## Dependencies Identified

- **Eigen3**: Already in use, no new dependency
- **GoogleTest 1.14.0**: Already configured via FetchContent
- **GP_TEST compile define**: Already defined for test builds (CMakeLists.txt:104)

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| LUT functions not exposed for testing | Use `GP_TEST` ifdef to expose, already patterned in codebase |
| Quaternion edge cases (gimbal lock) | Explicitly test pitch=±90° scenarios per spec |
| Test suite too slow | Batch similar tests, avoid redundant state setup |

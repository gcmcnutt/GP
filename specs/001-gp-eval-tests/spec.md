# Feature Specification: GP Evaluator Regression Tests

**Feature Branch**: `001-gp-eval-tests`
**Created**: 2026-03-06
**Status**: Draft
**Input**: Full regression testing of GP evaluator to ensure all nodes are tested with controlled inputs and verified outputs

## Clarifications

### Session 2026-03-06

- Q: What tolerance should be used for SIN/COS LUT comparisons? → A: None - bit-accurate. Pre-compute expected values using the same LUT, compare exactly.
- Q: Should bytecode interpreter tests be in scope for 100% node coverage? → A: Yes, in scope. This is about math accuracy - tests are clean slate validation, not regression freeze. Do NOT assume current eval is correct; expect impl may be wrong and tests will prove it.
- Q: What is the ground truth source for quaternion-dependent node expected values? → A: Eigen library quaternion operations - independent computation path using the same trusted math library.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Test Mathematical Operations (Priority: P1)

As a developer, I want comprehensive tests for all mathematical operator nodes (ADD, SUB, MUL, DIV, SIN, COS, ATAN2, CLAMP, ABS, SQRT, MIN, MAX) so that I can verify basic arithmetic correctness and catch regressions when modifying the evaluator.

**Why this priority**: Mathematical operations are foundational - all other nodes depend on math working correctly. These are also the simplest to test and provide immediate validation that the test framework itself works.

**Independent Test**: Run the math node tests in isolation - they require no path data or complex aircraft state setup.

**Acceptance Scenarios**:

1. **Given** two scalar inputs (a, b), **When** ADD is evaluated, **Then** result equals a + b exactly (FP32 bit-accurate)
2. **Given** inputs that would cause division by zero, **When** DIV is evaluated with divisor=0, **Then** result is 0 (not NaN or crash)
3. **Given** various angle inputs across all quadrants, **When** SIN/COS are evaluated, **Then** results match pre-computed LUT values exactly (bit-accurate)
4. **Given** inputs (y, x) in all four quadrants, **When** ATAN2 is evaluated, **Then** result matches expected angle with correct sign and quadrant

---

### User Story 2 - Test Control Side Effects (Priority: P1)

As a developer, I want tests for control output nodes (SETPITCH, SETROLL, SETTHROTTLE) and their corresponding getters (GETPITCH, GETROLL, GETTHROTTLE) so that I can verify the SET/GET cycle works correctly and values are properly clamped to [-1, 1].

**Why this priority**: Control outputs are the actual outputs of the GP controller. If these fail, the aircraft cannot be controlled regardless of how sophisticated the evolved program is.

**Independent Test**: Create an AircraftState, execute SET commands, verify GET returns clamped values.

**Acceptance Scenarios**:

1. **Given** an input value > 1.0, **When** SETPITCH is evaluated, **Then** pitchCommand is clamped to 1.0 and returned
2. **Given** an input value < -1.0, **When** SETROLL is evaluated, **Then** rollCommand is clamped to -1.0 and returned
3. **Given** SETTHROTTLE(0.5) followed by GETTHROTTLE, **When** both are evaluated in sequence, **Then** GETTHROTTLE returns 0.5

---

### User Story 3 - Test Path Following Nodes with Time Offsets and 3D Quaternions (Priority: P1)

As a developer, I want comprehensive tests for path-following nodes (GETDPHI, GETDTHETA, GETDTARGET) with various time offset arguments and diverse 3D quaternion orientations so that I can verify the time-based path indexing works correctly and angle calculations handle all flight attitudes.

**Why this priority**: Path following is the suspected regression area. The recent change from step offsets to time offsets needs thorough validation. Additionally, GETDPHI and GETDTHETA involve quaternion transforms from world to body frame that must be tested with non-identity orientations.

**Independent Test**: Create a controlled path with known positions and timestamps, set aircraft position/orientation using explicit quaternions, verify computed angles match expected values.

**Quaternion Test Coverage for Path-Following**:
- Identity orientation (baseline)
- Aircraft pitched up/down while tracking level target
- Aircraft banked while tracking target in various positions
- Aircraft in climbing/descending turn (combined rotations)
- Target positions in all octants (±X, ±Y, ±Z combinations)

**Acceptance Scenarios**:

1. **Given** aircraft at origin facing +X with target directly ahead on path, **When** GETDPHI(0) and GETDTHETA(0) are evaluated, **Then** both return 0 (no correction needed)
2. **Given** target above the aircraft nose, **When** GETDTHETA(0) is evaluated, **Then** result is positive (pitch up to reach target)
3. **Given** target below the aircraft nose, **When** GETDTHETA(0) is evaluated, **Then** result is negative (pitch down)
4. **Given** target to the right of aircraft, **When** GETDPHI(0) is evaluated, **Then** result is positive (roll right)
5. **Given** a path with multiple waypoints at 100ms intervals and current index at waypoint 5, **When** GETDPHI(2) is evaluated (200ms ahead), **Then** the function targets waypoint 7's position
6. **Given** a path with 10 waypoints and current index at waypoint 8, **When** GETDPHI(5) is evaluated (500ms ahead), **Then** the function clamps to the last waypoint (index 9)
7. **Given** current index at waypoint 2, **When** GETDPHI(-3) is evaluated (300ms behind), **Then** the function clamps to waypoint 0 (cannot go negative)
8. **Given** aircraft banked 45° right with target directly ahead in world frame, **When** GETDPHI(0) is evaluated, **Then** result accounts for current bank angle in body-frame calculation
9. **Given** aircraft pitched up 30° with target at same altitude, **When** GETDTHETA(0) is evaluated, **Then** result indicates pitch-down correction needed

---

### User Story 4 - Test Sensor Nodes with 3D Quaternion Orientations (Priority: P1)

As a developer, I want tests for aircraft state sensor nodes (GETVEL, GETVELX, GETVELY, GETVELZ, GETALPHA, GETBETA, GETROLL_RAD, GETPITCH_RAD, GETDHOME) with comprehensive 3D quaternion orientations so that I can verify sensors accurately reflect aircraft state in all flight attitudes.

**Why this priority**: ELEVATED to P1 because these sensors involve quaternion transforms that must work correctly in 3D space. Testing only identity quaternions misses critical rotation-dependent behavior in GETALPHA, GETBETA, and body-frame velocity transforms.

**Independent Test**: Set known velocity/orientation on AircraftState using explicit quaternion construction, evaluate each sensor, verify returned values.

**Quaternion Test Coverage**:
- Identity (level flight)
- Pure pitch rotations (nose up 30°, 60°, 90°; nose down 30°, 60°, 90°)
- Pure roll rotations (bank left/right 30°, 60°, 90°, 180°)
- Pure yaw rotations (heading 0°, 90°, 180°, 270°)
- Combined pitch+roll (climbing turn configurations)
- Combined pitch+yaw (heading change while climbing/descending)
- Combined all three (general attitude)

**Acceptance Scenarios**:

1. **Given** velocity vector (20, 0, 0) in NED frame, **When** GETVEL is evaluated, **Then** result equals 20 (magnitude)
2. **Given** velocity (10, 5, -3) in NED frame, **When** GETVELX, GETVELY, GETVELZ are evaluated, **Then** results are 10, 5, -3 respectively
3. **Given** aircraft flying straight and level (identity orientation) with forward velocity, **When** GETALPHA and GETBETA are evaluated, **Then** both return 0
4. **Given** aircraft pitched up 30° with horizontal velocity, **When** GETALPHA is evaluated, **Then** result equals -30° (angle of attack in body frame)
5. **Given** aircraft yawed 45° right with northward velocity, **When** GETBETA is evaluated, **Then** result reflects sideslip angle
6. **Given** aircraft rolled 90° with forward velocity and no vertical velocity, **When** GETALPHA and GETBETA are evaluated, **Then** results correctly transform to body frame
7. **Given** aircraft with known Euler angles, **When** GETROLL_RAD and GETPITCH_RAD are evaluated, **Then** results match expected angles including gimbal lock edge cases

---

### User Story 5 - Test Temporal History Nodes (Priority: P2)

As a developer, I want tests for temporal history nodes (GETDPHI_PREV, GETDTHETA_PREV, GETDPHI_RATE, GETDTHETA_RATE) so that I can verify historical error lookback and rate calculations work correctly.

**Why this priority**: These are newer nodes added for derivative-style control. They depend on correct history buffer management.

**Independent Test**: Populate AircraftState history buffer with known values, evaluate PREV and RATE nodes, verify results.

**Acceptance Scenarios**:

1. **Given** history buffer with dPhi values [0.1, 0.2, 0.3] (most recent first), **When** GETDPHI_PREV(0) is evaluated, **Then** result is 0.1
2. **Given** history buffer with dPhi values [0.1, 0.2, 0.3], **When** GETDPHI_PREV(2) is evaluated, **Then** result is 0.3
3. **Given** empty history buffer, **When** GETDPHI_PREV(0) is evaluated, **Then** result is 0.0
4. **Given** two consecutive history entries with 100ms time delta and dPhi change of 0.1, **When** GETDPHI_RATE is evaluated, **Then** result is approximately 1.0 rad/s

---

### User Story 6 - Test Logic and Constant Nodes (Priority: P3)

As a developer, I want tests for logic nodes (IF, EQ, GT, PROGN) and constant nodes (ZERO, ONE, TWO, OP_PI) so that I can verify conditional execution and constant values.

**Why this priority**: These are well-established, unlikely to regress, but should be covered for completeness.

**Independent Test**: Evaluate each node with known inputs, verify outputs.

**Acceptance Scenarios**:

1. **Given** condition=1.0, **When** IF(1.0, 7.0, -3.0) is evaluated, **Then** result is 7.0 (truthy branch)
2. **Given** condition=0.0, **When** IF(0.0, 7.0, -3.0) is evaluated, **Then** result is -3.0 (falsy branch)
3. **Given** two equal values, **When** EQ(5.0, 5.0) is evaluated, **Then** result is 1.0
4. **Given** PROGN(SETROLL(0.5), GETROLL), **When** evaluated, **Then** SETROLL side effect occurs and GETROLL value (0.5) is returned
5. **Given** OP_PI is evaluated, **When** compared to M_PI, **Then** values are equal within FP32 precision

---

### Edge Cases

- **Path boundary - before start**: When time offset would go before t=0, index clamps to 0
- **Path boundary - past end**: When time offset would exceed path length, index clamps to last waypoint
- **Division by zero**: DIV(x, 0) returns 0, not NaN
- **Negative sqrt**: SQRT(-1) returns 0, not NaN
- **NaN path argument**: If arg to GETDPHI is NaN, function returns current index (defensive)
- **All quadrants for ATAN2**: Test (+,+), (+,-), (-,+), (-,-) quadrant combinations
- **Gimbal lock region**: Test GETROLL_RAD and GETPITCH_RAD near pitch=±90°
- **Empty history buffer**: Temporal nodes return 0.0 when no history available
- **History index out of bounds**: GETDPHI_PREV with index > HISTORY_SIZE clamps to valid range
- **Rate calculation with zero dt**: GETDPHI_RATE uses default dt=0.1s if actual dt < 0.001s
- **Quaternion edge cases**:
  - 180° rotation (quaternion w≈0): Test body-frame transforms with inverted aircraft
  - Quaternion normalization: Ensure near-unit quaternions still produce valid transforms
  - Quaternion sign ambiguity: q and -q represent same rotation, transforms must handle both
  - Small angle approximations: Test that calculations remain accurate for very small rotation angles

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: Test suite MUST achieve 100% coverage of ALL operator node types defined in allNodes[] (autoc-eval.cc) - currently ~40 nodes including temporal state nodes. Coverage includes BOTH tree evaluator AND bytecode interpreter paths
- **FR-002**: Tests MUST use FP32 (gp_scalar) precision throughout to match xiao hardware constraints
- **FR-003**: Tests MUST be deterministic - same inputs produce bit-accurate same outputs
- **FR-004**: Tests MUST be clean slate validation - expected values derived from mathematical definitions, NOT from current evaluator behavior. Do not assume current implementation is correct; tests exist to prove correctness or expose bugs
- **FR-005**: Path-following tests MUST verify time-based indexing with the getPathIndex() function
- **FR-006**: Tests MUST cover 3D quaternion orientations for ALL orientation-dependent nodes - not just identity quaternions. This includes GETDPHI, GETDTHETA, GETALPHA, GETBETA, GETROLL_RAD, GETPITCH_RAD, and body-frame velocity transforms
- **FR-007**: Test framework MUST provide quaternion construction helpers for common orientations (pitch/roll/yaw combinations, Euler-to-quaternion, axis-angle)
- **FR-008**: Tests MUST integrate with existing GoogleTest infrastructure in autoc/tests/ and run via existing `make test` / `ctest` commands
- **FR-009**: Tests MUST run as part of the standard build verification (make test or ctest) - no separate test target required
- **FR-010**: When new nodes are added to allNodes[], there MUST be a mechanism to detect untested nodes (compile-time or test-time enumeration check)

### Key Entities

- **AircraftState**: Complete aircraft state including position, orientation (gp_quat quaternion), velocity, control commands, and temporal history buffer
- **gp_quat**: Quaternion type (Eigen::Quaternionf or portable equivalent) used for all 3D rotations and attitude representation
- **Path**: Waypoint with position, orientation, distance/radians from start, and simulation timestamp (simTimeMsec)
- **PathProvider**: Abstract interface for path access (TestPathProvider for mocking, VectorPathProvider for real paths)
- **gp_scalar**: FP32 type used consistently for all calculations
- **QuaternionTestHelper**: Test utility for constructing quaternions from Euler angles, axis-angle, or common flight attitudes

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: 100% of ALL operator node types in allNodes[] have at least one dedicated test case (enforced by enumeration check at test time)
- **SC-002**: All path-following nodes (GETDPHI, GETDTHETA, GETDTARGET) have tests covering positive, negative, and zero time offsets
- **SC-003**: All orientation-dependent nodes have tests with non-identity quaternions covering: pure rotations (pitch/roll/yaw), combined rotations, and edge cases (gimbal lock region, 180° rotations)
- **SC-004**: Test suite passes with zero failures on both x86_64 and ARM builds (self-referential comparison)
- **SC-005**: Adding a new node to allNodes[] without a corresponding test produces a clear warning or failure at `make test` time
- **SC-006**: Test suite completes in under 5 seconds
- **SC-007**: Test coverage validated by existing `make test` infrastructure - no separate test targets required

## Assumptions

- Existing GoogleTest infrastructure in autoc/tests/ is functional and will be extended (not replaced)
- The gp_evaluator_tests.cc file provides the foundation patterns (TestPathProvider, makeState, makePath) which will be enhanced with quaternion helpers
- FP32 precision is sufficient for all operations (matching xiao hardware)
- "Clean slate" testing means expected values are derived from mathematical definitions and geometric first principles, not from current implementation behavior
- Path waypoints are stored with simTimeMsec timestamps at 100ms intervals (SIM_TIME_STEP_MSEC)
- The existing `make test` / `ctest` infrastructure is the target execution environment - tests extend this, not replace it
- Quaternion math uses Eigen library quaternion operations (or gp_quat portable equivalent) for consistency with production code

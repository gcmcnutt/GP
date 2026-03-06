# Data Model: GP Evaluator Regression Tests

**Date**: 2026-03-06 | **Branch**: `001-gp-eval-tests`

This document defines the test framework entities and their relationships for the GP evaluator regression test suite.

## Test Framework Entities

### TestPathProvider

Mock implementation of `PathProvider` for controlled path testing.

```cpp
struct TestPathProvider : public PathProvider {
    std::vector<Path> paths;          // Configurable path waypoints
    int currentIndex{0};              // Current position along path

    // Interface implementation
    const Path& getPath(int index) const override;
    int getCurrentIndex() const override;
    int getPathSize() const override;

    // Test helpers
    void addPath(const Path& p);
    void setCurrentIndex(int idx);
    void clear();
};
```

**Relationships**:
- Implements `PathProvider` abstract interface
- Contains `std::vector<Path>` for waypoint storage
- Used by all navigation node tests (GETDPHI, GETDTHETA, GETDTARGET)

**Validation Rules**:
- `getPath(index)` clamps to valid range [0, paths.size()-1]
- Empty paths vector → single default path at origin

---

### Path (existing, extended for testing)

Waypoint structure representing a point along the flight path.

```cpp
struct Path {
    gp_vec3 start;              // Position in NED frame (meters)
    gp_vec3 orientation;        // Orientation vector (deprecated, use quaternion)
    gp_scalar distanceFromStart; // Cumulative distance from path start
    gp_scalar radiansFromStart;  // Cumulative heading change
    gp_scalar simTimeMsec;       // Simulation timestamp (milliseconds)
};
```

**Test Helper Factory**:
```cpp
Path makePath(gp_scalar x, gp_scalar y, gp_scalar z, gp_scalar timeMsec = 0.0f);
Path makePathWithTime(const gp_vec3& pos, gp_scalar timeMsec);
```

**Validation Rules**:
- `simTimeMsec` should be monotonically increasing along path
- Standard path uses 100ms intervals (`SIM_TIME_STEP_MSEC`)

---

### AircraftState (existing, test construction patterns)

Complete aircraft state for GP evaluation context.

**Key Fields for Testing**:
```cpp
class AircraftState {
    gp_vec3 position_;           // NED position (meters)
    gp_quat orientation_;        // Earth→body quaternion
    gp_vec3 velocity_;           // NED velocity (m/s)
    gp_scalar pitchCommand_;     // Control output [-1, 1]
    gp_scalar rollCommand_;      // Control output [-1, 1]
    gp_scalar throttleCommand_;  // Control output [-1, 1]
    gp_scalar relVel_;           // Relative airspeed (m/s)

    // Temporal history (for PREV/RATE nodes)
    gp_scalar dPhiHistory_[HISTORY_SIZE];
    gp_scalar dThetaHistory_[HISTORY_SIZE];
    unsigned long timeHistory_[HISTORY_SIZE];
    int historyIndex_;
    int historyCount_;
};
```

**Test Helper Factory**:
```cpp
// Default state: at origin, identity orientation, zero velocity
AircraftState makeState();

// Configurable state
AircraftState makeState(const gp_vec3& position, const gp_quat& orientation);
AircraftState makeState(const gp_vec3& position, const gp_quat& orientation,
                        const gp_vec3& velocity);
```

---

### QuatHelper (new)

Test utility namespace for quaternion construction.

```cpp
namespace QuatHelper {
    // Core construction
    gp_quat fromEuler(gp_scalar roll, gp_scalar pitch, gp_scalar yaw);
    gp_quat fromAxisAngle(const gp_vec3& axis, gp_scalar angle);

    // Flight attitude presets
    gp_quat level();                      // Identity
    gp_quat pitchedUp(gp_scalar angle);   // Nose up (positive pitch)
    gp_quat pitchedDown(gp_scalar angle); // Nose down (negative pitch)
    gp_quat bankedRight(gp_scalar angle); // Right wing down
    gp_quat bankedLeft(gp_scalar angle);  // Left wing down
    gp_quat yawed(gp_scalar angle);       // Heading change

    // Combined attitudes
    gp_quat climbingTurn(gp_scalar pitch, gp_scalar roll, gp_scalar yaw);

    // Utility
    constexpr gp_scalar deg(gp_scalar degrees);  // Degrees to radians
}
```

**Convention**: ZYX intrinsic sequence (yaw→pitch→roll) per COORDINATE_CONVENTIONS.md

---

### NodeCoverageTracker (new)

Runtime coverage tracking for 100% node coverage verification.

```cpp
class NodeCoverageTracker {
    static std::set<GPOperator> testedOps_;

public:
    static void markTested(GPOperator op);
    static bool isTested(GPOperator op);
    static std::vector<GPOperator> getUntestedOps();
    static int getTestedCount();
    static void reset();
};

// Usage macro for test registration
#define TEST_NODE(op) NodeCoverageTracker::markTested(op)
```

**Integration**:
- Each test calls `TEST_NODE(op)` for operators it tests
- Final coverage test iterates `allNodes[]` and checks against tracker
- Test failure if any node in `allNodes[]` is missing from tracker

---

### ExpectedValueComputers (new)

Helper functions for computing ground truth expected values.

```cpp
namespace ExpectedValue {
    // LUT-based (bit-accurate) - uses same LUT functions as evaluator
    gp_scalar sin(gp_scalar angle);
    gp_scalar cos(gp_scalar angle);
    gp_scalar atan2(gp_scalar y, gp_scalar x);

    // Quaternion-based (Eigen ground truth)
    gp_scalar dPhi(const PathProvider& path, const AircraftState& state, gp_scalar offset);
    gp_scalar dTheta(const PathProvider& path, const AircraftState& state, gp_scalar offset);
    gp_scalar alpha(const AircraftState& state);   // Angle of attack
    gp_scalar beta(const AircraftState& state);    // Sideslip angle
    gp_scalar rollRad(const AircraftState& state); // Roll from quaternion
    gp_scalar pitchRad(const AircraftState& state);// Pitch from quaternion
}
```

---

## Entity Relationships

```text
┌─────────────────────────────────────────────────────────────────┐
│                        Test Suite                                │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────────┐      ┌──────────────────┐                 │
│  │ TestPathProvider │◄─────│ Path (vector)    │                 │
│  └────────┬─────────┘      └──────────────────┘                 │
│           │                                                      │
│           │ provides path data to                               │
│           ▼                                                      │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │              evaluateGPOperator() / bytecode              │   │
│  │  ┌─────────────┐  ┌─────────────┐  ┌─────────────────┐   │   │
│  │  │ AircraftState│  │ QuatHelper  │  │ ExpectedValue   │   │   │
│  │  │ (under test) │  │ (construct) │  │ (ground truth)  │   │   │
│  │  └─────────────┘  └─────────────┘  └─────────────────┘   │   │
│  └──────────────────────────────────────────────────────────┘   │
│           │                                                      │
│           │ registers coverage                                  │
│           ▼                                                      │
│  ┌──────────────────┐                                           │
│  │ NodeCoverageTracker│──► Coverage test validates 100%         │
│  └──────────────────┘                                           │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## State Transitions

### AircraftState History Buffer

For temporal nodes (GETDPHI_PREV, GETDTHETA_PREV, GETDPHI_RATE, GETDTHETA_RATE):

```text
Empty → Recording → Full (ring buffer wraps)

States:
- Empty: historyCount_ = 0, all history arrays zeroed
- Recording: 0 < historyCount_ < HISTORY_SIZE
- Full: historyCount_ = HISTORY_SIZE, historyIndex_ wraps

Transitions:
- recordErrorHistory() increments historyCount_ until HISTORY_SIZE
- historyIndex_ = (historyIndex_ + 1) % HISTORY_SIZE after each record
```

**Test Coverage Required**:
1. Empty history → PREV nodes return 0.0
2. Partial history → PREV returns available values, 0.0 for unavailable
3. Full history → PREV correctly indexes ring buffer
4. Rate calculation with dt < 0.001s → uses default dt

---

## Type Definitions Reference

```cpp
// From gp_types.h
using gp_scalar = float;                           // FP32 throughout
using gp_vec3 = Eigen::Matrix<gp_scalar, 3, 1>;    // 3D vector
using gp_quat = Eigen::Quaternion<gp_scalar>;      // Quaternion

// Constants
constexpr gp_scalar GP_PI = 3.14159265358979323846f;
constexpr gp_scalar GP_TWO_PI = 2.0f * GP_PI;
constexpr gp_scalar GP_HALF_PI = GP_PI / 2.0f;
```

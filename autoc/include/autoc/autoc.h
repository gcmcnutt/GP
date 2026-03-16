#ifndef AUTOC_H
#define AUTOC_H

#include <vector>

#include "autoc/types.h"
#include "autoc/nn/topology.h"

// Forward declarations
struct WorkerContext;

// ========================================================================
// SIMPLIFIED FITNESS FUNCTION (v3 - two objectives, path-relative scaling)
// Goal: Stay close to rabbit, fly smooth. Strategy belongs in higher layers.
// See specs/FITNESS_SIMPLIFY_20260221.md for rationale.
// ========================================================================

// Nonlinear fitness uses pow(x/NORM, POWER).  Below NORM the cost compresses
// (tolerated), above NORM it amplifies (penalized disproportionately).
//
// Distance uses a V-shaped function centered on a TARGET following distance:
//   pow(fabs(distance - TARGET) / NORM, POWER)
// This directly encodes the desired 5-10m trailing position:
//   - At TARGET: zero cost (perfect following)
//   - Too close (<5m): penalty rises — discourages overshooting on turns
//   - Too far (>10m): penalty rises — pressures controller to close gap
//   - Fractional power keeps far-tail (55m+) from dominating the sum
// With TARGET=7.5, NORM=5, POWER=1.5:
//   5m → 0.35  |  7.5m → 0  |  10m → 0.35  |  15m → 1.84  |  24m → 5.49
#define DISTANCE_TARGET 7.5
#define DISTANCE_NORM 5.0
#define DISTANCE_POWER 1.5

// Attitude delta: 200 deg/s nominal at 10Hz → 0.349 rad/step
// Keep higher power to still crush tumbles/flips (was 1.5, tested 0.75)
#define ATTITUDE_NORM 0.349
#define ATTITUDE_POWER 1.5

// Attitude scaling: computed per-path from path geometry (no manual tuning)
// attitude_scale = path_distance / path_turn_rad (meters per radian)

// Intercept-budget fitness scaling (see specs/005-entry-fitness-ramp)
// Ramps distance+attitude penalty from floor to ceiling over estimated intercept time
#define INTERCEPT_SCALE_FLOOR 0.1           // Minimum scale factor at t=0
#define INTERCEPT_SCALE_CEILING 1.0         // Maximum scale factor (full penalty)
#define INTERCEPT_BUDGET_MAX 15.0           // Maximum budget cap in seconds
#define INTERCEPT_TURN_RATE (M_PI / 4.0)    // Estimated turn rate (rad/s) for budget calc (~45°/s)

// Entry position safe bounds (~15m margin from crash boundaries)
#define ENTRY_SAFE_RADIUS static_cast<gp_scalar>(55.0f)       // SIM_PATH_RADIUS_LIMIT - 15m
#define ENTRY_SAFE_ALT_MIN static_cast<gp_scalar>(-22.0f)     // SIM_MIN_ELEVATION - 15m (NED, shallowest)
#define ENTRY_SAFE_ALT_MAX static_cast<gp_scalar>(-105.0f)    // SIM_MAX_ELEVATION + 15m (NED, deepest)

// Crash penalty: soft lexicographic multiplier
// Completion dominates (1e6 scale), quality provides gradient within similar completion levels
#define CRASH_COMPLETION_WEIGHT 1e6         // Multiplier for (1 - fraction_completed)

struct WindScenarioConfig {
  unsigned int windSeed = 0;
  int windVariantIndex = 0;
};

struct ScenarioDescriptor {
  std::vector<std::vector<Path>> pathList;
  std::vector<WindScenarioConfig> windScenarios;
  unsigned int windSeed = 0;
  int pathVariantIndex = -1;   // -1 = unset/aggregated
  int windVariantIndex = -1;   // -1 = unset/aggregated
};

// Global aircraft state used by minisim evaluation loop
extern AircraftState aircraftState;

#endif

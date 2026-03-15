
// autoc.cc

/* -------------------------------------------------------------------
From skeleton/skeleton.cc
------------------------------------------------------------------- */

#include <boost/asio.hpp>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <new>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <thread>
#include <chrono>
#include <memory>
#include <cstdint>
#include <limits>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <cstring>
#include <mutex>

#include "gp.h"
#include "gp_bytecode.h"
#include "gpconfig.h"
#include "minisim.h"
#include "threadpool.h"
#include "autoc.h"
#include "logger.h"
#include "pathgen.h"
#include "config_manager.h"
#include "variation_generator.h"
#include "nn_population.h"
#include "nn_serialization.h"
#include "nn_evaluator_portable.h"

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/PutObjectRequest.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/client/ClientConfiguration.h>

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace std;
namespace logging = boost::log;

Logger logger;
std::atomic_bool printEval = false; // verbose (used for rendering best of population)
std::atomic_bool enableDeterministicTestLogging = false; // enable detailed logging during determinism test

// FNV-1a hash for verification (same as dtest_tracker)
static uint64_t hashData(const void* data, size_t len) {
    const uint8_t* bytes = static_cast<const uint8_t*>(data);
    uint64_t hash = 0xcbf29ce484222325ULL;
    for (size_t i = 0; i < len; i++) {
        hash ^= bytes[i];
        hash *= 0x100000001b3ULL;
    }
    return hash;
}

// Hash serialized AircraftState vector to avoid padding issues
static uint64_t hashAircraftStateVector(const std::vector<AircraftState>& vec) {
    std::ostringstream oss(std::ios::binary);
    boost::archive::binary_oarchive oa(oss);
    oa << vec;
    std::string serialized = oss.str();
    return hashData(serialized.data(), serialized.size());
}
std::vector<MyGP*> tasks = std::vector<MyGP*>();

// Forward declarations for global variables
extern std::ofstream fout;
extern std::ofstream bout;
extern std::atomic_ulong nanDetector;

std::vector<std::vector<Path>> generationPaths;
std::vector<ScenarioDescriptor> generationScenarios;

// VARIATIONS1: Global sigma parameters, initialized at startup from config
static VariationSigmas gVariationSigmas = {0.0, 0.0, 0.0, 0.0, 0.0};
// Individual variation enable flags (from config)
static bool gEnableEntryVariations = false;
static bool gEnableWindVariations = false;

// Variable rabbit speed: Global config, initialized at startup from config
static RabbitSpeedConfig gRabbitSpeedConfig = RabbitSpeedConfig::defaultConfig();

// RAMP_LANDSCAPE: Gradual variation scaling (see specs/RAMP_LANDSCAPE.md)
static int gCurrentGeneration = 0;      // Updated at start of each generation
static int gTotalGenerations = 1;       // Set from config at startup
static int gVariationRampStep = 0;      // Set from config at startup (0 = disabled)

/**
 * Compute variation scale for current generation.
 * Scale ramps from 0.0 to 1.0 over the course of training.
 *
 * @return Scale factor to apply to variation offsets (0.0 to 1.0)
 */
static double computeVariationScale() {
    if (gVariationRampStep <= 0) return 1.0;  // Disabled - full variations
    if (gTotalGenerations <= 0) return 1.0;   // Safety check

    // Stepped: quantize to step boundaries
    int stepIndex = gCurrentGeneration / gVariationRampStep;
    int totalSteps = gTotalGenerations / gVariationRampStep;
    if (totalSteps <= 0) return 1.0;

    return std::min(1.0, static_cast<double>(stepIndex) / static_cast<double>(totalSteps));
}

// ============================================================================
// Intercept-Budget Fitness Scaling (see specs/005-entry-fitness-ramp)
// ============================================================================

/**
 * Compute per-step intercept scale factor.
 * Quadratic ramp from FLOOR to CEILING over the intercept budget.
 *
 * @param stepTimeSec  Current simulation step time in seconds
 * @param budgetSec    Estimated intercept budget in seconds
 * @return Scale factor in [FLOOR, CEILING]
 */
static double computeInterceptScale(double stepTimeSec, double budgetSec) {
    if (budgetSec <= 0.0) return INTERCEPT_SCALE_CEILING;  // No budget → full penalty
    double t = std::min(1.0, stepTimeSec / budgetSec);
    return INTERCEPT_SCALE_FLOOR + (INTERCEPT_SCALE_CEILING - INTERCEPT_SCALE_FLOOR) * t * t;
}

/**
 * Estimate time-to-intercept from initial conditions.
 * Crude geometric estimate: turn_time + closure_time + rabbit_compensation.
 *
 * @param displacement    Horizontal distance to path start (meters)
 * @param headingOffset   Heading offset from path tangent (radians)
 * @param aircraftSpeed   Aircraft speed (m/s)
 * @param rabbitSpeed     Rabbit (target) speed (m/s)
 * @return Estimated budget in seconds, clamped to [0, INTERCEPT_BUDGET_MAX]
 */
static double computeInterceptBudget(double displacement, double headingOffset,
                                      double aircraftSpeed, double rabbitSpeed) {
    if (aircraftSpeed <= 0.0) return INTERCEPT_BUDGET_MAX;
    double turn_time = fabs(headingOffset) / INTERCEPT_TURN_RATE;
    double closure_time = displacement / aircraftSpeed;
    double rabbit_compensation = closure_time * (rabbitSpeed / aircraftSpeed);
    double budget = turn_time + closure_time + rabbit_compensation;
    return std::clamp(budget, 0.0, INTERCEPT_BUDGET_MAX);
}

// ============================================================================
// Single PRNG Architecture: Pre-fetched scenario variations
// See specs/SINGLE_PRNG.md for design
// ============================================================================

/**
 * Pre-computed variations for a single scenario (wind variant).
 * Generated from GPrand() at startup, reused every generation.
 */
struct ScenarioVariations {
    unsigned int windSeed;                      // Seed for crrcsim CRRC_Random
    VariationOffsets entryOffsets;              // Heading, roll, pitch, speed, windDir
    std::vector<RabbitSpeedPoint> rabbitProfile;// Time-tagged speeds
};

// Global pre-computed table (indexed by wind scenario index 0..N-1)
static std::vector<ScenarioVariations> gScenarioVariations;
static unsigned int gPathSeed = 0;        // Derived from GPrand() or RandomPathSeedB override
static bool gPathSeedFromOverride = false;// True if RandomPathSeedB was used

/**
 * Pre-fetch all scenario variations from GPrand() at startup.
 * Called once after GPInit(), before any GP evolution.
 *
 * When a variation type is disabled, we still consume GPrand values to maintain
 * deterministic PRNG sequence, but store defaults in the table instead.
 * This way the table contains the actual values to send to sims - no filtering needed.
 *
 * @param numScenarios       Number of wind scenarios (windScenarioCount from config)
 * @param sigmas             Variation sigma parameters
 * @param rabbitCfg          Rabbit speed configuration
 * @param randomPathSeedB    Override path seed (-1 = derive from GPrand)
 * @param enableEntry        If false, store default entry offsets (0.0, 1.0 for speed)
 * @param enableWind         If false, store default wind offset (0.0)
 */
static void prefetchAllVariations(int numScenarios, const VariationSigmas& sigmas,
                                   const RabbitSpeedConfig& rabbitCfg, int randomPathSeedB,
                                   bool enableEntry, bool enableWind) {
    gScenarioVariations.clear();
    gScenarioVariations.reserve(numScenarios);

    // Derive pathSeed (unless overridden by config)
    if (randomPathSeedB == -1) {
        gPathSeed = static_cast<unsigned int>(rng::randLong());
        gPathSeedFromOverride = false;
    } else {
        gPathSeed = static_cast<unsigned int>(randomPathSeedB);
        gPathSeedFromOverride = true;
    }

    double totalDurationSec = SIM_TOTAL_TIME_MSEC / 1000.0;

    // When wind variations disabled, use same seed for all scenarios (identical thermals/gusts)
    // Still consume GPrand to maintain deterministic sequence
    unsigned int baseWindSeed = static_cast<unsigned int>(rng::randLong());

    for (int i = 0; i < numScenarios; i++) {
        ScenarioVariations sv;

        // Wind seed for crrcsim: unique per scenario if enabled, same for all if disabled
        if (enableWind) {
            sv.windSeed = (i == 0) ? baseWindSeed : static_cast<unsigned int>(rng::randLong());
        } else {
            // Disabled: all scenarios use same seed (identical thermals/gusts)
            // Still consume GPrand to keep PRNG sequence deterministic
            if (i > 0) { (void)rng::randLong(); }
            sv.windSeed = baseWindSeed;
        }

        // Entry/wind variations: always consume GPrand to maintain deterministic sequence,
        // but store defaults if that variation type is disabled
        VariationOffsets generated = generateVariationsFromGPrand(sigmas);

        if (enableEntry) {
            sv.entryOffsets.entryHeadingOffset = generated.entryHeadingOffset;
            sv.entryOffsets.entryRollOffset = generated.entryRollOffset;
            sv.entryOffsets.entryPitchOffset = generated.entryPitchOffset;
            sv.entryOffsets.entrySpeedFactor = generated.entrySpeedFactor;
            sv.entryOffsets.entryNorthOffset = generated.entryNorthOffset;
            sv.entryOffsets.entryEastOffset = generated.entryEastOffset;
            sv.entryOffsets.entryAltOffset = generated.entryAltOffset;
        } else {
            // Defaults: no offset
            sv.entryOffsets.entryHeadingOffset = 0.0;
            sv.entryOffsets.entryRollOffset = 0.0;
            sv.entryOffsets.entryPitchOffset = 0.0;
            sv.entryOffsets.entrySpeedFactor = 1.0;
            sv.entryOffsets.entryNorthOffset = 0.0;
            sv.entryOffsets.entryEastOffset = 0.0;
            sv.entryOffsets.entryAltOffset = 0.0;
        }

        if (enableWind) {
            sv.entryOffsets.windDirectionOffset = generated.windDirectionOffset;
        } else {
            // Default: no offset
            sv.entryOffsets.windDirectionOffset = 0.0;
        }

        // Rabbit speed profile (consume ~50-200 GPrand values)
        sv.rabbitProfile = generateSpeedProfileFromGPrand(rabbitCfg, totalDurationSec);

        gScenarioVariations.push_back(std::move(sv));
    }
}

/**
 * Log pre-fetched variations at startup for verification.
 * Format matches spec in SINGLE_PRNG.md.
 */
static void logPrefetchedVariations(int numScenarios, long gpSeed) {
    *logger.info() << endl;
    *logger.info() << "=== Pre-fetched Scenario Variations (GPSeed=" << gpSeed << ") ===" << endl;
    *logger.info() << "PathSeed: " << gPathSeed
                   << " (override: " << (gPathSeedFromOverride ? "yes" : "no")
                   << ")" << endl;
    *logger.info() << "Scenarios: " << numScenarios << endl;
    *logger.info() << endl;

    *logger.info() << "Scenario  WindSeed    Heading°   Roll°   Pitch°  Speed%  WindDir°  RabbitSpd      North°  East°  Down°" << endl;
    *logger.info() << "--------  ----------  --------  ------  ------  ------  --------  ---------      ------  -----  -----" << endl;

    for (int i = 0; i < numScenarios; i++) {
        const auto& sv = gScenarioVariations[i];

        // Compute rabbit speed stats
        double minSpd = sv.rabbitProfile.empty() ? 0 : sv.rabbitProfile[0].speed;
        double maxSpd = minSpd, sumSpd = 0;
        for (const auto& pt : sv.rabbitProfile) {
            if (pt.speed < minSpd) minSpd = pt.speed;
            if (pt.speed > maxSpd) maxSpd = pt.speed;
            sumSpd += pt.speed;
        }
        double avgSpd = sv.rabbitProfile.empty() ? 0 : sumSpd / sv.rabbitProfile.size();

        // Format rabbit speed as min/avg/max with 1 decimal
        std::ostringstream rabbitStr;
        rabbitStr << std::fixed << std::setprecision(1) << minSpd << "/" << avgSpd << "/" << maxSpd;

        std::ostringstream line;
        line << std::setw(4) << i << "      "
             << "0x" << std::hex << std::setw(8) << std::setfill('0') << sv.windSeed
             << std::dec << std::setfill(' ')
             << "  " << std::setw(7) << std::fixed << std::setprecision(2) << radToDeg(sv.entryOffsets.entryHeadingOffset)
             << "  " << std::setw(6) << radToDeg(sv.entryOffsets.entryRollOffset)
             << "  " << std::setw(6) << radToDeg(sv.entryOffsets.entryPitchOffset)
             << "  " << std::setw(5) << std::setprecision(1) << ((sv.entryOffsets.entrySpeedFactor - 1.0) * 100) << "%"
             << "  " << std::setw(7) << std::setprecision(2) << radToDeg(sv.entryOffsets.windDirectionOffset)
             << "  " << std::setw(9) << rabbitStr.str();

        // Position offsets columns
        const auto& o = sv.entryOffsets;
        line << "    " << std::setprecision(1)
             << std::setw(6) << o.entryNorthOffset
             << std::setw(7) << o.entryEastOffset
             << std::setw(7) << o.entryAltOffset;
        *logger.info() << line.str() << endl;
    }
    *logger.info() << endl;
}

// Helper to populate variation offsets in ScenarioMetadata from pre-fetched table
// Table already contains correct values: defaults when disabled, variations when enabled
// (filtering is done at prefetch time, not here)
// RAMP_LANDSCAPE: Offsets are scaled by computeVariationScale() for gradual introduction
static void populateVariationOffsets(ScenarioMetadata& meta) {
  int windIdx = meta.windVariantIndex;
  if (windIdx >= 0 && windIdx < static_cast<int>(gScenarioVariations.size())) {
    const auto& v = gScenarioVariations[windIdx].entryOffsets;

    // Get current scale (0.0 to 1.0 over training)
    double scale = computeVariationScale();

    // Angular offsets: scale toward 0
    meta.entryHeadingOffset = v.entryHeadingOffset * scale;
    meta.entryRollOffset = v.entryRollOffset * scale;
    meta.entryPitchOffset = v.entryPitchOffset * scale;
    meta.windDirectionOffset = v.windDirectionOffset * scale;

    // Speed factor: interpolate from 1.0 toward pre-computed value
    // At scale=0: factor=1.0 (nominal), at scale=1: factor=pre-computed
    meta.entrySpeedFactor = 1.0 + scale * (v.entrySpeedFactor - 1.0);

    // Position offsets: scale and clamp to safe arena bounds
    // RAMP_LANDSCAPE scaling applies here too (gradual introduction)
    double northScaled = v.entryNorthOffset * scale;
    double eastScaled = v.entryEastOffset * scale;
    double altScaled = v.entryAltOffset * scale;

    // Clamp horizontal to safe radius
    double horizDist = sqrt(northScaled * northScaled + eastScaled * eastScaled);
    if (horizDist > ENTRY_SAFE_RADIUS) {
      double clampFactor = ENTRY_SAFE_RADIUS / horizDist;
      northScaled *= clampFactor;
      eastScaled *= clampFactor;
    }

    // Clamp altitude to safe NED bounds (ENTRY_SAFE_ALT_MAX < ENTRY_SAFE_ALT_MIN since NED Down)
    if (altScaled < ENTRY_SAFE_ALT_MAX) altScaled = ENTRY_SAFE_ALT_MAX;
    if (altScaled > ENTRY_SAFE_ALT_MIN) altScaled = ENTRY_SAFE_ALT_MIN;

    meta.entryNorthOffset = northScaled;
    meta.entryEastOffset = eastScaled;
    meta.entryAltOffset = altScaled;
  }
  // If windIdx out of range, leave defaults (shouldn't happen)
}

// Helper to apply variable rabbit speed to a path using pre-fetched profile
// Recomputes simTimeMsec for each path point based on the speed profile
static void applySpeedProfileToPath(std::vector<Path>& path, int windVariantIndex) {
  if (path.empty()) return;
  if (gRabbitSpeedConfig.sigma <= 0.0) return;  // Constant speed mode, no changes needed

  // Use pre-fetched speed profile (single PRNG architecture)
  if (windVariantIndex < 0 || windVariantIndex >= static_cast<int>(gScenarioVariations.size())) {
    return;  // Out of range, leave path unchanged
  }
  const auto& speedProfile = gScenarioVariations[windVariantIndex].rabbitProfile;

  // RAMP_LANDSCAPE: Get scale for rabbit speed variation
  double scale = computeVariationScale();
  double nominal = gRabbitSpeedConfig.nominal;

  // Recompute simTimeMsec for each path point based on variable speed profile
  gp_scalar accumTimeMsec = 0.0f;
  path[0].simTimeMsec = 0;

  for (size_t i = 1; i < path.size(); i++) {
    gp_scalar segmentDistance = path[i].distanceFromStart - path[i-1].distanceFromStart;
    double profileSpeed = getSpeedAtTime(speedProfile, static_cast<double>(accumTimeMsec) / 1000.0);
    // Scale speed toward nominal: at scale=0 use nominal, at scale=1 use profileSpeed
    double speed = nominal + scale * (profileSpeed - nominal);
    gp_scalar dt = (segmentDistance / static_cast<gp_scalar>(speed)) * 1000.0f;
    accumTimeMsec += dt;
    path[i].simTimeMsec = static_cast<int32_t>(accumTimeMsec);
  }
}
std::atomic_ulong nanDetector = 0;
std::ofstream fout;
std::ofstream bout;
EvalResults bestOfEvalResults;
EvalResults aggregatedEvalResults;
EvalResults* activeEvalCollector = nullptr;
boost::mutex evalCollectorMutex;  // Protects activeEvalCollector and aggregatedEvalResults
std::atomic<uint64_t> evaluationProgress{0};
std::atomic<uint64_t> globalScenarioCounter{0};
std::atomic<uint64_t> bakeoffPathCounter{0};
std::atomic<uint64_t> globalSimRunCounter{0};

// Global elite tracking across generations (dispatcher-side)
// Detect divergence when same gpHash + scenarios evaluates to different fitness
static gp_fitness gLastEliteFitness = std::numeric_limits<gp_fitness>::infinity();
static uint64_t gLastEliteGpHash = 0;
static int gLastEliteLength = 0;
static int gLastEliteDepth = 0;
static size_t gLastEliteNumPaths = 0;  // Number of paths in elite evaluation
static uint64_t gLastEliteScenarioHash = 0;  // Hash of scenario metadata
static int gEliteDivergenceCount = 0;  // Count of reevals where fitness differed (non-determinism)
static int gEliteReevalCount = 0;      // Count of times elite was re-evaluated (same gpHash across generations)
static EvalResults gLastEliteEvalResults;  // Full results including physics trace for divergence analysis
static EvalResults gPendingEliteEvalResults;  // Current gen's elite trace (moved to gLastEliteEvalResults at end of gen)

// Runtime flag to enable physics trace collection for elite reeval divergence analysis.
// When true, crrcsim will return detailed physics traces for elite reevals only.
// TODO: Consider adding to config system later.
static bool gEliteTraceEnabled = true;
// Flag to ensure we only flag ONE individual as elite reeval per generation
// (multiple individuals may have same hash due to cloning/crossover)
static bool gEliteReevalFlaggedThisGen = false;

static bool bitwiseEqual(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

static uint64_t computeGpHash(const EvalResults& res) {
  if (res.gpHash != 0) {
    return res.gpHash;
  }
  // gpHash must be populated; fallback only to maintain runtime safety in debug.
  if (!res.gp.empty()) {
    return hashByteVector(res.gp);
  }
  return 0;
}

// Hash scenario metadata to detect if path set changed
static uint64_t computeScenarioHash(const std::vector<ScenarioMetadata>& scenarios) {
  uint64_t hash = 0xcbf29ce484222325ULL; // FNV-1a offset basis
  for (const auto& s : scenarios) {
    // Hash: pathVariantIndex, windVariantIndex, windSeed
    hash ^= s.pathVariantIndex;
    hash *= 0x100000001b3ULL;
    hash ^= s.windVariantIndex;
    hash *= 0x100000001b3ULL;
    hash ^= s.windSeed;
    hash *= 0x100000001b3ULL;
  }
  return hash;
}


ThreadPool* threadPool;
std::string computedKeyName;

namespace {

unsigned int sanitizeStride(int stride) {
  int absStride = stride == 0 ? 1 : std::abs(stride);
  return static_cast<unsigned int>(absStride);
}

void rebuildGenerationScenarios(const std::vector<std::vector<Path>>& basePaths) {
  generationScenarios.clear();

  const AutocConfig& cfg = ConfigManager::getConfig();
  int windScenarioCount = std::max(cfg.windScenarioCount, 1);

  // Use pre-fetched windSeeds from gScenarioVariations (single PRNG architecture)
  // Fallback to 0 if pre-fetch hasn't run yet (shouldn't happen in normal flow)
  auto getWindSeed = [&](int windIdx) -> unsigned int {
    if (windIdx < static_cast<int>(gScenarioVariations.size())) {
      return gScenarioVariations[windIdx].windSeed;
    }
    return 0;  // Fallback (shouldn't happen)
  };

  if (basePaths.empty()) {
    ScenarioDescriptor scenario;
    scenario.pathList = basePaths;
    scenario.windVariantIndex = 0;
    scenario.windSeed = getWindSeed(0);
    scenario.windScenarios.push_back({getWindSeed(0), 0});
    generationScenarios.push_back(std::move(scenario));
    return;
  }

  if (cfg.demeticGrouping && cfg.demeSize > 0) {
    // DEMETIC MODE: Create one scenario per path variant.
    // Each scenario contains ONE path geometry evaluated across ALL wind conditions.
    // This allows demes to specialize on specific path geometries while being
    // robust to wind variations.
    for (size_t pathIdx = 0; pathIdx < basePaths.size(); ++pathIdx) {
      ScenarioDescriptor scenario;
      scenario.pathVariantIndex = static_cast<int>(pathIdx);
      for (int windIdx = 0; windIdx < windScenarioCount; ++windIdx) {
        scenario.pathList.push_back(basePaths[pathIdx]);
        WindScenarioConfig windScenario;
        windScenario.windSeed = getWindSeed(windIdx);
        windScenario.windVariantIndex = windIdx;
        scenario.windScenarios.push_back(windScenario);
      }
      if (!scenario.windScenarios.empty()) {
        scenario.windSeed = scenario.windScenarios.front().windSeed;
        scenario.windVariantIndex = scenario.windScenarios.front().windVariantIndex;
      } else {
        scenario.windSeed = getWindSeed(0);
        scenario.windVariantIndex = 0;
      }
      generationScenarios.push_back(std::move(scenario));
    }
  } else {
    // NON-DEMETIC MODE: Create ONE scenario containing ALL path variants × ALL winds.
    // Each individual evaluates on the same complete test suite for fair comparison.
    ScenarioDescriptor scenario;

    // Build wind scenarios list first
    for (int windIdx = 0; windIdx < windScenarioCount; ++windIdx) {
      WindScenarioConfig windScenario;
      windScenario.windSeed = getWindSeed(windIdx);
      windScenario.windVariantIndex = windIdx;
      scenario.windScenarios.push_back(windScenario);
    }

    // Add path×wind combinations in PATH-MAJOR order: [p0w0, p0w1, ..., p0w5, p1w0, ...]
    for (size_t pathIdx = 0; pathIdx < basePaths.size(); ++pathIdx) {
      for (int windIdx = 0; windIdx < windScenarioCount; ++windIdx) {
        scenario.pathList.push_back(basePaths[pathIdx]);
      }
    }

    if (!scenario.windScenarios.empty()) {
      scenario.windSeed = scenario.windScenarios.front().windSeed;
      scenario.windVariantIndex = scenario.windScenarios.front().windVariantIndex;
    } else {
      scenario.windSeed = getWindSeed(0);
      scenario.windVariantIndex = 0;
    }
    generationScenarios.push_back(std::move(scenario));
  }

  if (generationScenarios.empty()) {
    ScenarioDescriptor scenario;
    scenario.windVariantIndex = 0;
    scenario.windSeed = getWindSeed(0);
    scenario.pathList = basePaths;
    scenario.windScenarios.push_back({getWindSeed(0), 0});
    generationScenarios.push_back(std::move(scenario));
  }

  // Log scenario structure for verification (only once, not every generation)
  static bool logged = false;
  if (!logged) {
    logged = true;
    *logger.info() << "Scenario structure: " << generationScenarios.size() << " scenario(s)" << endl;
    for (size_t i = 0; i < generationScenarios.size(); ++i) {
      const auto& scenario = generationScenarios[i];
      if (scenario.pathVariantIndex == -1) {
        size_t numWinds = std::max(size_t(1), scenario.windScenarios.size());
        size_t numPaths = scenario.pathList.size() / numWinds;
        *logger.info() << "  Scenario " << i << ": ALL path variants"
                       << " (" << scenario.pathList.size() << " total flights = "
                       << numPaths << " paths × " << numWinds << " winds)" << endl;
      } else {
        *logger.info() << "  Scenario " << i << ": path variant " << scenario.pathVariantIndex
                       << " (" << scenario.pathList.size() << " flights across "
                       << scenario.windScenarios.size() << " winds)" << endl;
      }
    }
  }
}

struct PathFrame {
  gp_vec3 tangent = gp_vec3::UnitX();
  gp_vec3 insideNormal = gp_vec3::UnitZ();
  gp_vec3 binormal = gp_vec3::UnitY();
};

static gp_vec3 safeSegmentDirection(const std::vector<Path>& path, int idxA, int idxB) {
  if (path.empty()) {
    return gp_vec3::UnitX();
  }
  const int last = static_cast<int>(path.size()) - 1;
  idxA = std::clamp(idxA, 0, last);
  idxB = std::clamp(idxB, 0, last);
  if (idxA == idxB) {
    return gp_vec3::UnitX();
  }
  gp_vec3 delta = path[idxB].start - path[idxA].start;
  const gp_scalar norm = delta.norm();
  if (norm < static_cast<gp_scalar>(1e-5)) {
    return gp_vec3::UnitX();
  }
  return delta / norm;
}

static PathFrame computePathFrame(const std::vector<Path>& path, int index) {
  PathFrame frame;
  if (path.empty()) {
    return frame;
  }

  const int last = static_cast<int>(path.size()) - 1;
  const int clamped = std::clamp(index, 0, last);
  const int prevIdx = std::max(clamped - 1, 0);
  const int nextIdx = std::min(clamped + 1, last);

  frame.tangent = safeSegmentDirection(path, clamped, nextIdx);
  gp_vec3 prevDir = safeSegmentDirection(path, prevIdx, clamped);
  gp_vec3 nextDir = safeSegmentDirection(path, clamped, nextIdx);

  gp_vec3 curvature = nextDir - prevDir;
  // Remove any component that lies along the tangent
  curvature -= curvature.dot(frame.tangent) * frame.tangent;

  if (curvature.norm() < static_cast<gp_scalar>(1e-5)) {
    curvature = frame.tangent.cross(gp_vec3::UnitZ());
    if (curvature.norm() < static_cast<gp_scalar>(1e-5)) {
      curvature = frame.tangent.cross(gp_vec3::UnitY());
    }
  }
  if (curvature.norm() < static_cast<gp_scalar>(1e-5)) {
    curvature = gp_vec3::UnitZ();
  }
  frame.insideNormal = curvature.normalized();

  // Ensure the normal truly points "inside" the current turn (toward curvature)
  const gp_vec3 rawCurvature = nextDir - prevDir;
  if (rawCurvature.norm() > static_cast<gp_scalar>(1e-5) && rawCurvature.dot(frame.insideNormal) < 0.0f) {
    frame.insideNormal = -frame.insideNormal;
  }

  frame.binormal = frame.tangent.cross(frame.insideNormal);
  if (frame.binormal.norm() < static_cast<gp_scalar>(1e-5)) {
    frame.binormal = frame.tangent.cross(gp_vec3::UnitX());
  }
  if (frame.binormal.norm() < static_cast<gp_scalar>(1e-5)) {
    frame.binormal = gp_vec3::UnitY();
  }
  frame.binormal.normalize();

  // Re-orthogonalize insideNormal to eliminate accumulated numerical drift
  frame.insideNormal = (frame.binormal.cross(frame.tangent)).normalized();
  return frame;
}

int computeScenarioIndexForIndividual(int individualIndex) {
  const AutocConfig& cfg = ConfigManager::getConfig();
  int scenarioCount = std::max<int>(generationScenarios.size(), 1);

  if (cfg.demeticGrouping && cfg.demeSize > 0) {
    int demeIndex = individualIndex / cfg.demeSize;
    if (scenarioCount == 0) {
      return 0;
    }
    return demeIndex % scenarioCount;
  }

  // Without demetic grouping, all individuals evaluate on scenario 0
  // This ensures fair comparison within the population
  return 0;
}

const ScenarioDescriptor& scenarioForIndex(int scenarioIndex) {
  if (generationScenarios.empty()) {
    std::cerr << "FATAL: scenarioForIndex: generationScenarios is empty" << std::endl; exit(1);
  }
  int clampedIndex = ((scenarioIndex % static_cast<int>(generationScenarios.size())) + static_cast<int>(generationScenarios.size())) % static_cast<int>(generationScenarios.size());
  return generationScenarios[clampedIndex];
}

void warnIfScenarioMismatch() {
  const AutocConfig& cfg = ConfigManager::getConfig();
  if (cfg.demeticGrouping && cfg.demeSize > 0) {
    int demeCount = cfg.populationSize / cfg.demeSize;
    if (demeCount > 0 && generationScenarios.size() < static_cast<size_t>(demeCount)) {
      *logger.warn() << "Scenario count (" << generationScenarios.size()
                     << ") is smaller than deme count (" << demeCount
                     << "); scenarios will be reused across demes." << endl;
    }
  }
}

void appendEvalResults(EvalResults& dest, const EvalResults& src) {
  if (dest.gp.empty()) {
    dest.gp = src.gp;
  }
  // Copy gpHash from first evaluation (should be same for all scenarios of same GP)
  if (dest.gpHash == 0 && src.gpHash != 0) {
    dest.gpHash = src.gpHash;
  }
  dest.crashReasonList.insert(dest.crashReasonList.end(), src.crashReasonList.begin(), src.crashReasonList.end());
  dest.pathList.insert(dest.pathList.end(), src.pathList.begin(), src.pathList.end());
  dest.aircraftStateList.insert(dest.aircraftStateList.end(), src.aircraftStateList.begin(), src.aircraftStateList.end());
  dest.scenarioList.insert(dest.scenarioList.end(), src.scenarioList.begin(), src.scenarioList.end());
  dest.workerId = src.workerId;
  dest.workerPid = src.workerPid;
  dest.workerEvalCounter = src.workerEvalCounter;
}

} // namespace

std::string generate_iso8601_timestamp() {
  auto now = std::chrono::system_clock::now();
  auto ms_since_epoch = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count();
  auto itt = std::chrono::system_clock::to_time_t(now);
  std::ostringstream ss;
  ss << "autoc-" << INT64_MAX - ms_since_epoch << '-';
  ss << std::put_time(std::gmtime(&itt), "%FT%T");
  auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
  ss << '.' << std::setfill('0') << std::setw(3) << milliseconds.count() << 'Z';
  return ss.str();
}

std::shared_ptr<Aws::S3::S3Client> getS3Client() {
  return ConfigManager::getS3Client();
}

class MyPopulation : public GPPopulation
{
public:
  // Constructor (mandatory)
  MyPopulation(GPVariables& GPVar_, GPAdfNodeSet& adfNs_) :
    GPPopulation(GPVar_, adfNs_) {
  }

  // Construct from AutocConfig (bridge until GP code is deleted in T017)
  MyPopulation(AutocConfig& cfg, GPAdfNodeSet& adfNs_) :
    GPPopulation(gpVarsFromConfig(cfg), adfNs_) {
  }

  static GPVariables& gpVarsFromConfig(AutocConfig& cfg) {
    static GPVariables gv;
    gv.PopulationSize = cfg.populationSize;
    gv.NumberOfGenerations = cfg.numberOfGenerations;
    gv.CreationType = cfg.creationType;
    gv.CrossoverProbability = cfg.crossoverProbability;
    gv.CreationProbability = cfg.creationProbability;
    gv.MaximumDepthForCreation = cfg.maximumDepthForCreation;
    gv.MaximumDepthForCrossover = cfg.maximumDepthForCrossover;
    gv.SelectionType = cfg.selectionType;
    gv.TournamentSize = cfg.tournamentSize;
    gv.DemeticGrouping = cfg.demeticGrouping;
    gv.DemeSize = cfg.demeSize;
    gv.DemeticMigProbability = cfg.demeticMigProbability;
    gv.SwapMutationProbability = cfg.swapMutationProbability;
    gv.ShrinkMutationProbability = cfg.shrinkMutationProbability;
    gv.AddBestToNewPopulation = cfg.addBestToNewPopulation;
    gv.SteadyState = cfg.steadyState;
    return gv;
  }

  // Duplication (mandatory)
  MyPopulation(MyPopulation& gpo) : GPPopulation(gpo) {}
  virtual GPObject& duplicate() { return *(new MyPopulation(*this)); }

  // Creation of own class objects (mandatory)
  virtual GP* createGP(int numOfTrees) { return new MyGP(numOfTrees); }

  // Load and save (not mandatory)
  MyPopulation() {}
  virtual int isA() { return MyPopulationID; }
  virtual GPObject* createObject() { return new MyPopulation; }
  // virtual char* load (istream& is);
  // virtual void save (ostream& os);
  virtual void evaluate() override;

  virtual void endOfEvaluation() {

    evaluationProgress.store(0);

    // dispatch all the GPs now (TODO this may still work inline with evaluate)
    for (auto& task : tasks) {
      threadPool->enqueue([task](WorkerContext& context) {
        task->evalTask(context);
        });
    }

    // wait for all tasks to finish
    threadPool->wait_for_tasks();
    tasks.clear();

    // this argument should only be set if we are dumping best GP of a generation
    *logger.debug() << "printEval flag=" << printEval
                    << " activeEvalCollector=" << (activeEvalCollector ? 1 : 0) << endl;
    if (printEval) {
      bakeoffPathCounter.store(0, std::memory_order_relaxed);
      const AutocConfig& cfg = ConfigManager::getConfig();
      std::vector<int> candidateIndices;
      candidateIndices.reserve(static_cast<size_t>(std::max(1, containerSize())));

      if (cfg.demeticGrouping && cfg.demeSize > 0) {
        int demeSize = cfg.demeSize;
        if (demeSize <= 0) {
          demeSize = 1;
        }
        for (int start = 0; start < containerSize(); start += demeSize) {
          int end = std::min(start + demeSize, containerSize());
          gp_fitness bestFitness = std::numeric_limits<gp_fitness>::infinity();
          int bestIndex = start;
          for (int idx = start; idx < end; ++idx) {
            MyGP* candidate = NthMyGP(idx);
            if (!candidate) {
              continue;
            }
            gp_fitness fitness = candidate->getFitness();
            if (!std::isfinite(fitness)) {
              continue;
            }
            if (fitness < bestFitness) {
              bestFitness = fitness;
              bestIndex = idx;
            }
          }
          candidateIndices.push_back(bestIndex);
        }
      } else {
        // Non-demetic mode: just bakeoff the library's best individual
        candidateIndices.push_back(bestOfPopulation);
      }

      assert(!candidateIndices.empty() && "candidateIndices should never be empty after population scan");

      MyGP* best = nullptr;
      gp_fitness bestAggregatedFitness = std::numeric_limits<gp_fitness>::infinity();

      if (!generationScenarios.empty()) {
        EvalResults* previousCollector = activeEvalCollector;

        for (int candidateIndex : candidateIndices) {
          MyGP* candidate = NthMyGP(candidateIndex);
          if (!candidate) {
            continue;
          }

          // Enable detailed logging only when re-evaluating the previous elite.
          bool enableLogging = (candidateIndex == bestOfPopulation);
          bool prevLogging = enableDeterministicTestLogging.exchange(enableLogging, std::memory_order_relaxed);

          int originalScenarioIndex = candidate->getScenarioIndex();
          candidate->setBakeoffMode(true);
          aggregatedEvalResults.clear();
          activeEvalCollector = &aggregatedEvalResults;

          gp_fitness cumulativeFitness = 0.0;
          size_t scenariosEvaluated = 0;

          for (size_t scenarioIdx = 0; scenarioIdx < generationScenarios.size(); ++scenarioIdx) {
            candidate->setScenarioIndex(static_cast<int>(scenarioIdx));
            threadPool->enqueue([candidate](WorkerContext& context) {
              candidate->evalTask(context);
            });
            threadPool->wait_for_tasks();
            cumulativeFitness += candidate->getFitness();
            ++scenariosEvaluated;
          }

          candidate->setScenarioIndex(originalScenarioIndex);
          candidate->setBakeoffMode(false);

          gp_fitness averageFitness = scenariosEvaluated > 0
            ? cumulativeFitness / static_cast<gp_fitness>(scenariosEvaluated)
            : std::numeric_limits<gp_fitness>::infinity();

          // Restore previous logging state before moving on.
          enableDeterministicTestLogging.store(prevLogging, std::memory_order_relaxed);

          if (averageFitness < bestAggregatedFitness) {
            bestAggregatedFitness = averageFitness;
            bestOfEvalResults = aggregatedEvalResults;
            best = candidate;

          }
        }

        activeEvalCollector = previousCollector;

        // Update the best individual's fitness with the aggregated result
        // and update bestOfPopulation to point to the bakeoff winner
        if (best) {
          // Find the index of the best candidate in the population
          int bestIndex = -1;
          for (int idx = 0; idx < containerSize(); ++idx) {
            if (NthMyGP(idx) == best) {
              bestIndex = idx;
              break;
            }
          }

          // Set bakeoff winner's fitness to aggregated value
          // This ensures calculateStatistics() keeps it as best
          best->setFitness(bestAggregatedFitness);
          bestOfPopulation = bestIndex;

          // Demetic propagation (only in demetic mode)
          if (cfg.demeticGrouping && cfg.demeSize > 0 &&
              std::isfinite(bestAggregatedFitness)) {
            int demeSize = cfg.demeSize;
            int demesPropagated = 0;
            gp_scalar migProb = std::clamp(static_cast<gp_scalar>(cfg.demeticMigProbability) / static_cast<gp_scalar>(100.0f),
                                           static_cast<gp_scalar>(0.0f),
                                           static_cast<gp_scalar>(1.0f));

            for (int demeStart = 0; demeStart < containerSize(); demeStart += demeSize) {
              // Use GPrand() for deterministic migration decisions (single PRNG architecture)
              gp_scalar randVal = static_cast<gp_scalar>(rng::randDouble());
              if (randVal > migProb) {
                continue;
              }

              int demeEnd = std::min(demeStart + demeSize, containerSize());
              int worstInDeme = demeStart;
              gp_fitness worstFitness = std::numeric_limits<gp_fitness>::lowest();

              for (int idx = demeStart; idx < demeEnd; ++idx) {
                MyGP* candidate = NthMyGP(idx);
                if (!candidate) continue;
                gp_fitness fitness = candidate->getFitness();
                if (!std::isfinite(fitness) || fitness > worstFitness) {
                  worstFitness = fitness;
                  worstInDeme = idx;
                }
              }

              if (worstInDeme != bestOfPopulation) {
                MyGP* clone = (MyGP*)&best->duplicate();
                put(worstInDeme, *clone);
                demesPropagated++;
              }
            }

            bout << "# Propagated bakeoff winner (fitness="
                 << bestAggregatedFitness << ") to "
                 << demesPropagated << " demes" << endl;
            bout.flush();
          }
        }
      }

      // Use bakeoff winner for storage
      MyGP* eliteForStorage = best;
      if (!eliteForStorage) {
        eliteForStorage = NthMyGP(bestOfPopulation);
      }

      // Re-serialize the global elite for storage
      std::vector<char> updatedBuffer;
      boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> updatedOutStream(updatedBuffer);
      eliteForStorage->save(updatedOutStream);
      updatedOutStream.flush();
      bestOfEvalResults.gp = updatedBuffer;
      
      // now put the resulting elements into the S3 object
      Aws::S3::Model::PutObjectRequest request;
      request.SetBucket(cfg.s3Bucket);

      // path name is $base/RunDate/gen$gen.dmp
      request.SetKey(computedKeyName);

      std::ostringstream oss(std::ios::binary);
      boost::archive::binary_oarchive oa(oss);
      oa << bestOfEvalResults;

      std::shared_ptr<Aws::StringStream> ss = Aws::MakeShared<Aws::StringStream>("");
      *ss << oss.str();
      request.SetBody(ss);

      auto outcome = getS3Client()->PutObject(request);
      if (!outcome.IsSuccess()) {
        *logger.error() << "Error: " << outcome.GetError().GetMessage() << std::endl;
      }
    }
  }

  // Print (not mandatory)
  // virtual void printOn (ostream& os);

  // Access genetic programs (not mandatory)
  MyGP* NthMyGP(int n) {
    return (MyGP*)GPContainer::Nth(n);
  }
};

void MyPopulation::evaluate() {
  const AutocConfig& cfg = ConfigManager::getConfig();
  int scenarioCount = std::max<int>(generationScenarios.size(), 1);

  // Reset elite reeval flag for this generation (only flag first matching individual)
  gEliteReevalFlaggedThisGen = false;

  // Assign scenario index to each individual based on demetic grouping
  for (int n = 0; n < containerSize(); ++n) {
    MyGP* current = NthMyGP(n);
#if GPINTERNALCHECK
    if (!current) {
      std::cerr << "FATAL: MyPopulation::evaluate: Member of population is NULL" << std::endl; exit(1);
    }
#endif
    current->setScenarioIndex(computeScenarioIndexForIndividual(n));
  }

  // Evaluate all individuals using their assigned scenarios
  // In non-demetic mode: all individuals use scenario 0 (which contains all paths × all winds)
  // In demetic mode: each deme uses its own scenario (one path variant × all winds)
  GPPopulation::evaluate();
}

void MyGP::evaluate()
{
  tasks.push_back(this);
}

// Evaluate the fitness of a GP and save it into the class variable
// fitness.
void MyGP::evalTask(WorkerContext& context)
{
  stdFitness = 0;

  // TODO this save is probably cheaper when GP knows about boost archives...
  std::vector<char> buffer;
  boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> outStream(buffer);
  save(outStream);
  outStream.flush();

  const ScenarioDescriptor& scenario = scenarioForIndex(getScenarioIndex());
  uint64_t scenarioSequence = globalScenarioCounter.fetch_add(1, std::memory_order_relaxed) + 1;

  // Generate Lisp-form string of GP for exact comparison
  std::ostringstream gpStringStream;
  printOn(gpStringStream);
  std::string currentGpString = gpStringStream.str();

  EvalData evalData;
  evalData.gp = buffer;
  // Hash the Lisp string (deterministic) instead of binary serialization
  evalData.gpHash = hashString(currentGpString);
  evalData.isEliteReeval = false;  // Can be set by caller for crrcsim detailed logging
  evalData.pathList = scenario.pathList;
  evalData.scenario.scenarioSequence = scenarioSequence;
  evalData.scenario.bakeoffSequence = 0;
  evalData.scenarioList.clear();
  evalData.scenarioList.reserve(scenario.pathList.size());
  bool isBakeoff = bakeoffMode;
  bool enableLogging = enableDeterministicTestLogging.load(std::memory_order_relaxed);

  const AutocConfig& cfg = ConfigManager::getConfig();

  // Build scenario list FIRST before checking for elite re-eval
  for (size_t idx = 0; idx < scenario.pathList.size(); ++idx) {
    ScenarioMetadata meta;
    meta.scenarioSequence = scenarioSequence;
    meta.enableDeterministicLogging = enableLogging;

    if (isBakeoff) {
      meta.bakeoffSequence = bakeoffPathCounter.fetch_add(1, std::memory_order_relaxed) + 1;
    } else {
      meta.bakeoffSequence = 0;
    }

    if (cfg.demeticGrouping && cfg.demeSize > 0) {
      // Demetic mode: scenario has one path variant across multiple winds
      meta.pathVariantIndex = scenario.pathVariantIndex;
      if (idx < scenario.windScenarios.size()) {
        meta.windVariantIndex = scenario.windScenarios[idx].windVariantIndex;
        meta.windSeed = scenario.windScenarios[idx].windSeed;
      } else {
        meta.windVariantIndex = scenario.windVariantIndex;
        meta.windSeed = scenario.windSeed;
      }
    } else {
      // Non-demetic mode: scenario has all paths × all winds
      // pathList is organized in PATH-MAJOR order: [p0w0, p0w1, ..., p0w5, p1w0, p1w1, ...]
      size_t numWindScenarios = scenario.windScenarios.size();
      if (numWindScenarios > 0 && scenario.pathList.size() > 0) {
        size_t numBasePaths = scenario.pathList.size() / numWindScenarios;
        size_t pathIdx = idx / numWindScenarios;
        size_t windIdx = idx % numWindScenarios;

        meta.pathVariantIndex = static_cast<int>(pathIdx);
        if (windIdx < scenario.windScenarios.size()) {
          meta.windVariantIndex = scenario.windScenarios[windIdx].windVariantIndex;
          meta.windSeed = scenario.windScenarios[windIdx].windSeed;
        } else {
          meta.windVariantIndex = scenario.windVariantIndex;
          meta.windSeed = scenario.windSeed;
        }
      } else {
        // Fallback if structure is unexpected
        meta.pathVariantIndex = static_cast<int>(idx);
        meta.windVariantIndex = scenario.windVariantIndex;
        meta.windSeed = scenario.windSeed;
      }
    }

    // VARIATIONS1: Populate entry/wind variation offsets from windSeed
    populateVariationOffsets(meta);

    // Variable rabbit speed: Apply speed profile to path timing (uses pre-fetched profile)
    applySpeedProfileToPath(evalData.pathList[idx], meta.windVariantIndex);

    evalData.scenarioList.push_back(meta);
  }
  if (!evalData.scenarioList.empty()) {
    evalData.scenario = evalData.scenarioList.front();
  } else {
    evalData.scenario.scenarioSequence = scenarioSequence;
    evalData.scenario.bakeoffSequence = isBakeoff
      ? bakeoffPathCounter.fetch_add(1, std::memory_order_relaxed) + 1
      : 0;
  }

  evalData.sanitizePaths();

  // Check if this is a re-evaluation of the elite from previous generation
  // If so, flag it for detailed physics trace collection
  if (gEliteTraceEnabled && evalData.gpHash != 0 && evalData.gpHash == gLastEliteGpHash) {
    evalData.isEliteReeval = true;
  }

  sendRPC(*context.socket, evalData);

  // How did it go?
  context.evalResults = receiveRPC<EvalResults>(*context.socket);

  // Verify GP payload integrity (dispatcher vs worker)
  if (context.evalResults.gpHash == 0 && !context.evalResults.gp.empty()) {
    context.evalResults.gpHash = hashByteVector(context.evalResults.gp);
  }
  if (evalData.gpHash != 0 && context.evalResults.gpHash != 0 &&
      evalData.gpHash != context.evalResults.gpHash) {
    *logger.warn() << "[AUTOC_GP_HASH_MISMATCH] expected=0x"
                   << std::hex << evalData.gpHash
                   << " got=0x" << context.evalResults.gpHash
                   << std::dec
                   << " workerId=" << context.evalResults.workerId
                   << " workerPid=" << context.evalResults.workerPid
                   << " evalCounter=" << context.evalResults.workerEvalCounter
                   << " size=" << evalData.gp.size()
                   << endl;
  }


  if (context.evalResults.scenario.windSeed == 0 && evalData.scenario.windSeed != 0) {
    context.evalResults.scenario = evalData.scenario;
  }
  if (context.evalResults.pathList.empty() && !evalData.pathList.empty()) {
    context.evalResults.pathList = evalData.pathList;
  }
  if (context.evalResults.scenarioList.empty() && !evalData.scenarioList.empty()) {
    context.evalResults.scenarioList = evalData.scenarioList;
  }
  globalSimRunCounter.fetch_add(context.evalResults.pathList.size(), std::memory_order_relaxed);
  // context.evalResults.dump(std::cerr);

  if (printEval) {
    boost::unique_lock<boost::mutex> lock(evalCollectorMutex);
    if (activeEvalCollector) {
      *logger.debug() << "Incoming eval results has "
                      << context.evalResults.pathList.size() << " paths" << endl;
      appendEvalResults(*activeEvalCollector, context.evalResults);
      *logger.debug() << "Collector now has " << activeEvalCollector->pathList.size()
                      << " paths" << endl;
    } else {
      bestOfEvalResults = context.evalResults;
    }
  }

  // If this was an elite reeval with trace, store for comparison at generation end
  // (Don't try to preserve in bestOfEvalResults - it gets overwritten by subsequent evals)
  if (evalData.isEliteReeval && !context.evalResults.physicsTrace.empty()) {
    boost::unique_lock<boost::mutex> lock(evalCollectorMutex);

    // Store in pending - survives subsequent bestOfEvalResults overwrites
    // Gets moved to gLastEliteEvalResults at end of generation
    gPendingEliteEvalResults.physicsTrace = std::move(context.evalResults.physicsTrace);
    gPendingEliteEvalResults.debugSamples = std::move(context.evalResults.debugSamples);
    gPendingEliteEvalResults.workerId = context.evalResults.workerId;
    gPendingEliteEvalResults.workerPid = context.evalResults.workerPid;
    gPendingEliteEvalResults.workerEvalCounter = context.evalResults.workerEvalCounter;
  }

  // Compute fitness for each path and sum
  // Crash penalty (soft lexicographic) applied per-path

  for (int i = 0; i < context.evalResults.pathList.size(); i++) {
    bool printHeader = true;

    // get path, actual and aircraft state
    auto& path = context.evalResults.pathList.at(i);
    auto& aircraftState = context.evalResults.aircraftStateList.at(i);
    auto& crashReason = context.evalResults.crashReasonList.at(i);

    if (path.empty() || aircraftState.empty()) {
      continue;
    }

    // Get scenario metadata for this path (for logging)
    int pathVariantIndex = i;  // default to loop index
    int windVariantIndex = -1;
    unsigned int windSeed = 0;
    uint64_t scenarioSequence = context.evalResults.scenario.scenarioSequence;
    uint64_t bakeoffSequence = context.evalResults.scenario.bakeoffSequence;
    if (i < context.evalResults.scenarioList.size()) {
      const auto& meta = context.evalResults.scenarioList.at(i);
      pathVariantIndex = meta.pathVariantIndex;
      windVariantIndex = meta.windVariantIndex;
      windSeed = meta.windSeed;
      scenarioSequence = meta.scenarioSequence;
      bakeoffSequence = meta.bakeoffSequence;
    } else {
      if (windVariantIndex < 0) {
        windVariantIndex = context.evalResults.scenario.windVariantIndex;
      }
      if (windSeed == 0) {
        windSeed = context.evalResults.scenario.windSeed;
      }
    }

    // ========================================================================
    // SIMPLIFIED FITNESS (v3): Two objectives - distance + attitude delta
    // See specs/FITNESS_SIMPLIFY_20260221.md
    // ========================================================================
    gp_fitness localFitness = 0.0;
    int stepIndex = 0;

    // Accumulate raw values, scale attitude at end using path geometry
    gp_fitness distance_sum = 0.0;
    gp_fitness attitude_sum = 0.0;
    int simulation_steps = 0;

    // Intercept-budget scaling: compute budget from initial displacement + heading
    double interceptBudget = 0.0;
    {
      gp_vec3 initialPos = aircraftState.at(0).getPosition();
      gp_vec3 pathStart = path.at(0).start;
      gp_vec3 offset3d = initialPos - pathStart;
      double displacement = sqrt(offset3d[0] * offset3d[0] + offset3d[1] * offset3d[1]);
      double headingOffset = 0.0;
      if (i < context.evalResults.scenarioList.size()) {
        headingOffset = context.evalResults.scenarioList.at(i).entryHeadingOffset;
      }
      interceptBudget = computeInterceptBudget(displacement, headingOffset,
                                                SIM_INITIAL_VELOCITY, gRabbitSpeedConfig.nominal);
    }

    // Attitude delta tracking
    gp_fitness prev_roll = 0.0, prev_pitch = 0.0;
    bool first_attitude_sample = true;

    // now walk next steps of actual path
    while (++stepIndex < aircraftState.size()) {
      auto& stepAircraftState = aircraftState.at(stepIndex);
      int rawPathIndex = stepAircraftState.getThisPathIndex();
      int pathIndex = std::clamp(rawPathIndex, 0, static_cast<int>(path.size()) - 1);
      const Path& currentPathPoint = path.at(pathIndex);

      gp_vec3 aircraftPosition = stepAircraftState.getPosition();
      gp_quat craftOrientation = stepAircraftState.getOrientation();

      // ====================================================================
      // ATTITUDE DELTA: Per-step roll/pitch change (radians)
      // ====================================================================
      gp_fitness qw = craftOrientation.w();
      gp_fitness qx = craftOrientation.x();
      gp_fitness qy = craftOrientation.y();
      gp_fitness qz = craftOrientation.z();

      // Roll (phi): rotation around X axis
      gp_fitness roll = atan2(2.0 * (qw * qx + qy * qz),
                              1.0 - 2.0 * (qx * qx + qy * qy));
      // Pitch (theta): rotation around Y axis
      gp_fitness sinp = 2.0 * (qw * qy - qz * qx);
      sinp = std::clamp(sinp, -1.0, 1.0);
      gp_fitness pitch = asin(sinp);

      gp_fitness attitude_delta = 0.0;
      if (first_attitude_sample) {
        prev_roll = roll;
        prev_pitch = pitch;
        first_attitude_sample = false;
      } else {
        attitude_delta = fabs(roll - prev_roll) + fabs(pitch - prev_pitch);
        prev_roll = roll;
        prev_pitch = pitch;
      }

      // ====================================================================
      // DISTANCE: 3D distance to rabbit (meters)
      // ====================================================================
      gp_vec3 craftOffset = aircraftPosition - currentPathPoint.start;
      gp_scalar distance = craftOffset.norm();

      // Accumulate with normalized nonlinear penalty per step
      // Intercept-budget scaling: ramp penalty from floor→ceiling over budget
      double stepTimeSec = stepAircraftState.getSimTimeMsec() / 1000.0;
      double interceptScale = computeInterceptScale(stepTimeSec, interceptBudget);
      double distDelta = fabs(distance - DISTANCE_TARGET);
      distance_sum += pow(interceptScale * distDelta / DISTANCE_NORM, DISTANCE_POWER);
      attitude_sum += pow(interceptScale * attitude_delta / ATTITUDE_NORM, ATTITUDE_POWER);
      simulation_steps++;

      // Per-step logging to data.dat
      if (printEval) {
        if (printHeader) {
          fout << "Scn    Bake   Pth/Wnd:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power      qw      qx      qy      qz   vxBody   vyBody   vzBody    alpha     beta   dtheta    dphi   dhome     dist   attDlt  rabVel intScl\n";
          printHeader = false;
        }

        // Compute rabbit velocity from path position delta
        static gp_vec3 prevPathPos(0, 0, 0);
        static long prevPathTime = 0;
        static bool firstStep = true;
        gp_scalar rabbitVel = 0.0f;
        gp_vec3 currPathPos = path.at(pathIndex).start;
        long currPathTime = path.at(pathIndex).simTimeMsec;
        if (simulation_steps == 1) {
          firstStep = true;
        }
        if (!firstStep && currPathTime > prevPathTime) {
          gp_scalar dist = (currPathPos - prevPathPos).norm();
          gp_scalar dt = static_cast<gp_scalar>(currPathTime - prevPathTime) / 1000.0f;
          rabbitVel = dist / dt;
        }
        prevPathPos = currPathPos;
        prevPathTime = currPathTime;
        firstStep = false;

        // Convert orientation to euler angles
        Eigen::Matrix3f rotMatrix = craftOrientation.toRotationMatrix();
        gp_vec3 euler;
        if (std::abs(rotMatrix(2, 0)) > static_cast<gp_scalar>(0.99999f)) {
          euler[0] = 0;
          if (rotMatrix(2, 0) > 0) {
            euler[1] = -static_cast<gp_scalar>(M_PI / 2.0);
            euler[2] = -atan2(rotMatrix(1, 2), rotMatrix(0, 2));
          } else {
            euler[1] = static_cast<gp_scalar>(M_PI / 2.0);
            euler[2] = atan2(rotMatrix(1, 2), rotMatrix(0, 2));
          }
        } else {
          euler[0] = atan2(rotMatrix(2, 1), rotMatrix(2, 2));
          euler[1] = -asin(rotMatrix(2, 0));
          euler[2] = atan2(rotMatrix(1, 0), rotMatrix(0, 0));
        }

        // Body-frame velocity
        gp_vec3 velocity_body = stepAircraftState.getOrientation().inverse() *
                                stepAircraftState.getVelocity();

        // GP operator values
        gp_scalar alpha_deg = atan2(-velocity_body.z(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);
        gp_scalar beta_deg = atan2(velocity_body.y(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);

        // Angle to target in body frame
        gp_vec3 craftToTarget = path.at(pathIndex).start - stepAircraftState.getPosition();
        gp_vec3 target_body = stepAircraftState.getOrientation().inverse() * craftToTarget;
        gp_scalar dtheta_deg = atan2(-target_body.z(), target_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);
        gp_scalar dphi_deg = atan2(target_body.y(), -target_body.z()) * static_cast<gp_scalar>(180.0 / M_PI);

        // Distance to home
        gp_vec3 home(0, 0, SIM_INITIAL_ALTITUDE);
        gp_scalar dhome = (home - stepAircraftState.getPosition()).norm();

        gp_quat q = craftOrientation;

        char outbuf[1600];
        sprintf(outbuf, "%06llu %06llu %03d/%02d:%04d: %06ld %3d % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 7.4f % 7.4f % 7.4f % 7.4f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.3f % 6.4f % 6.1f %5.3f\n",
          static_cast<unsigned long long>(scenarioSequence),
          static_cast<unsigned long long>(bakeoffSequence),
          pathVariantIndex, windVariantIndex, simulation_steps,
          stepAircraftState.getSimTimeMsec(), pathIndex,
          path.at(pathIndex).distanceFromStart,
          path.at(pathIndex).start[0],
          path.at(pathIndex).start[1],
          path.at(pathIndex).start[2],
          stepAircraftState.getPosition()[0],
          stepAircraftState.getPosition()[1],
          stepAircraftState.getPosition()[2],
          euler[2], euler[1], euler[0],
          stepAircraftState.getRelVel(),
          stepAircraftState.getRollCommand(),
          stepAircraftState.getPitchCommand(),
          stepAircraftState.getThrottleCommand(),
          q.w(), q.x(), q.y(), q.z(),
          velocity_body.x(), velocity_body.y(), velocity_body.z(),
          alpha_deg, beta_deg,
          dtheta_deg, dphi_deg,
          dhome,
          distance,        // distance to rabbit
          attitude_delta,  // attitude change this step
          rabbitVel,
          interceptScale   // intercept-budget scaling factor
        );
        fout << outbuf;
      }
    }

    // Scale attitude by path geometry so 1 rad excess ≈ (path_dist/path_turn) meters
    gp_scalar path_distance = path.back().distanceFromStart;
    gp_scalar path_turn_rad = path.back().radiansFromStart;
    // Fallback for straight paths: treat as one full rotation worth
    gp_scalar attitude_scale = path_distance / std::max(path_turn_rad, static_cast<gp_scalar>(2.0 * M_PI));

    localFitness = distance_sum + attitude_sum * attitude_scale;

    if (isnan(localFitness)) {
      nanDetector++;
    }

    // Crash penalty: soft lexicographic - completion dominates, quality provides gradient
    if (crashReason != CrashReason::None) {
      gp_fitness total_path_distance = path.back().distanceFromStart;
      gp_fitness fraction_completed =
        static_cast<gp_fitness>(path.at(aircraftState.back().getThisPathIndex()).distanceFromStart) / total_path_distance;
      fraction_completed = std::max(fraction_completed, static_cast<gp_fitness>(0.001));
      gp_fitness completion_penalty = (1.0 - fraction_completed) * CRASH_COMPLETION_WEIGHT;
      localFitness = completion_penalty + localFitness;
    }

    // Sum path fitness
    stdFitness += localFitness;
  }

  // Mark fitness as valid for the GP library
  fitnessValid = 1;
}

void newHandler()
{
  cerr << "\nFatal error: Out of memory." << endl;
  exit(1);
}


// Get compiled topology as std::vector (needed by NNPopulation API)
static std::vector<int> getCompiledTopology() {
  return std::vector<int>(NN_TOPOLOGY, NN_TOPOLOGY + NN_NUM_LAYERS);
}

// Compute fitness for an NN individual from EvalResults
// Same formula as MyGP::evalTask() lines 1210-1417
static double computeNNFitness(EvalResults& evalResults) {
  double totalFitness = 0.0;

  for (size_t i = 0; i < evalResults.pathList.size(); i++) {
    auto& path = evalResults.pathList.at(i);
    auto& aircraftStates = evalResults.aircraftStateList.at(i);
    auto& crashReason = evalResults.crashReasonList.at(i);

    if (path.empty() || aircraftStates.empty()) continue;

    double distance_sum = 0.0;
    double attitude_sum = 0.0;
    int simulation_steps = 0;

    // Intercept-budget scaling
    double interceptBudget = 0.0;
    {
      gp_vec3 initialPos = aircraftStates.at(0).getPosition();
      gp_vec3 pathStart = path.at(0).start;
      gp_vec3 offset3d = initialPos - pathStart;
      double displacement = sqrt(offset3d[0] * offset3d[0] + offset3d[1] * offset3d[1]);
      double headingOffset = 0.0;
      if (i < evalResults.scenarioList.size()) {
        headingOffset = evalResults.scenarioList.at(i).entryHeadingOffset;
      }
      interceptBudget = computeInterceptBudget(displacement, headingOffset,
                                                SIM_INITIAL_VELOCITY, gRabbitSpeedConfig.nominal);
    }

    // Attitude delta tracking
    double prev_roll = 0.0, prev_pitch = 0.0;
    bool first_attitude_sample = true;
    int stepIndex = 0;

    while (++stepIndex < static_cast<int>(aircraftStates.size())) {
      auto& stepState = aircraftStates.at(stepIndex);
      int pathIndex = std::clamp(stepState.getThisPathIndex(), 0, static_cast<int>(path.size()) - 1);

      gp_vec3 aircraftPosition = stepState.getPosition();
      gp_quat craftOrientation = stepState.getOrientation();

      // Attitude delta
      double qw = craftOrientation.w(), qx = craftOrientation.x();
      double qy = craftOrientation.y(), qz = craftOrientation.z();
      double roll = atan2(2.0 * (qw * qx + qy * qz), 1.0 - 2.0 * (qx * qx + qy * qy));
      double sinp = std::clamp(2.0 * (qw * qy - qz * qx), -1.0, 1.0);
      double pitch = asin(sinp);

      double attitude_delta = 0.0;
      if (first_attitude_sample) {
        prev_roll = roll;
        prev_pitch = pitch;
        first_attitude_sample = false;
      } else {
        attitude_delta = fabs(roll - prev_roll) + fabs(pitch - prev_pitch);
        prev_roll = roll;
        prev_pitch = pitch;
      }

      // Distance to rabbit
      gp_vec3 craftOffset = aircraftPosition - path.at(pathIndex).start;
      double distance = craftOffset.norm();

      // Accumulate penalties
      double stepTimeSec = stepState.getSimTimeMsec() / 1000.0;
      double interceptScale = computeInterceptScale(stepTimeSec, interceptBudget);
      double distDelta = fabs(distance - DISTANCE_TARGET);
      distance_sum += pow(interceptScale * distDelta / DISTANCE_NORM, DISTANCE_POWER);
      attitude_sum += pow(interceptScale * attitude_delta / ATTITUDE_NORM, ATTITUDE_POWER);
      simulation_steps++;
    }

    // Scale attitude by path geometry
    double path_distance = path.back().distanceFromStart;
    double path_turn_rad = path.back().radiansFromStart;
    double attitude_scale = path_distance / std::max(path_turn_rad, 2.0 * M_PI);

    double localFitness = distance_sum + attitude_sum * attitude_scale;

    // Crash penalty
    if (crashReason != CrashReason::None) {
      double total_path_distance = path.back().distanceFromStart;
      double fraction_completed =
        static_cast<double>(path.at(aircraftStates.back().getThisPathIndex()).distanceFromStart) / total_path_distance;
      fraction_completed = std::max(fraction_completed, 0.001);
      double completion_penalty = (1.0 - fraction_completed) * CRASH_COMPLETION_WEIGHT;
      localFitness = completion_penalty + localFitness;
    }

    totalFitness += localFitness;
  }

  return totalFitness;
}

// Log per-step data from EvalResults to data.dat
// NN mode: actual NN inputs (normalized) and outputs captured in minisim, then diagnostics
static void logEvalResults(std::ofstream& fout, EvalResults& results) {
  bool printHeader = true;

  for (size_t i = 0; i < results.pathList.size(); i++) {
    auto& path = results.pathList.at(i);
    auto& aircraftStates = results.aircraftStateList.at(i);
    if (path.empty() || aircraftStates.empty()) continue;

    uint64_t scenarioSequence = 0;
    uint64_t bakeoffSequence = 0;
    int pathVariantIndex = 0;
    int windVariantIndex = 0;
    if (i < results.scenarioList.size()) {
      scenarioSequence = results.scenarioList.at(i).scenarioSequence;
      bakeoffSequence = results.scenarioList.at(i).bakeoffSequence;
      pathVariantIndex = results.scenarioList.at(i).pathVariantIndex;
      windVariantIndex = results.scenarioList.at(i).windVariantIndex;
    }

    // Intercept-budget scaling
    double interceptBudget = 0.0;
    {
      gp_vec3 initialPos = aircraftStates.at(0).getPosition();
      gp_vec3 pathStart = path.at(0).start;
      gp_vec3 offset3d = initialPos - pathStart;
      double displacement = sqrt(offset3d[0] * offset3d[0] + offset3d[1] * offset3d[1]);
      double headingOffset = 0.0;
      if (i < results.scenarioList.size()) {
        headingOffset = results.scenarioList.at(i).entryHeadingOffset;
      }
      interceptBudget = computeInterceptBudget(displacement, headingOffset,
                                                SIM_INITIAL_VELOCITY, gRabbitSpeedConfig.nominal);
    }

    double prev_roll = 0.0, prev_pitch = 0.0;
    bool first_attitude_sample = true;
    int simulation_steps = 0;

    gp_vec3 prevPathPos(0, 0, 0);
    long prevPathTime = 0;
    bool firstStep = true;

    int stepIndex = 0;
    while (++stepIndex < static_cast<int>(aircraftStates.size())) {
      auto& stepState = aircraftStates.at(stepIndex);
      int pathIndex = std::clamp(stepState.getThisPathIndex(), 0, static_cast<int>(path.size()) - 1);

      gp_vec3 aircraftPosition = stepState.getPosition();
      gp_quat craftOrientation = stepState.getOrientation();

      // Attitude delta
      double aqw = craftOrientation.w(), aqx = craftOrientation.x();
      double aqy = craftOrientation.y(), aqz = craftOrientation.z();
      double roll = atan2(2.0 * (aqw * aqx + aqy * aqz), 1.0 - 2.0 * (aqx * aqx + aqy * aqy));
      double sinp = std::clamp(2.0 * (aqw * aqy - aqz * aqx), -1.0, 1.0);
      double pitch = asin(sinp);

      double attitude_delta = 0.0;
      if (first_attitude_sample) {
        prev_roll = roll;
        prev_pitch = pitch;
        first_attitude_sample = false;
      } else {
        attitude_delta = fabs(roll - prev_roll) + fabs(pitch - prev_pitch);
        prev_roll = roll;
        prev_pitch = pitch;
      }

      // Distance to rabbit
      gp_vec3 craftOffset = aircraftPosition - path.at(pathIndex).start;
      double distance = craftOffset.norm();

      double stepTimeSec = stepState.getSimTimeMsec() / 1000.0;
      double interceptScale = computeInterceptScale(stepTimeSec, interceptBudget);
      simulation_steps++;

      // Rabbit velocity
      gp_scalar rabbitVel = 0.0f;
      gp_vec3 currPathPos = path.at(pathIndex).start;
      long currPathTime = path.at(pathIndex).simTimeMsec;
      if (simulation_steps == 1) {
        firstStep = true;
      }
      if (!firstStep && currPathTime > prevPathTime) {
        gp_scalar dist = (currPathPos - prevPathPos).norm();
        gp_scalar dt = static_cast<gp_scalar>(currPathTime - prevPathTime) / 1000.0f;
        rabbitVel = dist / dt;
      }
      prevPathPos = currPathPos;
      prevPathTime = currPathTime;
      firstStep = false;

      // Body-frame velocity
      gp_vec3 velocity_body = stepState.getOrientation().inverse() * stepState.getVelocity();

      // Distance to home
      gp_vec3 home(0, 0, SIM_INITIAL_ALTITUDE);
      gp_scalar dhome = (home - stepState.getPosition()).norm();

      // Per-step logging — all % 7.4f fields are 8 chars, % 8.2f fields are 9 chars
      if (printHeader) {
        fout << "Scn    Bake   Pth/Wnd:Step:  Time Idx"  // 37 chars
             << "   dPhi0   dPhi1   dPhi3   dPhi9"       // 4 × 8 = 32
             << "   dTht0   dTht1   dTht3   dTht9"       // 4 × 8 = 32
             << "    dst0    dst1    dst3    dst9"        // 4 × 8 = 32
             << "      qw      qx      qy      qz"      // 4 × 8 = 32
             << "     vel   alpha    beta"                // 3 × 8 = 24
             << "   outPt   outRl   outTh"                // 3 × 8 = 24
             << "    cmdP    cmdR    cmdT"                // 3 × 8 = 24
             << "    pathX    pathY    pathZ"              // 3 × 9 = 27
             << "        X        Y        Z"             // 3 × 9 = 27
             << "   vxBody   vyBody   vzBody"             // 3 × 9 = 27
             << "    dhome     dist  attDlt   rabVl   intSc" // 9+9+8+8+8=42
             << "\n";
        printHeader = false;
      }

      // Actual NN inputs/outputs captured in minisim
      const float* in = stepState.getNNInputs();
      const float* out = stepState.getNNOutputs();

      char outbuf[2048];
      sprintf(outbuf,
        "%06llu %06llu %03d/%02d:%04d: %06ld %3d"
        " % 7.4f % 7.4f % 7.4f % 7.4f"    // dPhi[0,1,3,9]
        " % 7.4f % 7.4f % 7.4f % 7.4f"    // dTheta[0,1,3,9]
        " % 7.4f % 7.4f % 7.4f % 7.4f"    // dist[0,1,3,9]
        " % 7.4f % 7.4f % 7.4f % 7.4f"    // qw,qx,qy,qz
        " % 7.4f % 7.4f % 7.4f"            // vel, alpha, beta
        " % 7.4f % 7.4f % 7.4f"            // NN outputs: pitch, roll, throttle
        " % 7.4f % 7.4f % 7.4f"            // cmd feedback: pitch, roll, throttle
        " % 8.2f % 8.2f % 8.2f"            // pathX, pathY, pathZ
        " % 8.2f % 8.2f % 8.2f"            // X, Y, Z
        " % 8.2f % 8.2f % 8.2f"            // vxBody, vyBody, vzBody
        " % 8.2f % 8.3f % 7.4f % 7.1f % 7.3f"  // dhome, dist, attDlt, rabVel, intScl
        "\n",
        static_cast<unsigned long long>(scenarioSequence),
        static_cast<unsigned long long>(bakeoffSequence),
        pathVariantIndex, windVariantIndex, simulation_steps,
        stepState.getSimTimeMsec(), pathIndex,
        // NN inputs [0-21] — actual normalized values presented to NN
        in[0], in[1], in[2], in[3],         // dPhi temporal history
        in[4], in[5], in[6], in[7],         // dTheta temporal history
        in[8], in[9], in[10], in[11],       // distance temporal history
        in[12], in[13], in[14], in[15],     // quaternion attitude
        in[16], in[17], in[18],             // vel, alpha, beta
        // NN outputs [0-2] — actual tanh outputs
        out[0], out[1], out[2],             // pitch, roll, throttle
        // cmd feedback inputs — previous tick's commands (P, R, T order)
        in[19], in[20], in[21],             // pitch, roll, throttle
        // Diagnostics
        path.at(pathIndex).start[0],
        path.at(pathIndex).start[1],
        path.at(pathIndex).start[2],
        stepState.getPosition()[0],
        stepState.getPosition()[1],
        stepState.getPosition()[2],
        velocity_body.x(), velocity_body.y(), velocity_body.z(),
        dhome,
        static_cast<gp_scalar>(distance),
        static_cast<gp_scalar>(attitude_delta),
        rabbitVel,
        static_cast<gp_scalar>(interceptScale)
      );
      fout << outbuf;
    }
  }
  fout.flush();
}

// NN evaluation mode: load weight file, run through scenarios, report fitness
static void runNNEvaluation(
    const std::string& startTime,
    std::ofstream& fout,
    std::ofstream& bout
) {
  const AutocConfig& cfg = ConfigManager::getConfig();

  *logger.info() << "NN Evaluation mode" << endl;
  *logger.info() << "  Weight file: " << cfg.nnWeightFile << endl;

  // Load NN weight file
  std::ifstream weightFile(cfg.nnWeightFile, std::ios::binary | std::ios::ate);
  if (!weightFile.is_open()) {
    *logger.error() << "Cannot open NN weight file: " << cfg.nnWeightFile << endl;
    exit(1);
  }
  std::streamsize fileSize = weightFile.tellg();
  weightFile.seekg(0, std::ios::beg);
  std::vector<uint8_t> fileData(fileSize);
  if (!weightFile.read(reinterpret_cast<char*>(fileData.data()), fileSize)) {
    *logger.error() << "Error reading NN weight file" << endl;
    exit(1);
  }
  weightFile.close();

  if (!nn_detect_format(fileData.data(), fileData.size())) {
    *logger.error() << "Weight file is not in NN01 format" << endl;
    exit(1);
  }

  NNGenome genome;
  if (!nn_deserialize(fileData.data(), fileData.size(), genome)) {
    *logger.error() << "Failed to deserialize NN genome" << endl;
    exit(1);
  }

  {
    std::ostringstream topo;
    for (size_t i = 0; i < genome.topology.size(); i++) {
      if (i > 0) topo << " -> ";
      topo << genome.topology[i];
    }
    *logger.info() << "  Topology: " << topo.str() << " (" << genome.weights.size() << " weights)" << endl;
  }
  *logger.info() << "  Stored fitness: " << std::fixed << std::setprecision(6) << genome.fitness << endl;
  *logger.info() << "  Stored generation: " << genome.generation << endl;

  // Serialize genome for RPC (same format minisim expects)
  std::vector<uint8_t> nnData;
  nn_serialize(genome, nnData);

  // Use scenario 0 (same as training)
  const ScenarioDescriptor& scenario = scenarioForIndex(0);
  uint64_t scenarioSequence = globalScenarioCounter.fetch_add(1, std::memory_order_relaxed) + 1;

  EvalData evalData;
  evalData.controllerType = ControllerType::NEURAL_NET;
  evalData.gp.assign(reinterpret_cast<const char*>(nnData.data()),
                     reinterpret_cast<const char*>(nnData.data() + nnData.size()));
  evalData.gpHash = hashByteVector(evalData.gp);
  evalData.isEliteReeval = false;
  evalData.pathList = scenario.pathList;
  evalData.scenario.scenarioSequence = scenarioSequence;
  evalData.scenario.bakeoffSequence = 0;

  evalData.scenarioList.clear();
  evalData.scenarioList.reserve(scenario.pathList.size());
  for (size_t idx = 0; idx < scenario.pathList.size(); ++idx) {
    ScenarioMetadata meta;
    meta.scenarioSequence = scenarioSequence;
    meta.enableDeterministicLogging = false;
    meta.bakeoffSequence = 0;
    size_t numWindScenarios = scenario.windScenarios.size();
    if (numWindScenarios > 0 && scenario.pathList.size() > 0) {
      size_t pathIdx = idx / numWindScenarios;
      size_t windIdx = idx % numWindScenarios;
      meta.pathVariantIndex = static_cast<int>(pathIdx);
      if (windIdx < scenario.windScenarios.size()) {
        meta.windVariantIndex = scenario.windScenarios[windIdx].windVariantIndex;
        meta.windSeed = scenario.windScenarios[windIdx].windSeed;
      }
    } else {
      meta.pathVariantIndex = static_cast<int>(idx);
    }
    populateVariationOffsets(meta);
    applySpeedProfileToPath(evalData.pathList[idx], meta.windVariantIndex);
    evalData.scenarioList.push_back(meta);
  }
  if (!evalData.scenarioList.empty()) {
    evalData.scenario = evalData.scenarioList.front();
  }
  evalData.sanitizePaths();

  // Send to minisim and get results
  auto evalDataPtr = std::make_shared<EvalData>(std::move(evalData));
  EvalResults evalResults;
  threadPool->enqueue([evalDataPtr, &evalResults](WorkerContext& context) {
    sendRPC(*context.socket, *evalDataPtr);
    evalResults = receiveRPC<EvalResults>(*context.socket);
  });
  threadPool->wait_for_tasks();

  // Compute fitness
  double fitness = computeNNFitness(evalResults);

  // Log per-step data to data.dat
  logEvalResults(fout, evalResults);

  *logger.info() << "NN Eval fitness: " << std::fixed << std::setprecision(6) << fitness << endl;
  *logger.info() << "Stored fitness:  " << std::fixed << std::setprecision(6) << genome.fitness << endl;

  // Log to statistics file
  bout << "#NNEval fitness=" << std::fixed << std::setprecision(6) << fitness
       << " storedFitness=" << genome.fitness
       << " weightFile=" << cfg.nnWeightFile
       << " scenarios=" << evalResults.pathList.size()
       << endl;
}

// NN evolution main loop — runs when ControllerType=NN
static void runNNEvolution(
    const std::string& startTime,
    const std::chrono::steady_clock::time_point& runStartTime,
    std::ofstream& fout,
    std::ofstream& bout,
    const std::function<void(int)>& logGenerationStats
) {
  const AutocConfig& cfg = ConfigManager::getConfig();

  std::vector<int> topology = getCompiledTopology();
  int popSize = cfg.populationSize;
  int numGens = cfg.numberOfGenerations;

  *logger.info() << "NN Evolution mode" << endl;
  *logger.info() << "  Topology: " << NN_TOPOLOGY_STRING
                 << " (" << NN_WEIGHT_COUNT << " weights)" << endl;
  *logger.info() << "  Population: " << popSize << endl;
  *logger.info() << "  Generations: " << numGens << endl;
  *logger.info() << "  MutationSigma: " << cfg.nnMutationSigma << endl;
  *logger.info() << "  CrossoverAlpha: " << cfg.nnCrossoverAlpha << endl;
  *logger.info() << "  TournamentSize: " << cfg.tournamentSize << endl;
  *logger.info() << "  CrossoverProb: " << cfg.crossoverProbability << "%" << endl;
  *logger.info() << "  CreationProb: " << cfg.creationProbability << "%" << endl;
  *logger.info() << "  MutationOnlyProb: " << cfg.swapMutationProbability << "%" << endl;
  *logger.info() << "  Elitism: " << cfg.addBestToNewPopulation << endl;

  // Initialize population
  NNPopulation pop;
  nn_init_population(pop, topology, popSize);

  // Set initial mutation sigma from config
  for (auto& ind : pop.individuals) {
    ind.mutation_sigma = static_cast<float>(cfg.nnMutationSigma);
  }

  *logger.info() << "Population initialized." << endl;

  // Generation loop (start at 1 like GP — gen0 is Xavier baseline, not stored)
  for (int gen = 1; gen <= numGens; gen++) {

    // RAMP_LANDSCAPE: Update current generation for variation scaling
    gCurrentGeneration = gen;

    // Generate paths for this generation
    generationPaths = generateSmoothPaths(const_cast<char*>(cfg.generatorMethod.c_str()),
                                          cfg.simNumPathsPerGen,
                                          SIM_PATH_BOUNDS, SIM_PATH_BOUNDS,
                                          gPathSeed);
    rebuildGenerationScenarios(generationPaths);

    // Evaluate each individual
    for (int ind = 0; ind < popSize; ind++) {
      NNGenome& genome = pop.individuals[ind];

      // Serialize genome to binary
      std::vector<uint8_t> nnData;
      nn_serialize(genome, nnData);

      // Build EvalData
      const ScenarioDescriptor& scenario = scenarioForIndex(ind % generationScenarios.size());
      uint64_t scenarioSequence = globalScenarioCounter.fetch_add(1, std::memory_order_relaxed) + 1;

      EvalData evalData;
      evalData.controllerType = ControllerType::NEURAL_NET;
      evalData.gp.assign(reinterpret_cast<const char*>(nnData.data()),
                         reinterpret_cast<const char*>(nnData.data() + nnData.size()));
      evalData.gpHash = hashByteVector(evalData.gp);
      evalData.isEliteReeval = false;
      evalData.pathList = scenario.pathList;
      evalData.scenario.scenarioSequence = scenarioSequence;
      evalData.scenario.bakeoffSequence = 0;

      // Build scenario list
      evalData.scenarioList.clear();
      evalData.scenarioList.reserve(scenario.pathList.size());
      for (size_t idx = 0; idx < scenario.pathList.size(); ++idx) {
        ScenarioMetadata meta;
        meta.scenarioSequence = scenarioSequence;
        meta.enableDeterministicLogging = false;
        meta.bakeoffSequence = 0;

        size_t numWindScenarios = scenario.windScenarios.size();
        if (numWindScenarios > 0 && scenario.pathList.size() > 0) {
          size_t numBasePaths = scenario.pathList.size() / numWindScenarios;
          size_t pathIdx = idx / numWindScenarios;
          size_t windIdx = idx % numWindScenarios;
          meta.pathVariantIndex = static_cast<int>(pathIdx);
          if (windIdx < scenario.windScenarios.size()) {
            meta.windVariantIndex = scenario.windScenarios[windIdx].windVariantIndex;
            meta.windSeed = scenario.windScenarios[windIdx].windSeed;
          }
        } else {
          meta.pathVariantIndex = static_cast<int>(idx);
        }

        populateVariationOffsets(meta);
        applySpeedProfileToPath(evalData.pathList[idx], meta.windVariantIndex);
        evalData.scenarioList.push_back(meta);
      }
      if (!evalData.scenarioList.empty()) {
        evalData.scenario = evalData.scenarioList.front();
      }
      evalData.sanitizePaths();

      // Send to minisim worker via ThreadPool
      auto evalDataPtr = std::make_shared<EvalData>(std::move(evalData));
      threadPool->enqueue([&genome, evalDataPtr](WorkerContext& context) {
        sendRPC(*context.socket, *evalDataPtr);
        context.evalResults = receiveRPC<EvalResults>(*context.socket);
        globalSimRunCounter.fetch_add(context.evalResults.pathList.size(), std::memory_order_relaxed);

        // Compute fitness
        genome.fitness = computeNNFitness(context.evalResults);
      });
    }

    // Wait for all evaluations to complete
    threadPool->wait_for_tasks();

    // Find best individual
    int bestIdx = 0;
    for (int i = 1; i < popSize; i++) {
      if (pop.individuals[i].fitness < pop.individuals[bestIdx].fitness) {
        bestIdx = i;
      }
    }
    pop.best_fitness = pop.individuals[bestIdx].fitness;
    pop.best_index = bestIdx;

    // Compute population fitness stats
    double sumFitness = 0, minFitness = pop.individuals[0].fitness, maxFitness = pop.individuals[0].fitness;
    for (const auto& ind : pop.individuals) {
      sumFitness += ind.fitness;
      if (ind.fitness < minFitness) minFitness = ind.fitness;
      if (ind.fitness > maxFitness) maxFitness = ind.fitness;
    }
    double avgFitness = sumFitness / popSize;

    // Re-evaluate best individual and log per-step data to data.dat
    {
      NNGenome& bestGenome = pop.individuals[bestIdx];
      std::vector<uint8_t> nnData;
      nn_serialize(bestGenome, nnData);

      const ScenarioDescriptor& scenario = scenarioForIndex(bestIdx % generationScenarios.size());
      uint64_t scenarioSequence = globalScenarioCounter.fetch_add(1, std::memory_order_relaxed) + 1;

      EvalData evalData;
      evalData.controllerType = ControllerType::NEURAL_NET;
      evalData.gp.assign(reinterpret_cast<const char*>(nnData.data()),
                         reinterpret_cast<const char*>(nnData.data() + nnData.size()));
      evalData.gpHash = hashByteVector(evalData.gp);
      evalData.isEliteReeval = false;
      evalData.pathList = scenario.pathList;
      evalData.scenario.scenarioSequence = scenarioSequence;
      evalData.scenario.bakeoffSequence = 0;

      evalData.scenarioList.clear();
      evalData.scenarioList.reserve(scenario.pathList.size());
      for (size_t idx = 0; idx < scenario.pathList.size(); ++idx) {
        ScenarioMetadata meta;
        meta.scenarioSequence = scenarioSequence;
        meta.enableDeterministicLogging = false;
        meta.bakeoffSequence = 0;
        size_t numWindScenarios = scenario.windScenarios.size();
        if (numWindScenarios > 0 && scenario.pathList.size() > 0) {
          size_t pathIdx = idx / numWindScenarios;
          size_t windIdx = idx % numWindScenarios;
          meta.pathVariantIndex = static_cast<int>(pathIdx);
          if (windIdx < scenario.windScenarios.size()) {
            meta.windVariantIndex = scenario.windScenarios[windIdx].windVariantIndex;
            meta.windSeed = scenario.windScenarios[windIdx].windSeed;
          }
        } else {
          meta.pathVariantIndex = static_cast<int>(idx);
        }
        populateVariationOffsets(meta);
        applySpeedProfileToPath(evalData.pathList[idx], meta.windVariantIndex);
        evalData.scenarioList.push_back(meta);
      }
      if (!evalData.scenarioList.empty()) {
        evalData.scenario = evalData.scenarioList.front();
      }
      evalData.sanitizePaths();

      auto evalDataPtr = std::make_shared<EvalData>(std::move(evalData));
      EvalResults bestResults;
      threadPool->enqueue([evalDataPtr, &bestResults](WorkerContext& context) {
        sendRPC(*context.socket, *evalDataPtr);
        bestResults = receiveRPC<EvalResults>(*context.socket);
      });
      threadPool->wait_for_tasks();

      logEvalResults(fout, bestResults);

      // Determinism check: re-eval fitness must match stored fitness exactly (T106)
      double reevalFitness = computeNNFitness(bestResults);
      double storedFitness = pop.individuals[bestIdx].fitness;
      if (!bitwiseEqual(reevalFitness, storedFitness)) {
        *logger.warn() << "NN_ELITE_DIVERGED: gen=" << gen
                       << " stored=" << std::fixed << std::setprecision(6) << storedFitness
                       << " reeval=" << reevalFitness
                       << " delta=" << std::scientific << (reevalFitness - storedFitness)
                       << endl;
      } else {
        *logger.info() << "NN_ELITE_SAME: gen=" << gen
                       << " fitness=" << std::fixed << std::setprecision(6) << storedFitness
                       << endl;
      }

      // Embed serialized NN genome in bestResults.gp for S3 storage
      // (renderer uses .gp field to detect format and extract fitness)
      bestResults.gp.assign(reinterpret_cast<const char*>(nnData.data()),
                            reinterpret_cast<const char*>(nnData.data() + nnData.size()));
      bestResults.gpHash = hashByteVector(bestResults.gp);

      // Save to S3 as Boost-serialized EvalResults (same format as GP)
      {
        std::string keyName = startTime + "/gen" + std::to_string(10000 - gen) + ".dmp";
        auto s3Client = ConfigManager::getS3Client();
        if (s3Client) {
          std::ostringstream oss(std::ios::binary);
          boost::archive::binary_oarchive oa(oss);
          oa << bestResults;

          auto stream = Aws::MakeShared<Aws::StringStream>("PutObject");
          *stream << oss.str();

          Aws::S3::Model::PutObjectRequest request;
          request.SetBucket(cfg.s3Bucket);
          request.SetKey(keyName);
          request.SetBody(stream);

          auto outcome = s3Client->PutObject(request);
          if (!outcome.IsSuccess()) {
            *logger.warn() << "S3 upload failed: " << outcome.GetError().GetMessage() << endl;
          }
        }
      }
    }

    // Log generation report
    *logger.info() << "Gen " << gen
                   << "  Best=" << std::fixed << std::setprecision(2) << minFitness
                   << "  Avg=" << avgFitness
                   << "  Worst=" << maxFitness
                   << "  Sigma=" << pop.individuals[bestIdx].mutation_sigma
                   << endl;

    // Log to statistics file
    bout << "#NNGen gen=" << gen
         << " best=" << std::fixed << std::setprecision(6) << minFitness
         << " avg=" << avgFitness
         << " worst=" << maxFitness
         << " bestSigma=" << pop.individuals[bestIdx].mutation_sigma
         << std::endl;
    bout.flush();

    logGenerationStats(gen);

    // Evolve next generation (skip on last gen)
    if (gen < numGens) {
      NNEvolveParams evoParams;
      evoParams.tournament_size = cfg.tournamentSize;          // GP & NN
      evoParams.crossover_prob = cfg.crossoverProbability;     // GP & NN
      evoParams.creation_prob = cfg.creationProbability;       // GP & NN
      evoParams.mutation_prob = cfg.swapMutationProbability;   // GP & NN
      evoParams.crossover_alpha = static_cast<float>(cfg.nnCrossoverAlpha); // NN only
      evoParams.elitism_count = cfg.addBestToNewPopulation;    // GP & NN
      nn_evolve_generation(pop, evoParams);
    }
  }

  *logger.info() << "NN evolution complete! Best fitness: "
                 << std::fixed << std::setprecision(2) << pop.best_fitness << endl;
}

int main(int argc, char** argv)
{
  // Parse command line arguments
  std::string configFile = "autoc.ini";
  
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-i") == 0) {
      if (i + 1 < argc) {
        configFile = argv[i + 1];
        i++; // Skip the next argument since we consumed it
      } else {
        std::cerr << "Error: -i option requires a filename argument" << std::endl;
        return 1;
      }
    } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      std::cout << "Usage: " << argv[0] << " [-i config_file]" << std::endl;
      std::cout << "  -i config_file  Use specified config file instead of autoc.ini" << std::endl;
      std::cout << "  -h, --help      Show this help message" << std::endl;
      return 0;
    }
  }

  logging::add_console_log(
    std::cout,
    boost::log::keywords::format = (
      boost::log::expressions::stream
      << boost::log::expressions::format_date_time<boost::posix_time::ptime>("TimeStamp", "%Y-%m-%d %H:%M:%S")
      << ": <" << logging::trivial::severity
      << "> " << boost::log::expressions::smessage
      )
  );
  logging::core::get()->set_filter(
    logging::trivial::severity >= logging::trivial::info
  );

  logging::add_common_attributes();

  logger = Logger();

  // Set up a new-handler, because we might need a lot of memory, and
  // we don't know it's there.
  set_new_handler(newHandler);

  // Initialize ConfigManager first so we can access gpSeed
  ConfigManager::initialize(configFile, *logger.info());
  const AutocConfig& cfg = ConfigManager::getConfig();

  // Init GP system with seed from config
  // Handle -1 as time-based seed (GP library doesn't recognize this convention)
  long gpSeed;
  if (cfg.gpSeed == -1) {
    gpSeed = static_cast<long>(time(NULL));
    *logger.info() << "GPSeed: -1 (auto) -> " << gpSeed << endl;
  } else {
    gpSeed = static_cast<long>(cfg.gpSeed);
  }
  rng::seed(static_cast<uint64_t>(gpSeed));

  // AWS setup
  Aws::SDKOptions options;
  Aws::InitAPI(options);

  // initialize workers
  threadPool = new ThreadPool(ConfigManager::getConfig());

  // Print the configuration
  *logger.info() << "Config: pop=" << cfg.populationSize
                 << " gens=" << cfg.numberOfGenerations << endl;
  *logger.info() << "SimNumPathsPerGen: " << cfg.simNumPathsPerGen << endl;
  *logger.info() << "EvalThreads: " << cfg.evalThreads << endl;
  *logger.info() << "MinisimProgram: " << cfg.minisimProgram << endl;
  *logger.info() << "MinisimPortOverride: " << cfg.minisimPortOverride << endl;
  *logger.info() << "EvaluateMode: " << cfg.evaluateMode << endl;
  *logger.info() << "ControllerType: " << cfg.controllerType << endl;
  *logger.info() << "WindScenarios: " << cfg.windScenarioCount << endl;
  *logger.info() << "GPSeed: " << cfg.gpSeed << endl;
  *logger.info() << "RandomPathSeedB: " << cfg.randomPathSeedB << endl;

  // Log VARIATIONS1 settings and initialize global sigmas
  *logger.info() << "EnableEntryVariations: " << cfg.enableEntryVariations << endl;
  *logger.info() << "EnableWindVariations: " << cfg.enableWindVariations << endl;

  // Initialize global variation parameters from config (degrees -> radians)
  // Store individual flags for selective application in populateVariationOffsets()
  gEnableEntryVariations = (cfg.enableEntryVariations != 0);
  gEnableWindVariations = (cfg.enableWindVariations != 0);
  gVariationSigmas = VariationSigmas::fromDegrees(
      cfg.entryHeadingSigma,
      cfg.entryRollSigma,
      cfg.entryPitchSigma,
      cfg.entrySpeedSigma,  // already a fraction
      cfg.windDirectionSigma,
      cfg.entryPositionRadiusSigma,  // meters (no conversion needed)
      cfg.entryPositionAltSigma      // meters (no conversion needed)
  );

  // Initialize global rabbit speed config
  gRabbitSpeedConfig = RabbitSpeedConfig{
      cfg.rabbitSpeedNominal,
      cfg.rabbitSpeedSigma,
      cfg.rabbitSpeedMin,
      cfg.rabbitSpeedMax,
      cfg.rabbitSpeedCycleMin,
      cfg.rabbitSpeedCycleMax
  };

  // Log rabbit speed configuration
  *logger.info() << "RabbitSpeed: nominal=" << cfg.rabbitSpeedNominal << " m/s"
                 << " sigma=" << cfg.rabbitSpeedSigma << " m/s"
                 << " range=[" << cfg.rabbitSpeedMin << ", " << cfg.rabbitSpeedMax << "] m/s"
                 << " cycles=[" << cfg.rabbitSpeedCycleMin << ", " << cfg.rabbitSpeedCycleMax << "] s"
                 << (cfg.rabbitSpeedSigma > 0 ? " (VARIABLE)" : " (CONSTANT)") << endl;

  // RAMP_LANDSCAPE: Initialize variation ramp globals
  gTotalGenerations = cfg.numberOfGenerations;
  gVariationRampStep = cfg.variationRampStep;
  gCurrentGeneration = 0;
  if (gVariationRampStep > 0) {
    int totalSteps = gTotalGenerations / gVariationRampStep;
    *logger.info() << "VariationRamp: step=" << gVariationRampStep
                   << " gens, " << totalSteps << " steps over " << gTotalGenerations
                   << " gens (scale +=" << std::fixed << std::setprecision(1)
                   << (100.0 / totalSteps) << "% per step)" << endl;
  }

  // SINGLE PRNG ARCHITECTURE: Pre-fetch all scenario variations from GPrand()
  // This consumes from GPrand() BEFORE GP evolution, ensuring deterministic sequence
  // When a variation type is disabled, defaults are stored in the table (not filtered later)
  int windScenarioCount = std::max(cfg.windScenarioCount, 1);
  prefetchAllVariations(windScenarioCount, gVariationSigmas, gRabbitSpeedConfig,
                        cfg.randomPathSeedB,
                        gEnableEntryVariations, gEnableWindVariations);

  // Log pre-fetched variations for verification
  *logger.info() << "Sigmas: heading=" << cfg.entryHeadingSigma << "° "
                 << "roll=" << cfg.entryRollSigma << "° "
                 << "pitch=" << cfg.entryPitchSigma << "° "
                 << "speed=" << (cfg.entrySpeedSigma * 100) << "% "
                 << "wind=" << cfg.windDirectionSigma << "°" << endl;
  logPrefetchedVariations(windScenarioCount, gpSeed);

  // Set training node mask before creating node set
  setTrainingNodesMask(cfg.trainingNodes.c_str());

  // Create the adf function/terminal set and print it out.
  createNodeSet(adfNs);
  *logger.info() << adfNs << endl;

  // Log training nodes (the subset used for tree creation/mutation)
  extern std::vector<std::string> trainingNodeNames;
  std::ostringstream nodeListStr;
  nodeListStr << "TrainingNodes (" << trainingNodeNames.size() << "): ";
  for (size_t i = 0; i < trainingNodeNames.size(); i++) {
    if (i > 0) nodeListStr << ",";
    nodeListStr << trainingNodeNames[i];
  }
  *logger.info() << nodeListStr.str() << endl;

  // Open the main output file for the data and statistics file.
  // First set up names for data file.  Remember we should delete the
  // string from the stream, well just a few bytes
  ostringstream strOutFile, strStatFile;
  const char* filePrefix = cfg.evaluateMode ? "eval-" : "";
  strOutFile << filePrefix << "data.dat" << ends;
  strStatFile << filePrefix << "data.stc" << ends;
  fout.open(strOutFile.str());
  bout.open(strStatFile.str());

  // Set fixed-point notation for statistics file
  bout << std::fixed << std::setprecision(6);

  std::string startTime = generate_iso8601_timestamp();
  auto runStartTime = std::chrono::steady_clock::now();
  auto lastThroughputTime = runStartTime;
  uint64_t lastSimRunCount = 0;
  auto logGenerationStats = [&](int genIndex) {
    auto now = std::chrono::steady_clock::now();
    gp_scalar deltaSec = std::chrono::duration<gp_scalar>(now - lastThroughputTime).count();
    uint64_t currentRuns = globalSimRunCounter.load(std::memory_order_relaxed);
    uint64_t deltaRuns = currentRuns - lastSimRunCount;
    gp_scalar rate = (deltaSec > static_cast<gp_scalar>(0.0f))
      ? static_cast<gp_scalar>(deltaRuns) / deltaSec
      : static_cast<gp_scalar>(0.0f);
    bout << std::fixed << std::setprecision(2)
         << "#GenSimStats gen=" << genIndex
         << " sims=" << deltaRuns
         << " total=" << currentRuns
         << " durationSec=" << deltaSec
         << " rate=" << rate
         << std::setprecision(6) << std::endl;
    bout.flush();
    lastThroughputTime = now;
    lastSimRunCount = currentRuns;
  };

  // Generate initial paths using pre-fetched gPathSeed (single PRNG architecture)
  generationPaths = generateSmoothPaths(const_cast<char*>(cfg.generatorMethod.c_str()),
                                        cfg.simNumPathsPerGen,
                                        SIM_PATH_BOUNDS, SIM_PATH_BOUNDS,
                                        gPathSeed);

  // DEBUG: Log segment counts for each path
  std::cout << "\n=== Path Segment Counts (method=" << cfg.generatorMethod
            << ", seed=" << gPathSeed << ") ===" << std::endl;
  int maxSegments = 0;
  for (size_t i = 0; i < generationPaths.size(); i++) {
    int segCount = generationPaths[i].size();
    std::cout << "  Path " << i << ": " << segCount << " segments" << std::endl;
    if (segCount > maxSegments) maxSegments = segCount;
  }
  std::cout << "  Maximum: " << maxSegments << " segments" << std::endl;
  std::cout << "===" << std::endl << std::endl;

  rebuildGenerationScenarios(generationPaths);
  const int windsPerPath = std::max(cfg.windScenarioCount, 1);
  *logger.debug() << "Wind scenarios this generation: paths="
                 << cfg.simNumPathsPerGen
                 << " windsPerPath=" << windsPerPath
                 << " dispatchScenarios=" << generationScenarios.size()
                 << " totalEvaluations=" << generationScenarios.size() * windsPerPath
                 << endl;
  warnIfScenarioMismatch();

  // Check controller type: NN evolution gets its own dedicated loop
  bool isNNMode = (cfg.controllerType == "NN");

  if (isNNMode && cfg.evaluateMode) {
    // NN evaluation mode: load weight file, evaluate, report fitness
    runNNEvaluation(startTime, fout, bout);
  } else if (isNNMode) {
    // NN evolution mode
    runNNEvolution(startTime, runStartTime, fout, bout, logGenerationStats);
  } else if (cfg.evaluateMode) {
    // Evaluation mode: use bytecode interpreter with external simulators
    *logger.info() << "Running in bytecode evaluation mode with external simulators..." << endl;
    *logger.info() << "Loading bytecode file: " << cfg.bytecodeFile << endl;
    
    // Load the bytecode program
    GPBytecodeInterpreter interpreter;
    if (!interpreter.loadProgram(cfg.bytecodeFile)) {
      *logger.error() << "Failed to load bytecode file: " << cfg.bytecodeFile << endl;
      exit(1);
    }
    
    // Create a custom evaluation GP class for bytecode
    class BytecodeEvaluationGP : public MyGP {
    public:
      std::vector<char> bytecodeBuffer;
      
      BytecodeEvaluationGP(const std::vector<char>& buffer) : MyGP(1), bytecodeBuffer(buffer) {}
      
      virtual void evaluate() override {
        tasks.push_back(this);
      }
      
      virtual void evalTask(WorkerContext& context) override {
        stdFitness = 0;
        
        // Send bytecode evaluation request to simulator
        const ScenarioDescriptor& scenario = scenarioForIndex(getScenarioIndex());
        EvalData evalData;
        evalData.gp = bytecodeBuffer;
        evalData.gpHash = hashByteVector(evalData.gp);
        evalData.pathList = scenario.pathList;
        uint64_t scenarioSequence = globalScenarioCounter.fetch_add(1, std::memory_order_relaxed) + 1;
        evalData.scenario.scenarioSequence = scenarioSequence;
        evalData.scenarioList.clear();
        evalData.scenarioList.reserve(scenario.pathList.size());
        bool isBakeoff = isBakeoffMode();
        const AutocConfig& cfg = ConfigManager::getConfig();

        for (size_t idx = 0; idx < scenario.pathList.size(); ++idx) {
          ScenarioMetadata meta;
          meta.scenarioSequence = scenarioSequence;
          if (isBakeoff) {
            meta.bakeoffSequence = bakeoffPathCounter.fetch_add(1, std::memory_order_relaxed) + 1;
          } else {
            meta.bakeoffSequence = 0;
          }

          if (cfg.demeticGrouping && cfg.demeSize > 0) {
            // Demetic mode: scenario has one path variant across multiple winds
            meta.pathVariantIndex = scenario.pathVariantIndex;
            if (idx < scenario.windScenarios.size()) {
              meta.windVariantIndex = scenario.windScenarios[idx].windVariantIndex;
              meta.windSeed = scenario.windScenarios[idx].windSeed;
            } else {
              meta.windVariantIndex = scenario.windVariantIndex;
              meta.windSeed = scenario.windSeed;
            }
          } else {
            // Non-demetic mode: scenario has all paths × all winds
            // pathList is organized in PATH-MAJOR order: [p0w0, p0w1, ..., p0w5, p1w0, p1w1, ...]
            size_t numWindScenarios = scenario.windScenarios.size();
            if (numWindScenarios > 0 && scenario.pathList.size() > 0) {
              size_t numBasePaths = scenario.pathList.size() / numWindScenarios;
              size_t pathIdx = idx / numWindScenarios;
              size_t windIdx = idx % numWindScenarios;

              meta.pathVariantIndex = static_cast<int>(pathIdx);
              if (windIdx < scenario.windScenarios.size()) {
                meta.windVariantIndex = scenario.windScenarios[windIdx].windVariantIndex;
                meta.windSeed = scenario.windScenarios[windIdx].windSeed;
              } else {
                meta.windVariantIndex = scenario.windVariantIndex;
                meta.windSeed = scenario.windSeed;
              }
            } else {
              // Fallback if structure is unexpected
              meta.pathVariantIndex = static_cast<int>(idx);
              meta.windVariantIndex = scenario.windVariantIndex;
              meta.windSeed = scenario.windSeed;
            }
          }
          // VARIATIONS1: Populate entry/wind variation offsets (uses pre-fetched values)
          populateVariationOffsets(meta);

          // Variable rabbit speed: Apply speed profile to path timing (uses pre-fetched profile)
          applySpeedProfileToPath(evalData.pathList[idx], meta.windVariantIndex);

          evalData.scenarioList.push_back(meta);
        }
        if (!evalData.scenarioList.empty()) {
          evalData.scenario = evalData.scenarioList.front();
        } else {
          evalData.scenario.scenarioSequence = scenarioSequence;
          evalData.scenario.bakeoffSequence = isBakeoff
            ? bakeoffPathCounter.fetch_add(1, std::memory_order_relaxed) + 1
            : 0;
        }
        evalData.sanitizePaths();
        const bool logDispatch = enableDeterministicTestLogging.load(std::memory_order_relaxed);
        if (logDispatch) {
          *logger.info() << "[AUTOC_SEND_GP] workerId=" << context.workerId
                         << " gpHash=0x" << std::hex << evalData.gpHash << std::dec
                         << " bytes=" << evalData.gp.size()
                         << " scenarioSeq=" << scenarioSequence
                         << " (bytecode eval)"
                         << endl;
        }
        sendRPC(*context.socket, evalData);
        
        // Get simulation results
        context.evalResults = receiveRPC<EvalResults>(*context.socket);

        if (context.evalResults.gpHash == 0 && !context.evalResults.gp.empty()) {
          context.evalResults.gpHash = hashByteVector(context.evalResults.gp);
        }

        if (context.evalResults.pathList.empty() && !evalData.pathList.empty()) {
          context.evalResults.pathList = evalData.pathList;
        }
        if (context.evalResults.scenario.windSeed == 0 && evalData.scenario.windSeed != 0) {
          context.evalResults.scenario = evalData.scenario;
        }
        if (context.evalResults.scenarioList.empty() && !evalData.scenarioList.empty()) {
          context.evalResults.scenarioList = evalData.scenarioList;
        }
        globalSimRunCounter.fetch_add(context.evalResults.pathList.size(), std::memory_order_relaxed);

        if (printEval) {
          boost::unique_lock<boost::mutex> lock(evalCollectorMutex);
          if (activeEvalCollector) {
            appendEvalResults(*activeEvalCollector, context.evalResults);
          } else {
            bestOfEvalResults = context.evalResults;
          }
        }
        
        // Compute fitness for each path and sum
        // Crash penalty (soft lexicographic) applied per-path
        for (int i = 0; i < context.evalResults.pathList.size(); i++) {
          auto& path = context.evalResults.pathList.at(i);
          auto& aircraftState = context.evalResults.aircraftStateList.at(i);
          auto& crashReason = context.evalResults.crashReasonList.at(i);

          if (path.empty() || aircraftState.empty()) {
            continue;
          }

          // Get scenario metadata for this path (for logging)
          int pathVariantIndex = i;  // default to loop index
          int windVariantIndex = -1;
          unsigned int windSeed = 0;
          if (i < context.evalResults.scenarioList.size()) {
            pathVariantIndex = context.evalResults.scenarioList.at(i).pathVariantIndex;
            windVariantIndex = context.evalResults.scenarioList.at(i).windVariantIndex;
            windSeed = context.evalResults.scenarioList.at(i).windSeed;
          }

          // ========================================================================
          // SIMPLIFIED FITNESS (v3): Two objectives - distance + attitude delta
          // See specs/FITNESS_SIMPLIFY_20260221.md
          // ========================================================================
          gp_fitness localFitness = 0.0;
          int stepIndex = 0;

          // Accumulate raw values, scale attitude at end using path geometry
          gp_fitness distance_sum = 0.0;
          gp_fitness attitude_sum = 0.0;
          int simulation_steps = 0;

          // Intercept-budget scaling: compute budget from initial displacement + heading
          double interceptBudget = 0.0;
          {
            gp_vec3 initialPos = aircraftState.at(0).getPosition();
            gp_vec3 pathStart = path.at(0).start;
            gp_vec3 offset3d = initialPos - pathStart;
            double displacement = sqrt(offset3d[0] * offset3d[0] + offset3d[1] * offset3d[1]);
            double headingOffset = 0.0;
            if (i < context.evalResults.scenarioList.size()) {
              headingOffset = context.evalResults.scenarioList.at(i).entryHeadingOffset;
            }
            interceptBudget = computeInterceptBudget(displacement, headingOffset,
                                                      SIM_INITIAL_VELOCITY, gRabbitSpeedConfig.nominal);
          }

          // Attitude delta tracking
          gp_fitness prev_roll = 0.0, prev_pitch = 0.0;
          bool first_attitude_sample = true;

          while (++stepIndex < aircraftState.size()) {
            auto& stepAircraftState = aircraftState.at(stepIndex);
            int rawPathIndex = stepAircraftState.getThisPathIndex();
            int pathIndex = std::clamp(rawPathIndex, 0, static_cast<int>(path.size()) - 1);
            const Path& currentPathPoint = path.at(pathIndex);

            gp_vec3 aircraftPosition = stepAircraftState.getPosition();
            gp_quat craftOrientation = stepAircraftState.getOrientation();

            // ====================================================================
            // ATTITUDE DELTA: Per-step roll/pitch change (radians)
            // ====================================================================
            gp_fitness qw = craftOrientation.w();
            gp_fitness qx = craftOrientation.x();
            gp_fitness qy = craftOrientation.y();
            gp_fitness qz = craftOrientation.z();

            // Roll (phi): rotation around X axis
            gp_fitness roll = atan2(2.0 * (qw * qx + qy * qz),
                                    1.0 - 2.0 * (qx * qx + qy * qy));
            // Pitch (theta): rotation around Y axis
            gp_fitness sinp = 2.0 * (qw * qy - qz * qx);
            sinp = std::clamp(sinp, -1.0, 1.0);
            gp_fitness pitch = asin(sinp);

            gp_fitness attitude_delta = 0.0;
            if (first_attitude_sample) {
              prev_roll = roll;
              prev_pitch = pitch;
              first_attitude_sample = false;
            } else {
              attitude_delta = fabs(roll - prev_roll) + fabs(pitch - prev_pitch);
              prev_roll = roll;
              prev_pitch = pitch;
            }

            // ====================================================================
            // DISTANCE: 3D distance to rabbit (meters)
            // ====================================================================
            gp_vec3 craftOffset = aircraftPosition - currentPathPoint.start;
            gp_scalar distance = craftOffset.norm();

            // Accumulate with normalized nonlinear penalty per step
            // Intercept-budget scaling: ramp penalty from floor→ceiling over budget
            double stepTimeSec = stepAircraftState.getSimTimeMsec() / 1000.0;
            double interceptScale = computeInterceptScale(stepTimeSec, interceptBudget);
            double distDelta = fabs(distance - DISTANCE_TARGET);
            distance_sum += pow(interceptScale * distDelta / DISTANCE_NORM, DISTANCE_POWER);
            attitude_sum += pow(interceptScale * attitude_delta / ATTITUDE_NORM, ATTITUDE_POWER);
            simulation_steps++;

            // Per-step logging
            if (printEval) {
              bestOfEvalResults = context.evalResults;

              if (stepIndex == 1) {
                fout << "Pth:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power      qw      qx      qy      qz   vxBody   vyBody   vzBody    alpha     beta   dtheta    dphi   dhome     dist   attDlt  rabVel\n";
              }

              // Compute rabbit velocity from path position delta
              static gp_vec3 prevPathPos2(0, 0, 0);
              static long prevPathTime2 = 0;
              static bool firstStep2 = true;
              gp_scalar rabbitVel = 0.0f;
              gp_vec3 currPathPos = path.at(pathIndex).start;
              long currPathTime = path.at(pathIndex).simTimeMsec;
              if (stepIndex == 1) {
                firstStep2 = true;
              }
              if (!firstStep2 && currPathTime > prevPathTime2) {
                gp_scalar dist = (currPathPos - prevPathPos2).norm();
                gp_scalar dt = static_cast<gp_scalar>(currPathTime - prevPathTime2) / 1000.0f;
                rabbitVel = dist / dt;
              }
              prevPathPos2 = currPathPos;
              prevPathTime2 = currPathTime;
              firstStep2 = false;

              // Convert orientation to euler angles
              Eigen::Matrix3f rotMatrix = craftOrientation.toRotationMatrix();
              gp_vec3 euler;
              if (std::abs(rotMatrix(2, 0)) > static_cast<gp_scalar>(0.99999f)) {
                euler[0] = 0;
                if (rotMatrix(2, 0) > 0) {
                  euler[1] = -static_cast<gp_scalar>(M_PI / 2.0);
                  euler[2] = -atan2(rotMatrix(1, 2), rotMatrix(0, 2));
                } else {
                  euler[1] = static_cast<gp_scalar>(M_PI / 2.0);
                  euler[2] = atan2(rotMatrix(1, 2), rotMatrix(0, 2));
                }
              } else {
                euler[0] = atan2(rotMatrix(2, 1), rotMatrix(2, 2));
                euler[1] = -asin(rotMatrix(2, 0));
                euler[2] = atan2(rotMatrix(1, 0), rotMatrix(0, 0));
              }

              // Body-frame velocity
              gp_vec3 velocity_body = stepAircraftState.getOrientation().inverse() *
                                      stepAircraftState.getVelocity();

              // GP operator values
              gp_scalar alpha_deg = atan2(-velocity_body.z(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);
              gp_scalar beta_deg = atan2(velocity_body.y(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);

              // Angle to target in body frame
              gp_vec3 craftToTarget = path.at(pathIndex).start - stepAircraftState.getPosition();
              gp_vec3 target_body = stepAircraftState.getOrientation().inverse() * craftToTarget;
              gp_scalar dtheta_deg = atan2(-target_body.z(), target_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);
              gp_scalar dphi_deg = atan2(target_body.y(), -target_body.z()) * static_cast<gp_scalar>(180.0 / M_PI);

              // Distance to home
              gp_vec3 home(0, 0, SIM_INITIAL_ALTITUDE);
              gp_scalar dhome = (home - stepAircraftState.getPosition()).norm();

              gp_quat q = craftOrientation;

              char outbuf[1500];
              sprintf(outbuf, "%03d:%04d: %06ld %3d % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 7.4f % 7.4f % 7.4f % 7.4f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.3f % 6.4f % 6.1f\n",
                pathVariantIndex, simulation_steps,
                stepAircraftState.getSimTimeMsec(), pathIndex,
                path.at(pathIndex).distanceFromStart,
                path.at(pathIndex).start[0],
                path.at(pathIndex).start[1],
                path.at(pathIndex).start[2],
                stepAircraftState.getPosition()[0],
                stepAircraftState.getPosition()[1],
                stepAircraftState.getPosition()[2],
                euler[2], euler[1], euler[0],
                stepAircraftState.getRelVel(),
                stepAircraftState.getRollCommand(),
                stepAircraftState.getPitchCommand(),
                stepAircraftState.getThrottleCommand(),
                q.w(), q.x(), q.y(), q.z(),
                velocity_body.x(), velocity_body.y(), velocity_body.z(),
                alpha_deg, beta_deg,
                dtheta_deg, dphi_deg,
                dhome,
                distance,        // New: distance to rabbit
                attitude_delta,  // New: attitude change this step
                rabbitVel
              );
              fout << outbuf;
            }
          }

          // Scale attitude by path geometry so 1 rad excess ≈ (path_dist/path_turn) meters
          gp_scalar path_distance = path.back().distanceFromStart;
          gp_scalar path_turn_rad = path.back().radiansFromStart;
          // Fallback for straight paths: treat as one full rotation worth
          gp_scalar attitude_scale = path_distance / std::max(path_turn_rad, static_cast<gp_scalar>(2.0 * M_PI));

          localFitness = distance_sum + attitude_sum * attitude_scale;

          if (isnan(localFitness)) {
            nanDetector++;
          }

          // Crash penalty: soft lexicographic - completion dominates, quality provides gradient
          if (crashReason != CrashReason::None) {
            gp_fitness total_path_distance = path.back().distanceFromStart;
            gp_fitness fraction_completed =
              static_cast<gp_fitness>(path.at(aircraftState.back().getThisPathIndex()).distanceFromStart) / total_path_distance;
            fraction_completed = std::max(fraction_completed, static_cast<gp_fitness>(0.001));
            gp_fitness completion_penalty = (1.0 - fraction_completed) * CRASH_COMPLETION_WEIGHT;
            localFitness = completion_penalty + localFitness;
          }

          // Sum path fitness
          stdFitness += localFitness;
        }
        
        // Mark fitness as valid for the GP library
        fitnessValid = 1;
      }
      
      gp_fitness getFinalFitness() const { return stdFitness; }
    };
    
    // Serialize the bytecode interpreter for transport to simulators
    std::vector<char> bytecodeBuffer;
    {
      boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> outStream(bytecodeBuffer);
      boost::archive::binary_oarchive archive(outStream);
      archive << interpreter;
      outStream.flush();
    }
    
    *logger.info() << "Running bytecode evaluation on " << generationScenarios.size() << " scenarios..." << endl;
    
    // Initialize evaluation
    printEval = true; // Enable detailed output
    
    // Create bytecode evaluation GP
    BytecodeEvaluationGP evalGP(bytecodeBuffer);
    
    aggregatedEvalResults.clear();
    EvalResults* previousCollector = activeEvalCollector;
    activeEvalCollector = &aggregatedEvalResults;
    
    evaluationProgress.store(0);
    gp_fitness cumulativeFitness = 0.0;
    size_t scenariosEvaluated = 0;

    for (size_t scenarioIdx = 0; scenarioIdx < generationScenarios.size(); ++scenarioIdx) {
      const auto& scenario = scenarioForIndex(static_cast<int>(scenarioIdx));
      *logger.debug() << "Evaluating scenario " << scenarioIdx
                      << " paths=" << scenario.pathList.size() << endl;
      evalGP.setScenarioIndex(static_cast<int>(scenarioIdx));
      evalGP.evaluate();

      for (auto& task : tasks) {
        threadPool->enqueue([task](WorkerContext& context) {
          task->evalTask(context);
          });
      }

      threadPool->wait_for_tasks();
      tasks.clear();

      cumulativeFitness += evalGP.getFinalFitness();
      scenariosEvaluated++;
    }

    activeEvalCollector = previousCollector;
    printEval = false;

    gp_fitness averageFitness = scenariosEvaluated > 0 ? cumulativeFitness / static_cast<gp_fitness>(scenariosEvaluated) : 0.0;
    *logger.info() << "Aggregated results: paths=" << aggregatedEvalResults.pathList.size()
                   << " states=" << aggregatedEvalResults.aircraftStateList.size() << endl;
    bestOfEvalResults = aggregatedEvalResults;

    // Serialize a GP object with the computed fitness for renderer display
    MyGP fitnessHolder(1);
    fitnessHolder.setFitness(averageFitness);
    std::vector<char> gpBuffer;
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> gpOutStream(gpBuffer);
    fitnessHolder.save(gpOutStream);
    gpOutStream.flush();
    bestOfEvalResults.gp = gpBuffer;

    *logger.info() << "Bytecode evaluation complete!" << endl;
    *logger.info() << "Computed fitness: " << averageFitness << endl;
    *logger.info() << "Original GP fitness: " << interpreter.getFitness() << endl;
    
    // Set up S3 storage for results (use gen-1.dmp so renderer can interpret it)
    computedKeyName = startTime + "/gen-1.dmp";
    
    // Store results to S3
    if (bestOfEvalResults.pathList.size() > 0) {
      Aws::S3::Model::PutObjectRequest request;
      request.SetBucket(cfg.s3Bucket);
      request.SetKey(computedKeyName);

      std::ostringstream oss(std::ios::binary);
      boost::archive::binary_oarchive oa(oss);
      oa << bestOfEvalResults;

      std::shared_ptr<Aws::StringStream> ss = Aws::MakeShared<Aws::StringStream>("");
      *ss << oss.str();
      request.SetBody(ss);

      auto outcome = getS3Client()->PutObject(request);
      if (!outcome.IsSuccess()) {
        *logger.error() << "Error storing results to S3: " << outcome.GetError().GetMessage() << std::endl;
      } else {
        *logger.info() << "Results stored to S3: " << computedKeyName << std::endl;
      }
    }
  } else {
    // Normal GP evolution mode
    // Create a population with this configuration
    *logger.info() << "Creating initial population ..." << endl;
    // TODO(T017-T018): Remove GPVariables bridge once GP code is deleted
    MyPopulation* pop = new MyPopulation(ConfigManager::getConfig(), adfNs);
    pop->create();
    *logger.info() << "Ok." << endl;
    pop->createGenerationReport(1, 0, fout, bout, *logger.info());
    logGenerationStats(0);

    // This next for statement is the actual genetic programming system
    // which is in essence just repeated reproduction and crossover loop
    // through all the generations ...
    MyPopulation* newPop = NULL;

    for (int gen = 1; gen <= cfg.numberOfGenerations; gen++)
    {
      // RAMP_LANDSCAPE: Update current generation for variation scaling
      int prevGen = gCurrentGeneration;
      gCurrentGeneration = gen;

      // Log when scale changes (at step boundaries)
      if (gVariationRampStep > 0) {
        int prevStep = prevGen / gVariationRampStep;
        int currStep = gen / gVariationRampStep;
        if (currStep != prevStep || gen == 1) {
          double scale = computeVariationScale();
          int totalSteps = gTotalGenerations / gVariationRampStep;
          *logger.info() << "VariationScale: " << std::fixed << std::setprecision(2)
                         << (scale * 100.0) << "% (step " << currStep << "/" << totalSteps << ")" << endl;
        }
      }

      // For this generation, build a smooth path goal using pre-fetched gPathSeed
      // (single PRNG architecture: consistent paths derived from GPSeed)
      generationPaths = generateSmoothPaths(const_cast<char*>(cfg.generatorMethod.c_str()),
                                            cfg.simNumPathsPerGen,
                                            SIM_PATH_BOUNDS, SIM_PATH_BOUNDS,
                                            gPathSeed);
      rebuildGenerationScenarios(generationPaths);
      warnIfScenarioMismatch();

      // Create a new generation from the old one by applying the genetic operators
      if (!cfg.steadyState)
        // TODO(T017-T018): Remove GPVariables bridge once GP code is deleted
        newPop = new MyPopulation(ConfigManager::getConfig(), adfNs);
      pop->generate(*newPop);

      // Enable evaluation output for this generation
      printEval = true;

      // Switch to new population first to get the correct best individual
      if (!cfg.steadyState)
      {
        MyPopulation* oldPop = pop;
        pop = newPop;
        delete oldPop;
      }

      // reverse order names for s3...
      computedKeyName = startTime + "/gen" + std::to_string(10000 - gen) + ".dmp";
      pop->endOfEvaluation();

      printEval = false;

      // Create a report of this generation and how well it is doing
      if (nanDetector > 0) {
        *logger.warn() << "NanDetector count: " << nanDetector << endl;
      }
      auto logStream = logger.info();
      pop->createGenerationReport(0, gen, fout, bout, *logStream);
      logGenerationStats(gen);

      // ELITE TRACE CAPTURE: Save full simulation trace when elite is evaluated
      // Compare traces when fitness worsens to identify root cause

      // Helper to compute ULP (Units in Last Place) difference for float32
      auto computeULP = [](float a, float b) -> int64_t {
        if (std::isnan(a) || std::isnan(b)) return INT64_MAX;
        if (a == b) return 0;

        uint32_t bits_a, bits_b;
        std::memcpy(&bits_a, &a, sizeof(float));
        std::memcpy(&bits_b, &b, sizeof(float));

        // Handle sign bit - convert to signed magnitude representation
        int32_t signed_a = (bits_a & 0x80000000) ? (0x80000000 - (bits_a & 0x7FFFFFFF)) : bits_a;
        int32_t signed_b = (bits_b & 0x80000000) ? (0x80000000 - (bits_b & 0x7FFFFFFF)) : bits_b;

        return std::abs(static_cast<int64_t>(signed_a) - static_cast<int64_t>(signed_b));
      };

      // In-memory trace comparison function (only called when gEliteTraceEnabled)
      auto compareTraces = [&](const EvalResults& trace1, const EvalResults& trace2,
                               const std::string& label1, const std::string& label2) {
        *logger.warn() << "=== IN-MEMORY TRACE COMPARISON ===" << endl;
        *logger.warn() << "Trace 1: " << label1
                      << " workerId=" << trace1.workerId
                      << " pid=" << trace1.workerPid
                      << " evalCounter=" << trace1.workerEvalCounter << endl;
        *logger.warn() << "Trace 2: " << label2
                      << " workerId=" << trace2.workerId
                      << " pid=" << trace2.workerPid
                      << " evalCounter=" << trace2.workerEvalCounter << endl;
        *logger.warn() << endl;

        // Helper for computing ULP distance between doubles (native physics precision)
        auto computeULPDouble = [](double a, double b) -> int64_t {
          if (a == b) return 0;
          if (std::isnan(a) || std::isnan(b)) return INT64_MAX;

          int64_t ia, ib;
          memcpy(&ia, &a, sizeof(double));
          memcpy(&ib, &b, sizeof(double));

          // Handle sign bit
          if (ia < 0) ia = INT64_MIN - ia;
          if (ib < 0) ib = INT64_MIN - ib;

          return std::abs(ia - ib);
        };

        // Helper for double vectors (3D)
        auto ulpDouble3 = [&](const double* a, const double* b, int64_t out[3]) {
          for (int i = 0; i < 3; ++i) {
            out[i] = computeULPDouble(a[i], b[i]);
          }
        };

        // Require physicsTrace - no fallback to debugSamples
        if (trace1.physicsTrace.empty() || trace2.physicsTrace.empty()) {
          *logger.error() << "*** ERROR: Physics trace data missing! ***" << endl;
          *logger.error() << "  Trace1 physicsTrace size: " << trace1.physicsTrace.size()
                         << " (expected non-empty)" << endl;
          *logger.error() << "  Trace2 physicsTrace size: " << trace2.physicsTrace.size()
                         << " (expected non-empty)" << endl;
          *logger.error() << "  This indicates crrcsim binary was not recompiled with physics trace code." << endl;
          *logger.error() << "  Trace comparison ABORTED - cannot proceed without native precision data." << endl;
          return;
        }

        *logger.warn() << "=== PHYSICS TRACE COMPARISON (native double precision) ===" << endl;

          size_t pathCount = std::min(trace1.physicsTrace.size(), trace2.physicsTrace.size());

          // First pass: scan for RNG divergence to find where non-determinism starts
          int firstRngDivergenceStep = -1;
          size_t firstRngDivergencePath = 0;
          for (size_t pathIdx = 0; pathIdx < pathCount && firstRngDivergenceStep < 0; ++pathIdx) {
            const auto& p1 = trace1.physicsTrace[pathIdx];
            const auto& p2 = trace2.physicsTrace[pathIdx];
            if (p1.empty() || p2.empty()) continue;

            size_t stepCount = std::min(p1.size(), p2.size());
            for (size_t stepIdx = 0; stepIdx < stepCount; ++stepIdx) {
              const auto& s1 = p1[stepIdx];
              const auto& s2 = p2[stepIdx];
              bool rngDiff = (s1.rngState16 != s2.rngState16) || (s1.rngState32 != s2.rngState32);
              if (rngDiff) {
                firstRngDivergenceStep = static_cast<int>(stepIdx);
                firstRngDivergencePath = pathIdx;
                break;
              }
            }
          }

          if (firstRngDivergenceStep >= 0) {
            *logger.warn() << "*** RNG DIVERGENCE DETECTED ***" << endl;
            *logger.warn() << "First RNG divergence at Path " << firstRngDivergencePath
                          << ", Step " << firstRngDivergenceStep
                          << " (t=" << (firstRngDivergenceStep * 3) << "ms physics step)" << endl;
            *logger.warn() << endl;

            // Print context window: 3 steps before and 2 steps after the divergence
            const auto& p1 = trace1.physicsTrace[firstRngDivergencePath];
            const auto& p2 = trace2.physicsTrace[firstRngDivergencePath];
            size_t stepCount = std::min(p1.size(), p2.size());

            int contextBefore = 3;
            int contextAfter = 2;
            int startStep = std::max(0, firstRngDivergenceStep - contextBefore);
            int endStep = std::min(static_cast<int>(stepCount) - 1, firstRngDivergenceStep + contextAfter);

            *logger.warn() << "Context window (steps " << startStep << "-" << endStep << "):" << endl;
            for (int stepIdx = startStep; stepIdx <= endStep; ++stepIdx) {
              const auto& s1 = p1[stepIdx];
              const auto& s2 = p2[stepIdx];
              bool rngDiff = (s1.rngState16 != s2.rngState16) || (s1.rngState32 != s2.rngState32);

              std::string marker = (stepIdx == firstRngDivergenceStep) ? " <<< FIRST DIVERGENCE" : "";
              *logger.warn() << "  Step " << stepIdx << " (t=" << (stepIdx * 3) << "ms):" << marker << endl;
              *logger.warn() << "    RNG state16: " << s1.rngState16 << " vs " << s2.rngState16;
              if (rngDiff) *logger.warn() << " [DIFFER]";
              *logger.warn() << endl;
              *logger.warn() << "    RNG state32: " << s1.rngState32 << " vs " << s2.rngState32;
              if (rngDiff) *logger.warn() << " [DIFFER]";
              *logger.warn() << endl;

              // Show gust values for context
              *logger.warn() << "    Gust[xyz]: [" << s1.gustBody[0] << ", " << s1.gustBody[1] << ", " << s1.gustBody[2] << "]"
                            << " vs [" << s2.gustBody[0] << ", " << s2.gustBody[1] << ", " << s2.gustBody[2] << "]" << endl;
            }
            *logger.warn() << endl;

            // Now print full detail for the first divergence step only
            const auto& s1 = p1[firstRngDivergenceStep];
            const auto& s2 = p2[firstRngDivergenceStep];

            *logger.warn() << "=== DETAILED DIVERGENCE REPORT (Step " << firstRngDivergenceStep << ") ===" << endl;
          } else {
            *logger.warn() << "No RNG divergence found in physics trace (all " << pathCount << " paths checked)" << endl;
            *logger.warn() << endl;
          }

          // Detailed reporting for first divergence only
          for (size_t pathIdx = 0; pathIdx < pathCount; ++pathIdx) {
            const auto& p1 = trace1.physicsTrace[pathIdx];
            const auto& p2 = trace2.physicsTrace[pathIdx];

            if (p1.empty() || p2.empty()) continue;
            size_t stepCount = std::min(p1.size(), p2.size());

            for (size_t stepIdx = 0; stepIdx < stepCount; ++stepIdx) {
              const auto& s1 = p1[stepIdx];
              const auto& s2 = p2[stepIdx];

              // Check for ANY divergence in key physics values
              bool posDiff = (s1.pos[0] != s2.pos[0]) || (s1.pos[1] != s2.pos[1]) || (s1.pos[2] != s2.pos[2]);
              bool velDiff = (s1.vel[0] != s2.vel[0]) || (s1.vel[1] != s2.vel[1]) || (s1.vel[2] != s2.vel[2]);
              bool accDiff = (s1.acc[0] != s2.acc[0]) || (s1.acc[1] != s2.acc[1]) || (s1.acc[2] != s2.acc[2]);
              bool quatDiff = (s1.quat[0] != s2.quat[0]) || (s1.quat[1] != s2.quat[1]) ||
                             (s1.quat[2] != s2.quat[2]) || (s1.quat[3] != s2.quat[3]);
              bool alphaDiff = (s1.alpha != s2.alpha);
              bool betaDiff = (s1.beta != s2.beta);
              bool CLDiff = (s1.CL != s2.CL);
              bool CDDiff = (s1.CD != s2.CD);
              bool cmdDiff = (s1.pitchCommand != s2.pitchCommand) ||
                            (s1.rollCommand != s2.rollCommand) ||
                            (s1.throttleCommand != s2.throttleCommand);
              bool gustDiff = false;
              for (int i = 0; i < 6; ++i) {
                if (s1.gustBody[i] != s2.gustBody[i]) {
                  gustDiff = true;
                  break;
                }
              }
              bool rngDiff = (s1.rngState16 != s2.rngState16) || (s1.rngState32 != s2.rngState32);

              if (posDiff || velDiff || accDiff || quatDiff || alphaDiff || betaDiff ||
                  CLDiff || CDDiff || cmdDiff || gustDiff || rngDiff) {

                *logger.warn() << "*** PHYSICS TRACE DIVERGENCE ***" << endl;
                *logger.warn() << "  Path " << pathIdx << ", Step " << stepIdx
                              << " (step=" << s1.step << " vs " << s2.step << ")" << endl;
                *logger.warn() << "  Time: " << s1.simTimeMsec << " vs " << s2.simTimeMsec << " ms" << endl;
                *logger.warn() << "  Worker1: id=" << s1.workerId << " pid=" << s1.workerPid
                              << " eval=" << s1.evalCounter << endl;
                *logger.warn() << "  Worker2: id=" << s2.workerId << " pid=" << s2.workerPid
                              << " eval=" << s2.evalCounter << endl;

                // Compute ULP differences for all key values
                int64_t posULP[3], velULP[3], accULP[3], quatULP[4];
                ulpDouble3(s1.pos, s2.pos, posULP);
                ulpDouble3(s1.vel, s2.vel, velULP);
                ulpDouble3(s1.acc, s2.acc, accULP);
                for (int i = 0; i < 4; ++i) quatULP[i] = computeULPDouble(s1.quat[i], s2.quat[i]);

                int64_t alphaULP = computeULPDouble(s1.alpha, s2.alpha);
                int64_t betaULP = computeULPDouble(s1.beta, s2.beta);
                int64_t vRelWindULP = computeULPDouble(s1.vRelWind, s2.vRelWind);
                int64_t CL_ULP = computeULPDouble(s1.CL, s2.CL);
                int64_t CD_ULP = computeULPDouble(s1.CD, s2.CD);

                *logger.warn() << "  Position: [" << s1.pos[0] << ", " << s1.pos[1] << ", " << s1.pos[2] << "]"
                              << " vs [" << s2.pos[0] << ", " << s2.pos[1] << ", " << s2.pos[2] << "]" << endl;
                *logger.warn() << "    ULP: [" << posULP[0] << ", " << posULP[1] << ", " << posULP[2] << "]" << endl;

                *logger.warn() << "  Velocity: [" << s1.vel[0] << ", " << s1.vel[1] << ", " << s1.vel[2] << "]"
                              << " vs [" << s2.vel[0] << ", " << s2.vel[1] << ", " << s2.vel[2] << "]" << endl;
                *logger.warn() << "    ULP: [" << velULP[0] << ", " << velULP[1] << ", " << velULP[2] << "]" << endl;

                *logger.warn() << "  Acceleration: [" << s1.acc[0] << ", " << s1.acc[1] << ", " << s1.acc[2] << "]"
                              << " vs [" << s2.acc[0] << ", " << s2.acc[1] << ", " << s2.acc[2] << "]" << endl;
                *logger.warn() << "    ULP: [" << accULP[0] << ", " << accULP[1] << ", " << accULP[2] << "]" << endl;

                *logger.warn() << "  Quaternion: [" << s1.quat[0] << ", " << s1.quat[1] << ", "
                              << s1.quat[2] << ", " << s1.quat[3] << "]" << endl;
                *logger.warn() << "         vs: [" << s2.quat[0] << ", " << s2.quat[1] << ", "
                              << s2.quat[2] << ", " << s2.quat[3] << "]" << endl;
                *logger.warn() << "    ULP: [" << quatULP[0] << ", " << quatULP[1] << ", "
                              << quatULP[2] << ", " << quatULP[3] << "]" << endl;

                *logger.warn() << "  Alpha: " << s1.alpha << " vs " << s2.alpha
                              << " (diff=" << std::scientific << (s1.alpha - s2.alpha)
                              << std::defaultfloat << ", " << alphaULP << " ULPs)" << endl;
                *logger.warn() << "  Beta: " << s1.beta << " vs " << s2.beta
                              << " (diff=" << std::scientific << (s1.beta - s2.beta)
                              << std::defaultfloat << ", " << betaULP << " ULPs)" << endl;
                *logger.warn() << "  VrelWind: " << s1.vRelWind << " vs " << s2.vRelWind
                              << " (diff=" << std::scientific << (s1.vRelWind - s2.vRelWind)
                              << std::defaultfloat << ", " << vRelWindULP << " ULPs)" << endl;

                *logger.warn() << "  CL: " << s1.CL << " vs " << s2.CL
                              << " (diff=" << std::scientific << (s1.CL - s2.CL)
                              << std::defaultfloat << ", " << CL_ULP << " ULPs)" << endl;
                *logger.warn() << "  CD: " << s1.CD << " vs " << s2.CD
                              << " (diff=" << std::scientific << (s1.CD - s2.CD)
                              << std::defaultfloat << ", " << CD_ULP << " ULPs)" << endl;

                *logger.warn() << "  CL_wing: " << s1.CL_wing << " (left=" << s1.CL_left
                              << " cent=" << s1.CL_cent << " right=" << s1.CL_right << ")" << endl;
                *logger.warn() << "      vs: " << s2.CL_wing << " (left=" << s2.CL_left
                              << " cent=" << s2.CL_cent << " right=" << s2.CL_right << ")" << endl;

                *logger.warn() << "  Moments (Cl,Cm,Cn): [" << s1.Cl << ", " << s1.Cm << ", " << s1.Cn << "]"
                              << " vs [" << s2.Cl << ", " << s2.Cm << ", " << s2.Cn << "]" << endl;

                *logger.warn() << "  QS (dynamic pressure × area): " << s1.QS << " vs " << s2.QS << endl;

                *logger.warn() << "  Trig values: cos(α)=" << s1.cosAlpha << " vs " << s2.cosAlpha
                              << ", sin(α)=" << s1.sinAlpha << " vs " << s2.sinAlpha
                              << ", cos(β)=" << s1.cosBeta << " vs " << s2.cosBeta << endl;

                *logger.warn() << "  Commands: pitch=" << s1.pitchCommand << " vs " << s2.pitchCommand
                              << ", roll=" << s1.rollCommand << " vs " << s2.rollCommand
                              << ", throttle=" << s1.throttleCommand << " vs " << s2.throttleCommand << endl;

                *logger.warn() << "  Sim inputs: elev=" << s1.elevator << " vs " << s2.elevator
                              << ", ail=" << s1.aileron << " vs " << s2.aileron
                              << ", rud=" << s1.rudder << " vs " << s2.rudder << endl;

                *logger.warn() << "  RNG: state16=" << s1.rngState16 << " vs " << s2.rngState16
                              << ", state32=" << s1.rngState32 << " vs " << s2.rngState32 << endl;

                // Gust data (6 elements: linear velocity [0-2] + rotational rates [3-5])
                int64_t gustULP[6];
                for (int i = 0; i < 6; ++i) {
                  gustULP[i] = computeULPDouble(s1.gustBody[i], s2.gustBody[i]);
                }
                *logger.warn() << "  Gust Linear [xyz]: [" << s1.gustBody[0] << ", " << s1.gustBody[1] << ", " << s1.gustBody[2] << "]"
                              << " vs [" << s2.gustBody[0] << ", " << s2.gustBody[1] << ", " << s2.gustBody[2] << "]" << endl;
                *logger.warn() << "    ULP: [" << gustULP[0] << ", " << gustULP[1] << ", " << gustULP[2] << "]" << endl;
                *logger.warn() << "  Gust Rotational [pqr]: [" << s1.gustBody[3] << ", " << s1.gustBody[4] << ", " << s1.gustBody[5] << "]"
                              << " vs [" << s2.gustBody[3] << ", " << s2.gustBody[4] << ", " << s2.gustBody[5] << "]" << endl;
                *logger.warn() << "    ULP: [" << gustULP[3] << ", " << gustULP[4] << ", " << gustULP[5] << "]" << endl;

                *logger.warn() << endl;

                // Report first divergence and stop
                return;
              }
            }
          }

          *logger.warn() << "No physics trace divergence found - traces are bitwise identical!" << endl;
          *logger.warn() << endl;

          // MACRO TRAJECTORY COMPARISON (100ms GP evaluation steps)
          *logger.warn() << "=== MACRO TRAJECTORY COMPARISON (100ms GP evaluation steps) ===" << endl;
          size_t macroPathCount = std::min(trace1.aircraftStateList.size(), trace2.aircraftStateList.size());
          for (size_t pathIdx = 0; pathIdx < macroPathCount; ++pathIdx) {
            const auto& traj1 = trace1.aircraftStateList[pathIdx];
            const auto& traj2 = trace2.aircraftStateList[pathIdx];
            size_t stepCount = std::min(traj1.size(), traj2.size());

            *logger.warn() << "Path " << pathIdx << ": " << stepCount << " macro steps (100ms intervals)" << endl;

            for (size_t stepIdx = 0; stepIdx < stepCount; ++stepIdx) {
              const auto& a1 = traj1[stepIdx];
              const auto& a2 = traj2[stepIdx];

              gp_vec3 pos1 = a1.getPosition();
              gp_vec3 pos2 = a2.getPosition();
              gp_vec3 vel1 = a1.getVelocity();
              gp_vec3 vel2 = a2.getVelocity();
              gp_quat q1 = a1.getOrientation();
              gp_quat q2 = a2.getOrientation();

              bool posDiff = (pos1.x() != pos2.x()) || (pos1.y() != pos2.y()) || (pos1.z() != pos2.z());
              bool velDiff = (vel1.x() != vel2.x()) || (vel1.y() != vel2.y()) || (vel1.z() != vel2.z());
              bool quatDiff = (q1.x() != q2.x()) || (q1.y() != q2.y()) || (q1.z() != q2.z()) || (q1.w() != q2.w());

              if (posDiff || velDiff || quatDiff) {
                *logger.warn() << "  DIVERGENCE at step " << stepIdx << " (t=" << (stepIdx * 100) << "ms):" << endl;
                if (posDiff) {
                  *logger.warn() << "    Position: [" << pos1.x() << "," << pos1.y() << "," << pos1.z() << "]"
                                << " vs [" << pos2.x() << "," << pos2.y() << "," << pos2.z() << "]" << endl;
                }
                if (velDiff) {
                  *logger.warn() << "    Velocity: [" << vel1.x() << "," << vel1.y() << "," << vel1.z() << "]"
                                << " vs [" << vel2.x() << "," << vel2.y() << "," << vel2.z() << "]" << endl;
                }
                if (quatDiff) {
                  *logger.warn() << "    Quat: [" << q1.w() << "," << q1.x() << "," << q1.y() << "," << q1.z() << "]"
                                << " vs [" << q2.w() << "," << q2.x() << "," << q2.y() << "," << q2.z() << "]" << endl;
                }
                // Only show first divergence per path to avoid log spam
                break;
              }
            }
          }
          *logger.warn() << endl;
      };

      // Elite tracking for non-demetic mode
      // Single ELITE_STORE log per generation with outcome:
      //   INITIAL - first generation
      //   REPLACED - new elite (different gpHash)
      //   SAME - same elite, same fitness (no-op)
      //   DIVERGED - same elite, different fitness (non-determinism)
      MyGP* currentElite = pop->NthMyGP(pop->bestOfPopulation);
      double currentFitness = currentElite->getFitness();  // Use double to match library's stdFitness
      int currentLength = currentElite->length();
      int currentDepth = currentElite->depth();

      // Generate Lisp-form string of elite GP for hash
      std::ostringstream eliteGpStringStream;
      currentElite->printOn(eliteGpStringStream);
      std::string currentEliteGpString = eliteGpStringStream.str();
      uint64_t currentEliteGpHash = hashString(currentEliteGpString);

      // Determine outcome
      bool isFirstGen = (gLastEliteGpHash == 0);
      bool sameStructure = !isFirstGen && (currentEliteGpHash == gLastEliteGpHash);
      const char* outcome;

      if (isFirstGen) {
        outcome = "INITIAL";
      } else if (!sameStructure) {
        outcome = "REPLACED";
      } else {
        // Same structure - check if fitness matches
        gEliteReevalCount++;
        if (!bitwiseEqual(currentFitness, gLastEliteFitness)) {
          gEliteDivergenceCount++;
          outcome = "DIVERGED";
        } else {
          outcome = "SAME";
        }
      }

      // Single log line per generation (warn level for divergence, info otherwise)
      // Compact format: gpHash=old->new fitness=old->new len=old->new depth=old->new
      if (strcmp(outcome, "DIVERGED") == 0) {
        double delta = currentFitness - gLastEliteFitness;
        *logger.warn() << "ELITE_STORE: gen=" << gen
                       << " gpHash=0x" << std::hex << gLastEliteGpHash << "->0x" << currentEliteGpHash << std::dec
                       << " fitness=" << std::fixed << std::setprecision(6) << gLastEliteFitness << "->" << currentFitness
                       << " len=" << gLastEliteLength << "->" << currentLength
                       << " depth=" << gLastEliteDepth << "->" << currentDepth
                       << " delta=" << std::scientific << delta << std::defaultfloat
                       << " " << outcome << endl;
        if (gEliteTraceEnabled) {
          // Compare physics traces to find divergence point
          compareTraces(gLastEliteEvalResults, gPendingEliteEvalResults,
                       "gen" + std::to_string(gen-1) + " elite",
                       "gen" + std::to_string(gen) + " reeval");
        }
      } else {
        *logger.info() << "ELITE_STORE: gen=" << gen
                       << " gpHash=0x" << std::hex << gLastEliteGpHash << "->0x" << currentEliteGpHash << std::dec
                       << " fitness=" << std::fixed << std::setprecision(6) << gLastEliteFitness << "->" << currentFitness
                       << " len=" << gLastEliteLength << "->" << currentLength
                       << " depth=" << gLastEliteDepth << "->" << currentDepth
                       << " " << outcome << endl;
      }

      // Update baseline
      gLastEliteFitness = currentFitness;
      gLastEliteGpHash = currentEliteGpHash;
      gLastEliteLength = currentLength;
      gLastEliteDepth = currentDepth;
      gLastEliteNumPaths = bestOfEvalResults.pathList.size();
      gLastEliteScenarioHash = computeScenarioHash(bestOfEvalResults.scenarioList);
      if (gEliteTraceEnabled) {
        // Move pending elite trace to gLastEliteEvalResults for next gen comparison
        // (gPendingEliteEvalResults is populated by elite reeval callbacks during the gen)
        if (gPendingEliteEvalResults.physicsTrace.empty()) {
          *logger.warn() << "Elite reeval returned no physics trace - divergence analysis unavailable" << endl;
          *logger.warn() << "  (crrcsim may need rebuild, or isEliteReeval flag not set properly)" << endl;
        }
        gLastEliteEvalResults = std::move(gPendingEliteEvalResults);
        // Clear pending for next generation
        gPendingEliteEvalResults = EvalResults();
      }

      // Per-generation counter summary (useful if run is interrupted)
      *logger.info() << "ELITE_STATUS: gen=" << gen
                     << " reevals=" << gEliteReevalCount
                     << " divergences=" << gEliteDivergenceCount << endl;
    }

    // TODO send exit message to workers

    // go ahead and dump out the best of the best
    // ofstream bestGP("best.dat");
    // pop->NthMyGP(pop->bestOfPopulation)->save(bestGP);

    *logger.info() << "GP complete!" << endl;

    // Print elite tracking summary
    *logger.info() << "Elite tracking summary:" << endl;
    *logger.info() << "  Elite re-evals: " << gEliteReevalCount << endl;
    *logger.info() << "  Divergences: " << gEliteDivergenceCount << endl;

    if (gEliteDivergenceCount > 0) {
      *logger.warn() << "WARNING: " << gEliteDivergenceCount << " divergences detected - "
                     << "non-determinism in simulation!" << endl;
    }

    // Clean up final population to prevent memory leak
    delete pop;
    pop = nullptr;

    // Clean up accumulated evaluation results to prevent memory leak
    bestOfEvalResults.clear();
    aggregatedEvalResults.clear();
  }

  uint64_t totalRuns = globalSimRunCounter.load(std::memory_order_relaxed);
  gp_scalar durationSec = std::chrono::duration<gp_scalar>(std::chrono::steady_clock::now() - runStartTime).count();
  gp_scalar simsPerSec = (durationSec > static_cast<gp_scalar>(0.0f)) ? static_cast<gp_scalar>(totalRuns) / durationSec : static_cast<gp_scalar>(0.0f);
  bout << "#SimRuns " << totalRuns << " DurationSec " << durationSec << std::endl;
  bout.flush();
  *logger.info() << std::fixed << std::setprecision(2)
                 << "Simulation throughput: " << totalRuns << " runs in "
                 << durationSec << "s (" << simsPerSec << " sims/s)" << std::defaultfloat << std::endl;
}

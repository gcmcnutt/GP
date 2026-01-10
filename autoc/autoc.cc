
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
extern std::vector<std::vector<Path>> generationPaths;
extern EvalResults bestOfEvalResults;
extern std::ofstream fout;
extern std::ofstream bout;
extern std::atomic_ulong nanDetector;

std::vector<std::vector<Path>> generationPaths;
std::vector<ScenarioDescriptor> generationScenarios;
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
static gp_scalar gLastEliteFitness = std::numeric_limits<gp_scalar>::infinity();
static std::string gLastEliteHash;
static EvalResults gLastEliteTrace;
static std::string gLastEliteKey;
static uint64_t gLastEliteGpHash = 0;
static std::string gLastEliteGpString;  // Lisp form of elite GP for exact comparison
static int gLastEliteLength = std::numeric_limits<int>::max();
static int gLastEliteDepth = std::numeric_limits<int>::max();
static int gEliteDivergenceCount = 0;
static int gEliteCheckCount = 0;
static int gEliteNotCapturedCount = 0;  // Count of generations where elite wasn't re-evaluated
static int gEliteReevalCount = 0;       // Count of successful elite re-evaluations
static int gEliteWorseReplacementCount = 0;  // Count of elite replaced with worse fitness (different structure)
static EvalResults gPriorEliteEval;
static gp_scalar gPriorEliteFitness = std::numeric_limits<gp_scalar>::infinity();
static bool gPriorEliteCaptured = false;
static size_t gLastEliteNumPaths = 0;  // Number of paths in elite evaluation
static uint64_t gLastEliteScenarioHash = 0;  // Hash of scenario metadata
static bool gLastEliteWasBakeoff = false;  // Whether baseline elite was evaluated in bakeoff mode
static std::mutex gPriorEliteMutex;

// Store ALL eval results during current generation, keyed by gpHash
// This allows us to look up the elite's trace after selection completes
static std::map<uint64_t, EvalResults> gCurrentGenEvalResults;
static std::mutex gCurrentGenEvalResultsMutex;
static std::mutex gEliteTrackingMutex;  // Protects all gLastElite* variables for reads/writes

// Global elite tracking - persists across generations (population objects are recreated)
// This is separate from library's bestOfPopulation which tracks single-scenario fitness
static int gGlobalEliteIndex = -1;
static gp_scalar gGlobalEliteFitness = std::numeric_limits<gp_scalar>::infinity();

static bool bitwiseEqual(gp_scalar a, gp_scalar b) {
  return std::memcmp(&a, &b, sizeof(gp_scalar)) == 0;
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

  const ExtraConfig& extraCfg = ConfigManager::getExtraConfig();
  int windScenarioCount = std::max(extraCfg.windScenarioCount, 1);
  unsigned int seedBase = static_cast<unsigned int>(extraCfg.windSeedBase);
  unsigned int seedStride = sanitizeStride(extraCfg.windSeedStride);

  if (basePaths.empty()) {
    ScenarioDescriptor scenario;
    scenario.pathList = basePaths;
    scenario.pathVariantIndex = -1;
    scenario.windVariantIndex = 0;
    scenario.windSeed = seedBase;
    scenario.windScenarios.push_back({seedBase, 0});
    generationScenarios.push_back(std::move(scenario));
    return;
  }

  const GPVariables& gpCfg = ConfigManager::getGPConfig();

  if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0) {
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
        windScenario.windSeed = seedBase + seedStride * static_cast<unsigned int>(windIdx);
        windScenario.windVariantIndex = windIdx;
        scenario.windScenarios.push_back(windScenario);
      }
      if (!scenario.windScenarios.empty()) {
        scenario.windSeed = scenario.windScenarios.front().windSeed;
        scenario.windVariantIndex = scenario.windScenarios.front().windVariantIndex;
      } else {
        scenario.windSeed = seedBase;
        scenario.windVariantIndex = 0;
      }
      generationScenarios.push_back(std::move(scenario));
    }
  } else {
    // NON-DEMETIC MODE: Create ONE scenario containing ALL path variants × ALL winds.
    // Each individual evaluates on the same complete test suite for fair comparison.
    ScenarioDescriptor scenario;
    scenario.pathVariantIndex = -1;  // Multiple paths, no single variant

    // Add all combinations of paths × winds
    for (int windIdx = 0; windIdx < windScenarioCount; ++windIdx) {
      WindScenarioConfig windScenario;
      windScenario.windSeed = seedBase + seedStride * static_cast<unsigned int>(windIdx);
      windScenario.windVariantIndex = windIdx;
      scenario.windScenarios.push_back(windScenario);

      for (size_t pathIdx = 0; pathIdx < basePaths.size(); ++pathIdx) {
        scenario.pathList.push_back(basePaths[pathIdx]);
      }
    }

    if (!scenario.windScenarios.empty()) {
      scenario.windSeed = scenario.windScenarios.front().windSeed;
      scenario.windVariantIndex = scenario.windScenarios.front().windVariantIndex;
    } else {
      scenario.windSeed = seedBase;
      scenario.windVariantIndex = 0;
    }
    generationScenarios.push_back(std::move(scenario));
  }

  if (generationScenarios.empty()) {
    ScenarioDescriptor scenario;
    scenario.pathVariantIndex = -1;
    scenario.windVariantIndex = 0;
    scenario.windSeed = seedBase;
    scenario.pathList = basePaths;
    scenario.windScenarios.push_back({seedBase, 0});
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
  const GPVariables& gpCfg = ConfigManager::getGPConfig();
  int scenarioCount = std::max<int>(generationScenarios.size(), 1);

  if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0) {
    int demeIndex = individualIndex / gpCfg.DemeSize;
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
    GPExitSystem("scenarioForIndex", "generationScenarios is empty - this should never happen!");
  }
  int clampedIndex = ((scenarioIndex % static_cast<int>(generationScenarios.size())) + static_cast<int>(generationScenarios.size())) % static_cast<int>(generationScenarios.size());
  return generationScenarios[clampedIndex];
}

void warnIfScenarioMismatch() {
  const GPVariables& gpCfg = ConfigManager::getGPConfig();
  if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0) {
    int demeCount = gpCfg.PopulationSize / gpCfg.DemeSize;
    if (demeCount > 0 && generationScenarios.size() < static_cast<size_t>(demeCount)) {
      *logger.warn() << "Scenario count (" << generationScenarios.size()
                     << ") is smaller than deme count (" << demeCount
                     << "); scenarios will be reused across demes." << endl;
    }
  }
}

void clearEvalResults(EvalResults& result) {
  result.gp.clear();
  result.gpHash = 0;
  result.crashReasonList.clear();
  result.pathList.clear();
  result.aircraftStateList.clear();
  result.scenario = ScenarioMetadata();
  result.scenarioList.clear();
#ifdef PHYSICS_TRACE_ENABLED
  result.debugSamples.clear();
  result.physicsTrace.clear();
#endif
  result.workerId = -1;
  result.workerPid = 0;
  result.workerEvalCounter = 0;
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
#ifdef PHYSICS_TRACE_ENABLED
  dest.debugSamples.insert(dest.debugSamples.end(), src.debugSamples.begin(), src.debugSamples.end());
  dest.physicsTrace.insert(dest.physicsTrace.end(), src.physicsTrace.begin(), src.physicsTrace.end());
#endif
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
      const GPVariables& gpCfg = ConfigManager::getGPConfig();
      std::vector<int> candidateIndices;
      candidateIndices.reserve(static_cast<size_t>(std::max(1, containerSize())));

      if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0) {
        int demeSize = gpCfg.DemeSize;
        if (demeSize <= 0) {
          demeSize = 1;
        }
        for (int start = 0; start < containerSize(); start += demeSize) {
          int end = std::min(start + demeSize, containerSize());
          gp_scalar bestFitness = std::numeric_limits<gp_scalar>::infinity();
          int bestIndex = start;
          for (int idx = start; idx < end; ++idx) {
            MyGP* candidate = NthMyGP(idx);
            if (!candidate) {
              continue;
            }
            gp_scalar fitness = static_cast<gp_scalar>(candidate->getFitness());
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
        candidateIndices.push_back(bestOfPopulation);
      }

      if (candidateIndices.empty()) {
        candidateIndices.push_back(bestOfPopulation);
      }

      MyGP* best = nullptr;
      gp_scalar bestAggregatedFitness = std::numeric_limits<gp_scalar>::infinity();

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
          clearEvalResults(aggregatedEvalResults);
          aggregatedEvalResults.scenario.pathVariantIndex = -1;
          aggregatedEvalResults.scenario.windVariantIndex = -1;
          aggregatedEvalResults.scenario.windSeed = 0;
          activeEvalCollector = &aggregatedEvalResults;

          gp_scalar cumulativeFitness = 0.0f;
          size_t scenariosEvaluated = 0;

          for (size_t scenarioIdx = 0; scenarioIdx < generationScenarios.size(); ++scenarioIdx) {
            candidate->setScenarioIndex(static_cast<int>(scenarioIdx));
            threadPool->enqueue([candidate](WorkerContext& context) {
              candidate->evalTask(context);
            });
            threadPool->wait_for_tasks();
            cumulativeFitness += static_cast<gp_scalar>(candidate->getFitness());
            ++scenariosEvaluated;
          }

          candidate->setScenarioIndex(originalScenarioIndex);
          candidate->setBakeoffMode(false);

          gp_scalar averageFitness = scenariosEvaluated > 0
            ? cumulativeFitness / static_cast<gp_scalar>(scenariosEvaluated)
            : std::numeric_limits<gp_scalar>::infinity();

          // If this candidate is the prior elite, detect drift even if it loses the top slot.
          std::ostringstream candStream;
          candidate->save(candStream);
          std::string candHash = candStream.str();
          if (!gLastEliteHash.empty() && candHash == gLastEliteHash && !bitwiseEqual(averageFitness, gLastEliteFitness)) {
            gp_scalar delta = averageFitness - gLastEliteFitness;
            *logger.warn() << "ELITE_TRACE: prior elite re-evaluated with different fitness before selection"
                           << " prev=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                           << " current=" << averageFitness
                           << " delta=" << std::scientific << delta << std::defaultfloat << endl;
          }

          // Restore previous logging state before moving on.
          enableDeterministicTestLogging.store(prevLogging, std::memory_order_relaxed);

          if (averageFitness < bestAggregatedFitness) {
            bestAggregatedFitness = averageFitness;
            bestOfEvalResults = aggregatedEvalResults;
            best = candidate;

            // Store bakeoff winner's trace for determinism checking
            // This is the only place where we capture elite's aggregated multi-scenario results
            if (bestOfEvalResults.gpHash != 0) {
              std::lock_guard<std::mutex> lock(gCurrentGenEvalResultsMutex);
              gCurrentGenEvalResults[bestOfEvalResults.gpHash] = bestOfEvalResults;
              *logger.info() << "ELITE_STORE: Stored bakeoff winner gpHash=0x" << std::hex << bestOfEvalResults.gpHash << std::dec
                             << " paths=" << bestOfEvalResults.pathList.size()
                             << " (map now has " << gCurrentGenEvalResults.size() << " entries)" << endl;
            }
          }
        }

        activeEvalCollector = previousCollector;

        // Compare bakeoff winner against global elite
        // Only update global elite if new winner has better aggregated fitness
        if (best) {
          // Find the index of the best candidate in the population
          int bestIndex = -1;
          for (int idx = 0; idx < containerSize(); ++idx) {
            if (NthMyGP(idx) == best) {
              bestIndex = idx;
              break;
            }
          }

          // Check if bakeoff winner beats the global elite
          if (bestAggregatedFitness < gGlobalEliteFitness) {
            gp_scalar improvement = gGlobalEliteFitness - bestAggregatedFitness;
            *logger.info() << "GLOBAL_ELITE: New champion! fitness=" << std::fixed << std::setprecision(6) << bestAggregatedFitness
                           << " (prev=" << gGlobalEliteFitness << " improvement=" << std::scientific << improvement << std::defaultfloat << ")"
                           << " index=" << bestIndex << endl;
            gGlobalEliteFitness = bestAggregatedFitness;
            gGlobalEliteIndex = bestIndex;
            bestOfPopulation = bestIndex;
          } else if (bestIndex == gGlobalEliteIndex) {
            // Elite was re-evaluated in bakeoff - check for determinism
            if (!bitwiseEqual(bestAggregatedFitness, gGlobalEliteFitness)) {
              gp_scalar delta = bestAggregatedFitness - gGlobalEliteFitness;
              *logger.warn() << "GLOBAL_ELITE_DRIFT: Elite re-evaluated with different aggregated fitness!"
                             << " prev=" << std::fixed << std::setprecision(6) << gGlobalEliteFitness
                             << " current=" << bestAggregatedFitness
                             << " delta=" << std::scientific << delta << std::defaultfloat << endl;
            } else {
              *logger.info() << "GLOBAL_ELITE: Elite re-evaluated, fitness unchanged at " << std::fixed << std::setprecision(6) << gGlobalEliteFitness << endl;
            }
            // Elite remains, keep bestOfPopulation pointing to it
            bestOfPopulation = gGlobalEliteIndex;
          } else {
            // New bakeoff winner didn't beat global elite
            *logger.info() << "GLOBAL_ELITE: Bakeoff winner (fitness=" << std::fixed << std::setprecision(6) << bestAggregatedFitness
                           << ") did not beat global elite (fitness=" << gGlobalEliteFitness << ")" << endl;
            // Keep bestOfPopulation pointing to global elite for migration purposes
            bestOfPopulation = gGlobalEliteIndex;
          }

          // Propagate the global elite to all demes by replacing the worst
          // individual in each deme. This spreads the generalist's genes across
          // all path specializations for the next generation.
          // Use gGlobalEliteFitness since that's what bestOfPopulation now points to.
          if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0 &&
              std::isfinite(gGlobalEliteFitness) && gGlobalEliteIndex >= 0) {
            int demeSize = gpCfg.DemeSize;
            int demesPropagated = 0;
            // Config is expressed as percent (0..100)
            gp_scalar migProb = std::clamp(static_cast<gp_scalar>(gpCfg.DemeticMigProbability) / static_cast<gp_scalar>(100.0f),
                                           static_cast<gp_scalar>(0.0f),
                                           static_cast<gp_scalar>(1.0f));
            std::uniform_real_distribution<gp_scalar> dist(static_cast<gp_scalar>(0.0f), static_cast<gp_scalar>(1.0f));
            static thread_local std::mt19937 rng(std::random_device{}());

            for (int demeStart = 0; demeStart < containerSize(); demeStart += demeSize) {
              // Per-deme coin flip: only migrate if roll is below probability
              if (dist(rng) > migProb) {
                continue;
              }

              int demeEnd = std::min(demeStart + demeSize, containerSize());

              // Find worst individual in this deme (highest fitness = worst)
              int worstInDeme = demeStart;
              gp_scalar worstFitness = std::numeric_limits<gp_scalar>::lowest();

              for (int idx = demeStart; idx < demeEnd; ++idx) {
                MyGP* candidate = NthMyGP(idx);
                if (!candidate) {
                  continue;
                }
                gp_scalar fitness = static_cast<gp_scalar>(candidate->getFitness());
                if (!std::isfinite(fitness) || fitness > worstFitness) {
                  worstFitness = fitness;
                  worstInDeme = idx;
                }
              }

              // Don't replace the global elite itself (elitism protection)
              if (worstInDeme != gGlobalEliteIndex) {
                // Clone the global elite and replace the worst
                MyGP* globalElite = NthMyGP(gGlobalEliteIndex);
                if (globalElite) {
                  MyGP* clone = (MyGP*)&globalElite->duplicate();
                  // The clone will have the same fitness as the elite, but it's
                  // not THE elite - it will be re-evaluated next generation
                  put(worstInDeme, *clone);
                  demesPropagated++;
                }
              }
            }

            bout << "# Propagated global elite (fitness="
                 << gGlobalEliteFitness << ") to "
                 << demesPropagated << " demes" << endl;
            bout.flush();
          }
        }
      }

      // Use global elite for storage (not necessarily the bakeoff winner if it didn't beat prior elite)
      MyGP* eliteForStorage = (gGlobalEliteIndex >= 0) ? NthMyGP(gGlobalEliteIndex) : best;
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
      request.SetBucket(ConfigManager::getExtraConfig().s3Bucket);

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
  const GPVariables& gpCfg = ConfigManager::getGPConfig();
  int scenarioCount = std::max<int>(generationScenarios.size(), 1);

  // Assign scenario index to each individual based on demetic grouping
  for (int n = 0; n < containerSize(); ++n) {
    MyGP* current = NthMyGP(n);
#if GPINTERNALCHECK
    if (!current) {
      GPExitSystem("MyPopulation::evaluate", "Member of population is NULL");
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
  evalData.isEliteReeval = false;  // Will be set to true below if this is an elite re-evaluation
  evalData.pathList = scenario.pathList;
  evalData.scenario.scenarioSequence = scenarioSequence;
  evalData.scenario.bakeoffSequence = 0;
  evalData.scenarioList.clear();
  evalData.scenarioList.reserve(scenario.pathList.size());
  bool isBakeoff = bakeoffMode;
  bool enableLogging = enableDeterministicTestLogging.load(std::memory_order_relaxed);

  const GPVariables& gpCfg = ConfigManager::getGPConfig();

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

    if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0) {
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
      // pathList is organized as [wind0:path0, wind0:path1, ..., wind1:path0, wind1:path1, ...]
      size_t numWindScenarios = scenario.windScenarios.size();
      if (numWindScenarios > 0 && scenario.pathList.size() > 0) {
        size_t numBasePaths = scenario.pathList.size() / numWindScenarios;
        size_t windIdx = idx / numBasePaths;
        size_t pathIdx = idx % numBasePaths;

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

  // Check if this is an elite re-evaluation BEFORE sending to worker
  // This allows crrcsim to enable detailed aero logging during the evaluation
  // IMPORTANT: Must check BOTH GP hash AND scenario configuration (for deme mode)
  {
    std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
    size_t currentNumPaths = evalData.pathList.size();
    uint64_t currentScenarioHash = computeScenarioHash(evalData.scenarioList);

    bool gpHashMatches = (gLastEliteGpHash != 0) && (evalData.gpHash == gLastEliteGpHash);
    bool pathCountMatches = (currentNumPaths == gLastEliteNumPaths);
    bool scenarioMatches = (currentScenarioHash == gLastEliteScenarioHash);
    bool bakeoffMatches = (isBakeoff == gLastEliteWasBakeoff);

    // Only flag as elite re-eval if ALL conditions match AND both are bakeoff evaluations
    if (gpHashMatches && pathCountMatches && scenarioMatches && bakeoffMatches && isBakeoff) {
      evalData.isEliteReeval = true;
    }
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

  // NOTE: We no longer store every evaluation in gCurrentGenEvalResults.
  // This was causing massive memory usage (1.5GB peak from copying physicsTrace 758k times).
  // Instead, we only store the elite's EvalResults during bakeoff (see line ~1855).
  // The map is used solely for elite determinism checking, not for all individuals.

  // Capture re-evaluation of prior generation elite (same GP structure AND same path set)
  // CRITICAL: Must verify GP structure, scenario/path set, AND bakeoff mode match
  // We only compare bakeoff-to-bakeoff evaluations (same fitness calculation code path)
  bool isEliteReeval = false;
  {
    std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
    size_t currentNumPaths = context.evalResults.pathList.size();
    uint64_t currentScenarioHash = computeScenarioHash(context.evalResults.scenarioList);

    bool gpMatches = !gLastEliteGpString.empty() && (currentGpString == gLastEliteGpString);
    bool pathCountMatches = (currentNumPaths == gLastEliteNumPaths);
    bool scenarioMatches = (currentScenarioHash == gLastEliteScenarioHash);
    bool bakeoffMatches = (isBakeoff == gLastEliteWasBakeoff);

    // Only flag as elite re-eval if ALL conditions match AND both are bakeoff evaluations
    if (gpMatches && pathCountMatches && scenarioMatches && bakeoffMatches && isBakeoff) {
      isEliteReeval = true;
      evalData.isEliteReeval = true;  // Enable detailed logging in worker
    } else if (gpMatches && (!pathCountMatches || !scenarioMatches)) {
      // Same GP but different path set - expected in demetic mode where elite is
      // evaluated on single deme scenario (3 paths) but baseline is from bakeoff (9 paths)
      // Only log at debug level since this is normal behavior
      *logger.debug() << "ELITE_DEME_EVAL: Elite evaluated on deme scenario (not bakeoff)"
                     << " paths=" << currentNumPaths << " vs " << gLastEliteNumPaths
                     << " scenarioHash=0x" << std::hex << currentScenarioHash
                     << " vs 0x" << gLastEliteScenarioHash << std::dec << endl;
    } else if (gpMatches && pathCountMatches && scenarioMatches && !bakeoffMatches) {
      // Same GP and path set, but different evaluation mode (tournament vs bakeoff)
      *logger.debug() << "ELITE_REEVAL_SKIP: Same GP but different eval mode"
                      << " current=" << (isBakeoff ? "bakeoff" : "tournament")
                      << " baseline=" << (gLastEliteWasBakeoff ? "bakeoff" : "tournament") << endl;
    }
  }

  // Note: Can't capture fitness here yet - it's computed later in the loop
  // We'll capture it after stdFitness is finalized (after line 1141)

  // Compute the fitness results for each path
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

    // compute this path fitness
    gp_scalar localFitness = 0.0f;
    int stepIndex = 0; // where are we on the flight path?

    // SIMPLIFIED error accumulators
    gp_scalar waypoint_distance_sum = 0.0f;        // Distance from current waypoint
    gp_scalar cross_track_error_sum = 0.0f;        // Lateral deviation from path
    gp_scalar movement_direction_error_sum = 0.0f; // Move along path tangent
    gp_scalar throttle_energy_sum = 0.0f;          // Minimize throttle usage
    int simulation_steps = 0;

    // now walk next steps of actual path
    while (++stepIndex < aircraftState.size()) {
      auto& stepAircraftState = aircraftState.at(stepIndex);
      int rawPathIndex = stepAircraftState.getThisPathIndex();
      int pathIndex = std::clamp(rawPathIndex, 0, static_cast<int>(path.size()) - 1);
      const Path& currentPathPoint = path.at(pathIndex);
      PathFrame frame = computePathFrame(path, pathIndex);
      int nextIndex = std::min(pathIndex + 1, static_cast<int>(path.size()) - 1);

      gp_vec3 aircraftPosition = stepAircraftState.getPosition();
      gp_quat craftOrientation = stepAircraftState.getOrientation();
      gp_vec3 aircraftUp = craftOrientation * gp_vec3::UnitZ();

      // ====================================================================
      // METRIC 1: Waypoint distance (reach current target on time)
      // ====================================================================
      gp_scalar waypointDistance = (currentPathPoint.start - aircraftPosition).norm();
      gp_scalar waypointDistancePercent = (waypointDistance * static_cast<gp_scalar>(100.0f)) /
                                           (static_cast<gp_scalar>(2.0f) * SIM_PATH_RADIUS_LIMIT);
      waypoint_distance_sum += pow(waypointDistancePercent, WAYPOINT_DISTANCE_WEIGHT);

      // ====================================================================
      // METRIC 2: Cross-track error (stay near path centerline)
      // ====================================================================
      gp_vec3 craftOffset = aircraftPosition - currentPathPoint.start;
      gp_vec3 lateral = craftOffset - frame.tangent * craftOffset.dot(frame.tangent);
      gp_scalar crossTrackMagnitude = lateral.norm();
      gp_scalar crossTrackPercent = 0.0f;
      if (SIM_PATH_RADIUS_LIMIT > 0.0f) {
        crossTrackPercent = (crossTrackMagnitude * static_cast<gp_scalar>(100.0f)) / SIM_PATH_RADIUS_LIMIT;
        cross_track_error_sum += pow(crossTrackPercent, CROSS_TRACK_WEIGHT);
      }

      // ====================================================================
      // METRIC 2: Movement direction alignment (move along path efficiently)
      // ====================================================================
      gp_scalar movementDirectionError = 0.0f;
      if (stepIndex > 1) {
        gp_vec3 aircraft_movement = stepAircraftState.getPosition() - aircraftState.at(stepIndex-1).getPosition();
        gp_scalar aircraft_distance = aircraft_movement.norm();

        if (aircraft_distance > static_cast<gp_scalar>(0.1f)) {
          // Path direction is the tangent
          gp_scalar direction_alignment = aircraft_movement.normalized().dot(frame.tangent);
          // Clamp to valid range to handle floating-point precision errors at boundaries
          direction_alignment = std::clamp(direction_alignment, static_cast<gp_scalar>(-1.0f), static_cast<gp_scalar>(1.0f));
          // Convert to 0-100 scale: perfect=0, opposite=100
          movementDirectionError = (static_cast<gp_scalar>(1.0f) - direction_alignment) * static_cast<gp_scalar>(50.0f);
          movement_direction_error_sum += pow(movementDirectionError, MOVEMENT_DIRECTION_WEIGHT);
        }
      }

      // ====================================================================
      // METRIC 3: Throttle efficiency (minimize energy consumption)
      // ====================================================================
      gp_scalar current_throttle = stepAircraftState.getThrottleCommand();
      // Throttle is -1 to 1 range, map to 0-1 range for energy (motor off at -1, full power at +1)
      gp_scalar throttleNormalized = (current_throttle + 1.0f) * 0.5f;  // Maps [-1,1] -> [0,1]
      gp_scalar throttlePercent = throttleNormalized * static_cast<gp_scalar>(100.0f);
      throttle_energy_sum += pow(throttlePercent, THROTTLE_EFFICIENCY_WEIGHT);

      // ====================================================================
      // Logging variables for compatibility (not used in fitness)
      // ====================================================================
      gp_scalar distanceFromGoal = (currentPathPoint.start - aircraftPosition).norm();
      distanceFromGoal = distanceFromGoal * static_cast<gp_scalar>(100.0f) / (static_cast<gp_scalar>(2.0f) * SIM_PATH_RADIUS_LIMIT);

      gp_vec3 target_direction = path.at(nextIndex).start - currentPathPoint.start;
      if (target_direction.norm() < static_cast<gp_scalar>(1e-5)) {
        target_direction = frame.tangent;
      }
      gp_vec3 aircraft_to_target = currentPathPoint.start - aircraftPosition;
      gp_scalar angle_rad = 0.0f;
      gp_scalar angle_denom = target_direction.norm() * aircraft_to_target.norm();
      if (angle_denom > static_cast<gp_scalar>(1e-5)) {
        angle_rad = std::acos(std::clamp(target_direction.dot(aircraft_to_target) / angle_denom, -1.0f, 1.0f));
      }
      angle_rad = angle_rad * static_cast<gp_scalar>(100.0f) / static_cast<gp_scalar>(M_PI);

      gp_scalar crossTrackPercentForLog = crossTrackPercent;
      gp_scalar oscillationPercentForLog = 0.0f;  // No longer used
      gp_scalar orientationPenaltyForLog = 0.0f;  // No longer used
      gp_scalar movement_efficiency = movementDirectionError;

      simulation_steps++;

      // use the ugly global to communicate best of gen
      if (printEval) {

        if (printHeader) {
          fout << "Scn    Bake   Pth/Wnd:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP  movEffP      qw      qx      qy      qz   vxBody   vyBody   vzBody    alpha     beta   dtheta    dphi   dhome   xtrkP   xOscP orientP\n";
          printHeader = false;
        }

        // convert aircraft_orientaton to euler
        Eigen::Matrix3f rotMatrix = craftOrientation.toRotationMatrix();

        // Extract Euler angles
        // Note: atan2 returns angle in range [-pi, pi]
        gp_vec3 euler;

        // Handle special case near pitch = ±90° (gimbal lock)
        if (std::abs(rotMatrix(2, 0)) > static_cast<gp_scalar>(0.99999f)) {
          // Gimbal lock case
          euler[0] = 0; // Roll becomes undefined, set to zero

          // Determine pitch based on r31 sign
          if (rotMatrix(2, 0) > 0) {
            euler[1] = -static_cast<gp_scalar>(M_PI / 2.0);
            euler[2] = -atan2(rotMatrix(1, 2), rotMatrix(0, 2)); // yaw
          }
          else {
            euler[1] = static_cast<gp_scalar>(M_PI / 2.0);
            euler[2] = atan2(rotMatrix(1, 2), rotMatrix(0, 2)); // yaw
          }
        }
        else {
          // Normal case
          euler[0] = atan2(rotMatrix(2, 1), rotMatrix(2, 2)); // roll (phi)
          euler[1] = -asin(rotMatrix(2, 0));                 // pitch (theta)
          euler[2] = atan2(rotMatrix(1, 0), rotMatrix(0, 0)); // yaw (psi)
        }

        // Calculate body-frame velocity for GP operators
        gp_vec3 velocity_body = stepAircraftState.getOrientation().inverse() *
                                 stepAircraftState.getVelocity();

        // Calculate GP operator values (matching GP/autoc/gp_operators.h)
        gp_scalar alpha_deg = atan2(-velocity_body.z(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);  // GETALPHA
        gp_scalar beta_deg = atan2(velocity_body.y(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);    // GETBETA

        // Calculate angle to target in body frame (GETDTHETA, GETDPHI)
        gp_vec3 craftToTarget = path.at(pathIndex).start - stepAircraftState.getPosition();
        gp_vec3 target_body = stepAircraftState.getOrientation().inverse() * craftToTarget;
        gp_scalar dtheta_deg = atan2(-target_body.z(), target_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);    // GETDTHETA
        gp_scalar dphi_deg = atan2(target_body.y(), target_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);       // GETDPHI

        // Distance to home (GETDHOME)
        gp_vec3 home(0, 0, SIM_INITIAL_ALTITUDE);
        gp_scalar dhome = (home - stepAircraftState.getPosition()).norm();

        // Get raw quaternion
        gp_quat q = craftOrientation;

        char outbuf[1600]; // XXX use c++20
        sprintf(outbuf, "%06llu %06llu %03d/%02d:%04d: %06ld %3d % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 7.4f % 7.4f % 7.4f % 7.4f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 7.2f % 7.2f % 7.2f\n",
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
          euler[2],
          euler[1],
          euler[0],
          stepAircraftState.getRelVel(),
          stepAircraftState.getRollCommand(),
          stepAircraftState.getPitchCommand(),
          stepAircraftState.getThrottleCommand(),
          distanceFromGoal,
          angle_rad,
          movement_efficiency,
          // NEW FIELDS:
          q.w(), q.x(), q.y(), q.z(),                       // Quaternion
          velocity_body.x(), velocity_body.y(), velocity_body.z(),  // Body velocity
          alpha_deg, beta_deg,                              // GETALPHA, GETBETA
          dtheta_deg, dphi_deg,                             // GETDTHETA, GETDPHI
          dhome,                                            // GETDHOME
          crossTrackPercentForLog,
          oscillationPercentForLog,
          orientationPenaltyForLog
        );
        fout << outbuf;
      }
    }

    // ========================================================================
    // SIMPLIFIED FITNESS COMPUTATION
    // ========================================================================
    // Get total path distance (odometer) for normalization
    gp_scalar total_path_distance = path.back().distanceFromStart;

    // Normalize all metrics by total path distance (odometer reading)
    // This ensures paths of different lengths/granularity are compared fairly
    // Guard against near-zero to prevent divide-by-near-zero overflow
    gp_scalar normalization_factor = (total_path_distance > 0.1f) ? total_path_distance : 1.0f;

    gp_scalar normalized_waypoint_distance = (waypoint_distance_sum / normalization_factor);
    gp_scalar normalized_cross_track = (cross_track_error_sum / normalization_factor);
    gp_scalar normalized_movement_direction = (movement_direction_error_sum / normalization_factor);
    gp_scalar normalized_throttle_energy = (throttle_energy_sum / normalization_factor);

    // Sum all components (lower is better)
    localFitness = normalized_waypoint_distance +     // Reach waypoints on time
                   normalized_cross_track +           // Stay near path centerline
                   normalized_movement_direction +    // Move along path tangent
                   normalized_throttle_energy;        // Minimize throttle usage

    if (isnan(localFitness)) {
      nanDetector++;
    }

    // Crash penalty: add large penalty scaled by how much path remains
    if (crashReason != CrashReason::None) {
      gp_scalar fractional_distance_remaining = static_cast<gp_scalar>(1.0f) -
        path.at(aircraftState.back().getThisPathIndex()).distanceFromStart / path.back().distanceFromStart;
      localFitness += SIM_CRASH_PENALTY * fractional_distance_remaining;
    }

    // accumulate the local fitness
    stdFitness += localFitness;
  }

  // normalize
  stdFitness /= context.evalResults.pathList.size();

  // Mark fitness as valid for the GP library
  fitnessValid = 1;

  // NOW capture elite re-evaluation (after fitness is computed)
  if (isEliteReeval) {
    std::lock_guard<std::mutex> lock(gPriorEliteMutex);
    gPriorEliteEval = context.evalResults;
    gPriorEliteFitness = stdFitness;
    gPriorEliteCaptured = true;
    gEliteReevalCount++;

    *logger.info() << "ELITE_REEVAL: Detected re-evaluation of prior elite"
                   << " gpHash=0x" << std::hex << evalData.gpHash << std::dec
                   << " fitness=" << std::fixed << std::setprecision(6) << stdFitness << endl;

    // Check for divergence immediately (don't wait until next elite change)
    {
      std::lock_guard<std::mutex> eliteLock(gEliteTrackingMutex);
      if (!bitwiseEqual(stdFitness, gLastEliteFitness)) {
        gEliteDivergenceCount++;
        *logger.warn() << "ELITE_REEVAL: DIVERGENCE DETECTED during re-evaluation!"
                       << " baseline=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                       << " reeval=" << stdFitness
                       << " delta=" << std::scientific << (stdFitness - gLastEliteFitness)
                       << std::defaultfloat << endl;
      }
    }
  }

}

void newHandler()
{
  cerr << "\nFatal error: Out of memory." << endl;
  exit(1);
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

  // Init GP system with seed from config
  // Handle -1 as time-based seed (GP library doesn't recognize this convention)
  long gpSeed;
  if (ConfigManager::getExtraConfig().gpSeed == -1) {
    gpSeed = static_cast<long>(time(NULL));
  } else {
    gpSeed = static_cast<long>(ConfigManager::getExtraConfig().gpSeed);
  }
  GPInit(1, gpSeed);

  // AWS setup
  Aws::SDKOptions options;
  Aws::InitAPI(options);

  // initialize workers
  threadPool = new ThreadPool(ConfigManager::getExtraConfig());

  // Print the configuration
  *logger.info() << ConfigManager::getGPConfig() << endl;
  *logger.info() << "SimNumPathsPerGen: " << ConfigManager::getExtraConfig().simNumPathsPerGen << endl;
  *logger.info() << "EvalThreads: " << ConfigManager::getExtraConfig().evalThreads << endl;
  *logger.info() << "MinisimProgram: " << ConfigManager::getExtraConfig().minisimProgram << endl;
  *logger.info() << "MinisimPortOverride: " << ConfigManager::getExtraConfig().minisimPortOverride << endl;
  *logger.info() << "EvaluateMode: " << ConfigManager::getExtraConfig().evaluateMode << endl;
  *logger.info() << "WindScenarios: " << ConfigManager::getExtraConfig().windScenarioCount << endl;
  *logger.info() << "WindSeedBase: " << ConfigManager::getExtraConfig().windSeedBase << endl;
  *logger.info() << "WindSeedStride: " << ConfigManager::getExtraConfig().windSeedStride << endl;
  *logger.info() << "GPSeed: " << ConfigManager::getExtraConfig().gpSeed << endl;
  *logger.info() << "RandomPathSeedB: " << ConfigManager::getExtraConfig().randomPathSeedB << endl << endl;

  // Create the adf function/terminal set and print it out.
  createNodeSet(adfNs);
  *logger.info() << adfNs << endl;

  // Open the main output file for the data and statistics file.
  // First set up names for data file.  Remember we should delete the
  // string from the stream, well just a few bytes
  ostringstream strOutFile, strStatFile;
  strOutFile << "data.dat" << ends;
  strStatFile << "data.stc" << ends;
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

  // prime the paths?
  // Handle -1 as time-based seed for path generation
  unsigned int pathSeed;
  if (ConfigManager::getExtraConfig().randomPathSeedB == -1) {
    pathSeed = static_cast<unsigned int>(time(NULL));
  } else {
    pathSeed = static_cast<unsigned int>(ConfigManager::getExtraConfig().randomPathSeedB);
  }

  generationPaths = generateSmoothPaths(ConfigManager::getExtraConfig().generatorMethod,
                                        ConfigManager::getExtraConfig().simNumPathsPerGen,
                                        SIM_PATH_BOUNDS, SIM_PATH_BOUNDS,
                                        pathSeed);

  // DEBUG: Log segment counts for each path
  std::cout << "\n=== Path Segment Counts (method=" << ConfigManager::getExtraConfig().generatorMethod
            << ", seed=" << pathSeed << ") ===" << std::endl;
  int maxSegments = 0;
  for (size_t i = 0; i < generationPaths.size(); i++) {
    int segCount = generationPaths[i].size();
    std::cout << "  Path " << i << ": " << segCount << " segments" << std::endl;
    if (segCount > maxSegments) maxSegments = segCount;
  }
  std::cout << "  Maximum: " << maxSegments << " segments" << std::endl;
  std::cout << "===" << std::endl << std::endl;

  rebuildGenerationScenarios(generationPaths);
  const int windsPerPath = std::max(ConfigManager::getExtraConfig().windScenarioCount, 1);
  *logger.debug() << "Wind scenarios this generation: paths="
                 << ConfigManager::getExtraConfig().simNumPathsPerGen
                 << " windsPerPath=" << windsPerPath
                 << " dispatchScenarios=" << generationScenarios.size()
                 << " totalEvaluations=" << generationScenarios.size() * windsPerPath
                 << endl;
  warnIfScenarioMismatch();

  if (ConfigManager::getExtraConfig().evaluateMode) {
    // Evaluation mode: use bytecode interpreter with external simulators
    *logger.info() << "Running in bytecode evaluation mode with external simulators..." << endl;
    *logger.info() << "Loading bytecode file: " << ConfigManager::getExtraConfig().bytecodeFile << endl;
    
    // Load the bytecode program
    GPBytecodeInterpreter interpreter;
    if (!interpreter.loadProgram(ConfigManager::getExtraConfig().bytecodeFile)) {
      *logger.error() << "Failed to load bytecode file: " << ConfigManager::getExtraConfig().bytecodeFile << endl;
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
        for (size_t idx = 0; idx < scenario.pathList.size(); ++idx) {
          ScenarioMetadata meta;
          meta.pathVariantIndex = scenario.pathVariantIndex;
          meta.scenarioSequence = scenarioSequence;
          if (isBakeoff) {
            meta.bakeoffSequence = bakeoffPathCounter.fetch_add(1, std::memory_order_relaxed) + 1;
          } else {
            meta.bakeoffSequence = 0;
          }
          if (idx < scenario.windScenarios.size()) {
            meta.windVariantIndex = scenario.windScenarios[idx].windVariantIndex;
            meta.windSeed = scenario.windScenarios[idx].windSeed;
          } else {
            meta.windVariantIndex = scenario.windVariantIndex;
            meta.windSeed = scenario.windSeed;
          }
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
        
        // Use same fitness computation as normal GP evaluation
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

          gp_scalar localFitness = 0;
          int stepIndex = 0;

          // SIMPLIFIED error accumulators
          gp_scalar waypoint_distance_sum = 0.0f;
          gp_scalar cross_track_error_sum = 0.0f;
          gp_scalar movement_direction_error_sum = 0.0f;
          gp_scalar throttle_energy_sum = 0.0f;
          int simulation_steps = 0;
          
          while (++stepIndex < aircraftState.size()) {
            auto& stepAircraftState = aircraftState.at(stepIndex);
            int rawPathIndex = stepAircraftState.getThisPathIndex();
            int pathIndex = std::clamp(rawPathIndex, 0, static_cast<int>(path.size()) - 1);
            const Path& currentPathPoint = path.at(pathIndex);
            PathFrame frame = computePathFrame(path, pathIndex);
            int nextIndex = std::min(pathIndex + 1, static_cast<int>(path.size()) - 1);
            
            gp_vec3 aircraftPosition = stepAircraftState.getPosition();
            gp_quat craftOrientation = stepAircraftState.getOrientation();
            gp_vec3 aircraftUp = craftOrientation * gp_vec3::UnitZ();

            // Waypoint distance
            gp_scalar waypointDistance = (currentPathPoint.start - aircraftPosition).norm();
            gp_scalar waypointDistancePercent = (waypointDistance * static_cast<gp_scalar>(100.0f)) /
                                                 (static_cast<gp_scalar>(2.0f) * SIM_PATH_RADIUS_LIMIT);
            waypoint_distance_sum += pow(waypointDistancePercent, WAYPOINT_DISTANCE_WEIGHT);

            // Cross-track error
            gp_vec3 craftOffset = aircraftPosition - currentPathPoint.start;
            gp_vec3 lateral = craftOffset - frame.tangent * craftOffset.dot(frame.tangent);
            gp_scalar crossTrackMagnitude = lateral.norm();
            gp_scalar crossTrackPercent = 0.0f;
            if (SIM_PATH_RADIUS_LIMIT > static_cast<gp_scalar>(0.0f)) {
              crossTrackPercent = (crossTrackMagnitude * static_cast<gp_scalar>(100.0f)) / SIM_PATH_RADIUS_LIMIT;
              cross_track_error_sum += pow(crossTrackPercent, CROSS_TRACK_WEIGHT);
            }

            // Movement direction alignment
            gp_scalar movementDirectionError = 0.0f;
            if (stepIndex > 1) {
              gp_vec3 aircraft_movement = stepAircraftState.getPosition() - aircraftState.at(stepIndex-1).getPosition();
              gp_scalar aircraft_distance = aircraft_movement.norm();
              if (aircraft_distance > static_cast<gp_scalar>(0.1f)) {
                gp_scalar direction_alignment = aircraft_movement.normalized().dot(frame.tangent);
                // Clamp to valid range to handle floating-point precision errors at boundaries
                direction_alignment = std::clamp(direction_alignment, static_cast<gp_scalar>(-1.0f), static_cast<gp_scalar>(1.0f));
                movementDirectionError = (static_cast<gp_scalar>(1.0f) - direction_alignment) * static_cast<gp_scalar>(50.0f);
                movement_direction_error_sum += pow(movementDirectionError, MOVEMENT_DIRECTION_WEIGHT);
              }
            }

            // Throttle efficiency
            gp_scalar current_throttle = stepAircraftState.getThrottleCommand();
            // Throttle is -1 to 1 range, map to 0-1 range for energy (motor off at -1, full power at +1)
            gp_scalar throttleNormalized = (current_throttle + 1.0f) * 0.5f;  // Maps [-1,1] -> [0,1]
            gp_scalar throttlePercent = throttleNormalized * static_cast<gp_scalar>(100.0f);
            throttle_energy_sum += pow(throttlePercent, THROTTLE_EFFICIENCY_WEIGHT);

            // Logging variables for compatibility (not used in fitness)
            gp_scalar distanceFromGoal = (currentPathPoint.start - aircraftPosition).norm();
            distanceFromGoal = distanceFromGoal * static_cast<gp_scalar>(100.0f) / (static_cast<gp_scalar>(2.0f) * SIM_PATH_RADIUS_LIMIT);

            gp_vec3 target_direction = path.at(nextIndex).start - currentPathPoint.start;
            if (target_direction.norm() < static_cast<gp_scalar>(1e-5f)) {
              target_direction = frame.tangent;
            }
            gp_vec3 aircraft_to_target = currentPathPoint.start - aircraftPosition;
            gp_scalar angle_rad = 0.0f;
            gp_scalar angle_denom = target_direction.norm() * aircraft_to_target.norm();
            if (angle_denom > static_cast<gp_scalar>(1e-5f)) {
              angle_rad = static_cast<gp_scalar>(
                std::acos(std::clamp(target_direction.dot(aircraft_to_target) / angle_denom,
                                     static_cast<gp_scalar>(-1.0f),
                                     static_cast<gp_scalar>(1.0f))));
            }
            angle_rad = angle_rad * static_cast<gp_scalar>(100.0f) / static_cast<gp_scalar>(M_PI);

            gp_scalar crossTrackPercentForLog = crossTrackPercent;
            gp_scalar oscillationPercentForLog = 0.0f;
            gp_scalar orientationPenaltyForLog = 0.0f;
            gp_scalar movement_efficiency = movementDirectionError;

            simulation_steps++;
            
            if (printEval) {
              bestOfEvalResults = context.evalResults;
              
              if (stepIndex == 1) {
                fout << "Pth:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP  movEffP      qw      qx      qy      qz   vxBody   vyBody   vzBody    alpha     beta   dtheta    dphi   dhome   xtrkP   xOscP orientP\n";
              }
              
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

              // Calculate body-frame velocity for GP operators
              gp_vec3 velocity_body = stepAircraftState.getOrientation().inverse() *
                                       stepAircraftState.getVelocity();

              // Calculate GP operator values (matching GP/autoc/gp_operators.h)
              gp_scalar alpha_deg = atan2(-velocity_body.z(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);  // GETALPHA
              gp_scalar beta_deg = atan2(velocity_body.y(), velocity_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);    // GETBETA

              // Calculate angle to target in body frame (GETDTHETA, GETDPHI)
              gp_vec3 craftToTarget = path.at(pathIndex).start - stepAircraftState.getPosition();
              gp_vec3 target_body = stepAircraftState.getOrientation().inverse() * craftToTarget;
              gp_scalar dtheta_deg = atan2(-target_body.z(), target_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);    // GETDTHETA
              gp_scalar dphi_deg = atan2(target_body.y(), target_body.x()) * static_cast<gp_scalar>(180.0 / M_PI);       // GETDPHI

              // Distance to home (GETDHOME)
              gp_vec3 home(0, 0, SIM_INITIAL_ALTITUDE);
              gp_scalar dhome = (home - stepAircraftState.getPosition()).norm();

              // Get raw quaternion
              gp_quat q = craftOrientation;

              char outbuf[1500];
              sprintf(outbuf, "%03d:%04d: %06ld %3d % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 7.4f % 7.4f % 7.4f % 7.4f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 7.2f % 7.2f % 7.2f\n",
                pathVariantIndex, simulation_steps,
                stepAircraftState.getSimTimeMsec(), pathIndex,
                path.at(pathIndex).distanceFromStart,
                path.at(pathIndex).start[0],
                path.at(pathIndex).start[1],
                path.at(pathIndex).start[2],
                stepAircraftState.getPosition()[0],
                stepAircraftState.getPosition()[1],
                stepAircraftState.getPosition()[2],
                euler[2],
                euler[1],
                euler[0],
                stepAircraftState.getRelVel(),
                stepAircraftState.getRollCommand(),
                stepAircraftState.getPitchCommand(),
                stepAircraftState.getThrottleCommand(),
                distanceFromGoal,
                angle_rad,
                movement_efficiency,
                // NEW FIELDS:
                q.w(), q.x(), q.y(), q.z(),                       // Quaternion
                velocity_body.x(), velocity_body.y(), velocity_body.z(),  // Body velocity
                alpha_deg, beta_deg,                              // GETALPHA, GETBETA
                dtheta_deg, dphi_deg,                             // GETDTHETA, GETDPHI
                dhome,                                            // GETDHOME
                crossTrackPercentForLog,
                oscillationPercentForLog,
                orientationPenaltyForLog
              );
              fout << outbuf;
            }
          }
          
          // SIMPLIFIED FITNESS COMPUTATION
          // Get total path distance (odometer) for normalization
          gp_scalar total_path_distance = path.back().distanceFromStart;

          // Normalize all metrics by total path distance (odometer reading)
          // This ensures paths of different lengths/granularity are compared fairly
          // Guard against near-zero to prevent divide-by-near-zero overflow
          gp_scalar normalization_factor = (total_path_distance > 0.1f) ? total_path_distance : 1.0f;

          gp_scalar normalized_waypoint_distance = (waypoint_distance_sum / normalization_factor);
          gp_scalar normalized_cross_track = (cross_track_error_sum / normalization_factor);
          gp_scalar normalized_movement_direction = (movement_direction_error_sum / normalization_factor);
          gp_scalar normalized_throttle_energy = (throttle_energy_sum / normalization_factor);

          localFitness = normalized_waypoint_distance +
                         normalized_cross_track +
                         normalized_movement_direction +
                         normalized_throttle_energy;

          if (isnan(localFitness)) {
            nanDetector++;
          }

          if (crashReason != CrashReason::None) {
            gp_scalar fractional_distance_remaining = static_cast<gp_scalar>(1.0f) -
              path.at(aircraftState.back().getThisPathIndex()).distanceFromStart / path.back().distanceFromStart;
            localFitness += SIM_CRASH_PENALTY * fractional_distance_remaining;
          }
          
          stdFitness += localFitness;
        }
        
        stdFitness /= context.evalResults.pathList.size();
        
        // Mark fitness as valid for the GP library
        fitnessValid = 1;
      }
      
      gp_scalar getFinalFitness() const { return static_cast<gp_scalar>(stdFitness); }
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
    
    clearEvalResults(aggregatedEvalResults);
    aggregatedEvalResults.scenario.pathVariantIndex = -1;
    aggregatedEvalResults.scenario.windVariantIndex = -1;
    aggregatedEvalResults.scenario.windSeed = 0;
    EvalResults* previousCollector = activeEvalCollector;
    activeEvalCollector = &aggregatedEvalResults;
    
    evaluationProgress.store(0);
    gp_scalar cumulativeFitness = 0.0f;
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
    
    gp_scalar averageFitness = scenariosEvaluated > 0 ? cumulativeFitness / static_cast<gp_scalar>(scenariosEvaluated) : static_cast<gp_scalar>(0.0f);
    *logger.info() << "Aggregated results: paths=" << aggregatedEvalResults.pathList.size()
                   << " states=" << aggregatedEvalResults.aircraftStateList.size() << endl;
    bestOfEvalResults = aggregatedEvalResults;
    
    *logger.info() << "Bytecode evaluation complete!" << endl;
    *logger.info() << "Computed fitness: " << averageFitness << endl;
    *logger.info() << "Original GP fitness: " << interpreter.getFitness() << endl;
    
    // Set up S3 storage for results (use gen-1.dmp so renderer can interpret it)
    computedKeyName = startTime + "/gen-1.dmp";
    
    // Store results to S3
    if (bestOfEvalResults.pathList.size() > 0) {
      Aws::S3::Model::PutObjectRequest request;
      request.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
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
    MyPopulation* pop = new MyPopulation(ConfigManager::getGPConfig(), adfNs);
    pop->create();
    *logger.info() << "Ok." << endl;
    pop->createGenerationReport(1, 0, fout, bout, *logger.info());
    logGenerationStats(0);

    // This next for statement is the actual genetic programming system
    // which is in essence just repeated reproduction and crossover loop
    // through all the generations ...
    MyPopulation* newPop = NULL;

    for (int gen = 1; gen <= ConfigManager::getGPConfig().NumberOfGenerations; gen++)
    {
      // Clear eval results from previous generation
      {
        std::lock_guard<std::mutex> lock(gCurrentGenEvalResultsMutex);
        gCurrentGenEvalResults.clear();
      }

      // For this generation, build a smooth path goal
      // Handle -1 as time-based seed for path generation
      unsigned int genPathSeed;
      if (ConfigManager::getExtraConfig().randomPathSeedB == -1) {
        genPathSeed = static_cast<unsigned int>(time(NULL));
      } else {
        genPathSeed = static_cast<unsigned int>(ConfigManager::getExtraConfig().randomPathSeedB);
      }

      generationPaths = generateSmoothPaths(ConfigManager::getExtraConfig().generatorMethod,
                                            ConfigManager::getExtraConfig().simNumPathsPerGen,
                                            SIM_PATH_BOUNDS, SIM_PATH_BOUNDS,
                                            genPathSeed);
      rebuildGenerationScenarios(generationPaths);
      warnIfScenarioMismatch();

      // Create a new generation from the old one by applying the genetic operators
      if (!ConfigManager::getGPConfig().SteadyState)
        newPop = new MyPopulation(ConfigManager::getGPConfig(), adfNs);
      pop->generate(*newPop);

      // Enable evaluation output for this generation
      printEval = true;

      // Switch to new population first to get the correct best individual
      if (!ConfigManager::getGPConfig().SteadyState)
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

      // In-memory trace comparison function
      auto compareTraces = [&](const EvalResults& trace1, const EvalResults& trace2,
                               const std::string& label1, const std::string& label2) {
#ifndef PHYSICS_TRACE_ENABLED
        // Physics trace disabled - skip detailed comparison
        (void)trace1; (void)trace2; (void)label1; (void)label2;
        return;
#else
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
#endif  // PHYSICS_TRACE_ENABLED
      };

      MyGP* currentElite = pop->NthMyGP(pop->bestOfPopulation);
      gp_scalar currentFitness = currentElite->getFitness();
      int currentEliteLength = currentElite->length();
      int currentEliteDepth = currentElite->depth();

      // Compute hash of elite's tree structure to detect if it's truly the same individual
      std::ostringstream eliteStream;
      currentElite->save(eliteStream);
      std::string currentEliteHash = eliteStream.str();

      // Generate Lisp-form string of elite GP for exact comparison
      std::ostringstream eliteGpStringStream;
      currentElite->printOn(eliteGpStringStream);
      std::string currentEliteGpString = eliteGpStringStream.str();

      // Check if elite changed (structure changed, not just population index)
      bool eliteStructureChanged = (currentEliteHash != gLastEliteHash);

      // If we captured a re-eval of the prior elite during this generation, check for drift now.
      if (gPriorEliteCaptured) {
        std::lock_guard<std::mutex> lock(gPriorEliteMutex);
        if (!bitwiseEqual(gPriorEliteFitness, gLastEliteFitness)) {
          gEliteDivergenceCount++;
          *logger.warn() << "ELITE_TRACE gen=" << gen << ": Prior elite drifted on re-eval before selection"
                         << " prev=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                         << " current=" << gPriorEliteFitness
                         << " delta=" << std::scientific << (gPriorEliteFitness - gLastEliteFitness)
                         << " [" << gEliteDivergenceCount << ":" << gEliteCheckCount << " divergences]"
                         << std::defaultfloat << endl;
          compareTraces(gLastEliteTrace, gPriorEliteEval, gLastEliteKey, "<prior-elite-reeval>");
        }
        gPriorEliteCaptured = false;
      } else if (gLastEliteGpHash != 0) {
        gEliteNotCapturedCount++;
        *logger.warn() << "ELITE_TRACE gen=" << gen << ": Prior elite hash "
                       << std::hex << "0x" << gLastEliteGpHash << std::dec
                       << " not captured during re-eval (possible mismatch or missing gp payload)"
                       << " [" << gEliteDivergenceCount << ":" << gEliteCheckCount << " divergences]" << endl;
      }

      // Compute elite GP hash from Lisp string (deterministic, same as evalTask)
      uint64_t currentEliteGpHash = hashString(currentEliteGpString);

      // Look up elite's eval results from current generation
      EvalResults* eliteEvalResults = nullptr;
      {
        std::lock_guard<std::mutex> lock(gCurrentGenEvalResultsMutex);
        auto it = gCurrentGenEvalResults.find(currentEliteGpHash);
        if (it != gCurrentGenEvalResults.end()) {
          eliteEvalResults = &it->second;
        }

        // Free memory: clear all non-elite traces immediately after lookup
        // Keep only the elite's trace to avoid memory leak from physicsTrace data
        if (eliteEvalResults && gCurrentGenEvalResults.size() > 1) {
          EvalResults eliteCopy = *eliteEvalResults;  // Copy before clearing map
          gCurrentGenEvalResults.clear();
          gCurrentGenEvalResults[currentEliteGpHash] = std::move(eliteCopy);
          eliteEvalResults = &gCurrentGenEvalResults[currentEliteGpHash];  // Update pointer
        }
      }

      // Generate trace label for logging (no longer writing files)
      std::string traceLabel;
      if (eliteEvalResults) {
        traceLabel = "gen" + std::to_string(gen)
                   + "-worker" + std::to_string(eliteEvalResults->workerId)
                   + "-pid" + std::to_string(eliteEvalResults->workerPid)
                   + "-eval" + std::to_string(eliteEvalResults->workerEvalCounter);
      } else {
        traceLabel = "gen" + std::to_string(gen) + "-NOT-FOUND";
        *logger.warn() << "ELITE_TRACE gen=" << gen << ": Elite gpHash 0x" << std::hex
                       << currentEliteGpHash << std::dec << " not found in eval results map!" << endl;
      }

      if (gen == 1) {
        // First generation: establish baseline
        std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
        gLastEliteFitness = currentFitness;
        gLastEliteHash = currentEliteHash;
        if (eliteEvalResults) {
          gLastEliteTrace = *eliteEvalResults;
        } else {
          // Fallback to bestOfEvalResults (bakeoff) if not in map
          gLastEliteTrace = bestOfEvalResults;
        }
        gLastEliteKey = traceLabel;
        gLastEliteGpHash = currentEliteGpHash;
        gLastEliteGpString = currentEliteGpString;
        if (eliteEvalResults) {
          gLastEliteNumPaths = eliteEvalResults->pathList.size();
          gLastEliteScenarioHash = computeScenarioHash(eliteEvalResults->scenarioList);
        } else {
          gLastEliteNumPaths = bestOfEvalResults.pathList.size();
          gLastEliteScenarioHash = computeScenarioHash(bestOfEvalResults.scenarioList);
        }
        gLastEliteWasBakeoff = true;  // Gen 1 elite always evaluated in bakeoff
        gLastEliteLength = currentEliteLength;
        gLastEliteDepth = currentEliteDepth;

        *logger.info() << "ELITE_TRACE gen=" << gen << ": Baseline elite captured"
                      << " fitness=" << std::fixed << std::setprecision(6) << currentFitness
                      << " paths=" << gLastEliteNumPaths
                      << " scenarioHash=0x" << std::hex << gLastEliteScenarioHash << std::dec
                      << " trace=" << traceLabel
                      << " [" << gEliteDivergenceCount << ":" << gEliteCheckCount << " divergences]" << endl;
      } else {
        // Subsequent generations: check for changes
        gEliteCheckCount++;

        auto logCounters = [&]() {
          *logger.warn() << "  EliteChecks=" << gEliteCheckCount << " Divergences=" << gEliteDivergenceCount
                        << " (" << std::fixed << std::setprecision(2)
                        << (gEliteCheckCount > 0 ? (100.0 * gEliteDivergenceCount / gEliteCheckCount) : 0.0)
                        << "%)" << endl;
        };

        bool fitnessImproved = currentFitness < gLastEliteFitness;
        bool fitnessWorsened = currentFitness > gLastEliteFitness;
        bool fitnessEqual = bitwiseEqual(currentFitness, gLastEliteFitness);

        if (eliteStructureChanged && fitnessImproved) {
          // New elite is better - update baseline
          *logger.info() << "ELITE_TRACE gen=" << gen << ": Elite replaced (better)"
                        << " old_fitness=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                        << " new_fitness=" << currentFitness
                        << " improvement=" << std::scientific << (gLastEliteFitness - currentFitness)
                        << " [" << gEliteDivergenceCount << ":" << gEliteCheckCount << " divergences]"
                        << std::defaultfloat << endl;

          {
            std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
            gLastEliteFitness = currentFitness;
            gLastEliteHash = currentEliteHash;
            gLastEliteTrace = eliteEvalResults ? *eliteEvalResults : bestOfEvalResults;
            gLastEliteKey = traceLabel;
            gLastEliteGpHash = currentEliteGpHash;
            gLastEliteGpString = currentEliteGpString;
            gLastEliteNumPaths = bestOfEvalResults.pathList.size();
            gLastEliteScenarioHash = computeScenarioHash(bestOfEvalResults.scenarioList);
            gLastEliteWasBakeoff = true;  // Elite replacement always in bakeoff
            gLastEliteLength = currentEliteLength;
            gLastEliteDepth = currentEliteDepth;
          }

          *logger.info() << "ELITE_TRACE gen=" << gen << ": New baseline captured"
                        << " trace=" << traceLabel
                        << " [" << gEliteDivergenceCount << ":" << gEliteCheckCount << " divergences]" << endl;
        } else if (!eliteStructureChanged && fitnessWorsened) {
          // Same elite structure but worse fitness - ELITISM VIOLATION!
          gEliteDivergenceCount++;
          *logger.warn() << "ELITE_TRACE gen=" << gen << ": ELITISM VIOLATION DETECTED!" << endl;
          *logger.warn() << "  Same elite individual re-evaluated with different fitness:" << endl;
          *logger.warn() << "  Previous: fitness=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                        << " trace=" << gLastEliteKey << endl;
          *logger.warn() << "  Current:  fitness=" << std::fixed << std::setprecision(6) << currentFitness
                        << " trace=" << traceLabel << endl;
          *logger.warn() << "  Delta: " << std::scientific << (currentFitness - gLastEliteFitness)
                        << " (" << std::setprecision(2) << (100.0 * (currentFitness - gLastEliteFitness) / gLastEliteFitness)
                        << "%)" << std::defaultfloat << endl;
          *logger.warn() << "  Compare: " << gLastEliteKey << " vs " << traceLabel << endl;
          logCounters();
          *logger.warn() << endl;

          // In-memory trace comparison
          compareTraces(gLastEliteTrace, bestOfEvalResults, gLastEliteKey, traceLabel);

          // Update baseline to current (violated) state
          {
            std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
            gLastEliteFitness = currentFitness;
            gLastEliteHash = currentEliteHash;
            gLastEliteTrace = eliteEvalResults ? *eliteEvalResults : bestOfEvalResults;
            gLastEliteKey = traceLabel;
            gLastEliteGpHash = currentEliteGpHash;
            gLastEliteGpString = currentEliteGpString;
            gLastEliteNumPaths = bestOfEvalResults.pathList.size();
            gLastEliteScenarioHash = computeScenarioHash(bestOfEvalResults.scenarioList);
            gLastEliteWasBakeoff = true;  // Violation detected in bakeoff
            gLastEliteLength = currentEliteLength;
            gLastEliteDepth = currentEliteDepth;
          }
        } else if (!eliteStructureChanged && fitnessImproved) {
          // Same elite structure but better fitness - NON-DETERMINISTIC!
          gEliteDivergenceCount++;
          *logger.warn() << "ELITE_TRACE gen=" << gen << ": ELITISM VIOLATION (fitness paradox)!" << endl;
          *logger.warn() << "  Same elite structure, but improved fitness (non-deterministic!)" << endl;
          *logger.warn() << "  Previous: fitness=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                        << " trace=" << gLastEliteKey << endl;
          *logger.warn() << "  Current:  fitness=" << std::fixed << std::setprecision(6) << currentFitness
                        << " trace=" << traceLabel << endl;
          *logger.warn() << "  Delta: " << std::scientific << (gLastEliteFitness - currentFitness)
                        << " (" << std::setprecision(2) << (100.0 * (currentFitness - gLastEliteFitness) / gLastEliteFitness)
                        << "%)" << std::defaultfloat << endl;
          *logger.warn() << "  Compare: " << gLastEliteKey << " vs " << traceLabel << endl;
          logCounters();
          *logger.warn() << endl;

          // In-memory trace comparison
          compareTraces(gLastEliteTrace, bestOfEvalResults, gLastEliteKey, traceLabel);

          // Update to current (improved) state
          {
            std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
            gLastEliteFitness = currentFitness;
            gLastEliteHash = currentEliteHash;
            gLastEliteTrace = eliteEvalResults ? *eliteEvalResults : bestOfEvalResults;
            gLastEliteKey = traceLabel;
            gLastEliteGpHash = currentEliteGpHash;
            gLastEliteGpString = currentEliteGpString;
            gLastEliteNumPaths = bestOfEvalResults.pathList.size();
            gLastEliteScenarioHash = computeScenarioHash(bestOfEvalResults.scenarioList);
            gLastEliteWasBakeoff = true;  // Paradox case in bakeoff
            gLastEliteLength = currentEliteLength;
            gLastEliteDepth = currentEliteDepth;
          }
        } else if (!eliteStructureChanged && !fitnessWorsened && !fitnessImproved) {
          // Same elite, same fitness - NORMAL DETERMINISTIC CASE
          // Silent update (no logging needed for normal case)
          {
            std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
            gLastEliteTrace = bestOfEvalResults;
            gLastEliteKey = traceLabel;
          }
        } else if (eliteStructureChanged && fitnessEqual) {
          bool shorter = currentEliteLength < gLastEliteLength;
          bool sameLength = currentEliteLength == gLastEliteLength;
          bool shallower = currentEliteDepth < gLastEliteDepth;

          if (shorter || (sameLength && shallower)) {
            *logger.info() << "ELITE_TRACE gen=" << gen << ": Elite replaced (tie-break)"
                          << " fitness=" << std::fixed << std::setprecision(6) << currentFitness
                          << " length " << gLastEliteLength << "->" << currentEliteLength
                          << " depth " << gLastEliteDepth << "->" << currentEliteDepth
                          << " [" << gEliteDivergenceCount << ":" << gEliteCheckCount << " divergences]" << endl;
          } else {
            *logger.warn() << "ELITE_TRACE gen=" << gen << ": Elite replaced with equal fitness but no tie-break win" << endl;
            *logger.warn() << "  Previous: fitness=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                          << " length=" << gLastEliteLength << " depth=" << gLastEliteDepth << endl;
            *logger.warn() << "  Current:  fitness=" << std::fixed << std::setprecision(6) << currentFitness
                          << " length=" << currentEliteLength << " depth=" << currentEliteDepth << endl;
            logCounters();
          }

          // Update baseline (structure changed due to tie or unexpected swap)
          {
            std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
            gLastEliteFitness = currentFitness;
            gLastEliteHash = currentEliteHash;
            gLastEliteTrace = eliteEvalResults ? *eliteEvalResults : bestOfEvalResults;
            gLastEliteKey = traceLabel;
            gLastEliteGpHash = currentEliteGpHash;
            gLastEliteGpString = currentEliteGpString;
            gLastEliteNumPaths = bestOfEvalResults.pathList.size();
            gLastEliteScenarioHash = computeScenarioHash(bestOfEvalResults.scenarioList);
            gLastEliteWasBakeoff = true;  // Elite swap in bakeoff
            gLastEliteLength = currentEliteLength;
            gLastEliteDepth = currentEliteDepth;
          }
        } else if (eliteStructureChanged && fitnessWorsened) {
          // Different individual with worse fitness - in demetic mode this indicates non-determinism
          // because the bakeoff winner should always be at least as good as before
          gEliteWorseReplacementCount++;
          *logger.warn() << "ELITE_TRACE gen=" << gen << ": Elite replaced with worse fitness (different structure)"
                        << " old_fitness=" << std::fixed << std::setprecision(6) << gLastEliteFitness
                        << " new_fitness=" << currentFitness
                        << " delta=" << std::scientific << (currentFitness - gLastEliteFitness)
                        << " length " << gLastEliteLength << "->" << currentEliteLength
                        << " depth " << gLastEliteDepth << "->" << currentEliteDepth
                        << std::defaultfloat << endl;

          {
            std::lock_guard<std::mutex> lock(gEliteTrackingMutex);
            gLastEliteFitness = currentFitness;
            gLastEliteHash = currentEliteHash;
            gLastEliteTrace = eliteEvalResults ? *eliteEvalResults : bestOfEvalResults;
            gLastEliteKey = traceLabel;
            gLastEliteGpHash = currentEliteGpHash;
            gLastEliteGpString = currentEliteGpString;
            gLastEliteNumPaths = bestOfEvalResults.pathList.size();
            gLastEliteScenarioHash = computeScenarioHash(bestOfEvalResults.scenarioList);
            gLastEliteWasBakeoff = true;  // Elite replacement in bakeoff
            if (gLastEliteGpHash == 0) {
              *logger.warn() << "ELITE_TRACE: baseline GP hash is zero; drift detection will be disabled" << endl;
            }
            gLastEliteLength = currentEliteLength;
            gLastEliteDepth = currentEliteDepth;
          }
        }
      }

      // Per-generation elite tracking summary (always output for visibility)
      *logger.info() << "ELITE_STATUS gen=" << gen
                     << " checks=" << gEliteCheckCount
                     << " reevals=" << gEliteReevalCount
                     << " notcaptured=" << gEliteNotCapturedCount
                     << " divergences=" << gEliteDivergenceCount
                     << " worse_replacements=" << gEliteWorseReplacementCount << endl;
    }

    // TODO send exit message to workers

    // go ahead and dump out the best of the best
    // ofstream bestGP("best.dat");
    // pop->NthMyGP(pop->bestOfPopulation)->save(bestGP);

    *logger.info() << "GP complete!" << endl;

    // Print elite tracking summary
    *logger.info() << "Elite tracking summary:" << endl;
    *logger.info() << "  Generations checked: " << gEliteCheckCount << endl;
    *logger.info() << "  Elite re-evals captured: " << gEliteReevalCount << endl;
    *logger.info() << "  Elite not captured: " << gEliteNotCapturedCount << endl;
    *logger.info() << "  Divergences detected: " << gEliteDivergenceCount << endl;
    *logger.info() << "  Worse replacements: " << gEliteWorseReplacementCount << endl;

    // Calculate total anomalies (any sign of non-determinism)
    int totalAnomalies = gEliteDivergenceCount + gEliteWorseReplacementCount;
    if (totalAnomalies > 0) {
      *logger.warn() << "WARNING: " << totalAnomalies << " elite anomalies detected ("
                     << gEliteDivergenceCount << " divergences, "
                     << gEliteWorseReplacementCount << " worse replacements) - "
                     << "check for non-determinism in simulation!" << endl;
    }
    if (gEliteReevalCount == 0 && gEliteCheckCount > 0 && gEliteWorseReplacementCount == 0) {
      *logger.warn() << "WARNING: No elite re-evaluations captured in " << gEliteCheckCount
                     << " generations - determinism verification was limited!" << endl;
    }

    // Clean up final population to prevent memory leak
    delete pop;
    pop = nullptr;

    // Clean up accumulated evaluation results to prevent memory leak
    clearEvalResults(bestOfEvalResults);
    clearEvalResults(aggregatedEvalResults);
    gCurrentGenEvalResults.clear();
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

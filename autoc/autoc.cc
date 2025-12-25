
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
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace std;
namespace logging = boost::log;

Logger logger;
std::atomic_bool printEval = false; // verbose (used for rendering best of population)
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
std::atomic<uint64_t> evaluationProgress{0};
std::atomic<uint64_t> globalScenarioCounter{0};
std::atomic<uint64_t> bakeoffPathCounter{0};
std::atomic<uint64_t> globalSimRunCounter{0};



ThreadPool* threadPool;
std::string computedKeyName;

namespace {

unsigned int sanitizeStride(int stride) {
  int absStride = stride == 0 ? 1 : std::abs(stride);
  return static_cast<unsigned int>(absStride);
}

ScenarioDescriptor& defaultScenario() {
  static ScenarioDescriptor fallback;
  return fallback;
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

  // Create one scenario per path. Each scenario contains ONE path evaluated
  // across ALL wind conditions. This allows demes to specialize on specific
  // path geometries while being robust to wind variations.
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

  if (generationScenarios.empty()) {
    ScenarioDescriptor scenario;
    scenario.pathVariantIndex = -1;
    scenario.windVariantIndex = 0;
    scenario.windSeed = seedBase;
    scenario.pathList = basePaths;
    scenario.windScenarios.push_back({seedBase, 0});
    generationScenarios.push_back(std::move(scenario));
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

  return individualIndex % scenarioCount;
}

const ScenarioDescriptor& scenarioForIndex(int scenarioIndex) {
  if (generationScenarios.empty()) {
    ScenarioDescriptor& fallback = defaultScenario();
    fallback.pathList = generationPaths;
    fallback.windSeed = static_cast<unsigned int>(ConfigManager::getExtraConfig().windSeedBase);
    fallback.pathVariantIndex = -1;
    fallback.windVariantIndex = 0;
    fallback.windScenarios.clear();
    fallback.windScenarios.resize(fallback.pathList.size());
    for (size_t idx = 0; idx < fallback.windScenarios.size(); ++idx) {
      fallback.windScenarios[idx].windSeed = fallback.windSeed;
      fallback.windScenarios[idx].windVariantIndex = static_cast<int>(idx);
    }
    return fallback;
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
  result.crashReasonList.clear();
  result.pathList.clear();
  result.aircraftStateList.clear();
  result.scenario = ScenarioMetadata();
  result.scenarioList.clear();
}

void appendEvalResults(EvalResults& dest, const EvalResults& src) {
  if (dest.gp.empty()) {
    dest.gp = src.gp;
  }
  dest.crashReasonList.insert(dest.crashReasonList.end(), src.crashReasonList.begin(), src.crashReasonList.end());
  dest.pathList.insert(dest.pathList.end(), src.pathList.begin(), src.pathList.end());
  dest.aircraftStateList.insert(dest.aircraftStateList.end(), src.aircraftStateList.begin(), src.aircraftStateList.end());
  dest.scenarioList.insert(dest.scenarioList.end(), src.scenarioList.begin(), src.scenarioList.end());
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
          best->setFitness(bestAggregatedFitness);
          // Find the index of the best candidate in the population
          for (int idx = 0; idx < containerSize(); ++idx) {
            if (NthMyGP(idx) == best) {
              bestOfPopulation = idx;
              break;
            }
          }

          // Propagate the bakeoff winner to all demes by replacing the worst
          // individual in each deme. This spreads the generalist's genes across
          // all path specializations for the next generation.
          if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0 &&
              std::isfinite(bestAggregatedFitness)) {
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

              // Don't replace the bakeoff winner itself (elitism protection)
              if (worstInDeme != bestOfPopulation) {
                // Clone the winner and replace the worst
                put(worstInDeme, best->duplicate());
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

      if (!best) {
        best = NthMyGP(bestOfPopulation);
      }

      // Re-serialize the best GP with updated fitness for storage
      std::vector<char> updatedBuffer;
      boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> updatedOutStream(updatedBuffer);
      best->save(updatedOutStream);
      updatedOutStream.flush();
      bestOfEvalResults.gp = updatedBuffer;
      
      // now put the resulting elements into the S3 object
      Aws::S3::Model::PutObjectRequest request;
      request.SetBucket(ConfigManager::getExtraConfig().s3Bucket);

      // path name is $base/RunDate/gen$gen.dmp
      request.SetKey(computedKeyName);

      std::ostringstream oss;
      boost::archive::text_oarchive oa(oss);
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
  for (int n = 0; n < containerSize(); ++n) {
    MyGP* current = NthMyGP(n);
#if GPINTERNALCHECK
    if (!current) {
      GPExitSystem("MyPopulation::evaluate", "Member of population is NULL");
    }
#endif
    current->setScenarioIndex(computeScenarioIndexForIndividual(n));
  }
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

  EvalData evalData;
  evalData.gp = buffer;
  evalData.pathList = scenario.pathList;
  evalData.scenario.scenarioSequence = scenarioSequence;
  evalData.scenario.bakeoffSequence = 0;
  evalData.scenarioList.clear();
  evalData.scenarioList.reserve(scenario.pathList.size());
  bool isBakeoff = bakeoffMode;
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
  sendRPC(*context.socket, evalData);

  // How did it go?
  context.evalResults = receiveRPC<EvalResults>(*context.socket);
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

    // error accumulators
    gp_scalar distance_error_sum = 0.0f;
    gp_scalar angle_error_sum = 0.0f;
    gp_scalar movement_efficiency_sum = 0.0f;
    gp_scalar cross_track_error_sum = 0.0f;
    gp_scalar cross_track_osc_sum = 0.0f;
    int cross_track_osc_events = 0;
    gp_scalar orientation_alignment_sum = 0.0f;
    gp_scalar prev_cross_track_signed = 0.0f;
    bool has_prev_cross_track = false;
    int simulation_steps = 0;

    // initial states
    gp_scalar roll_prev = aircraftState.at(stepIndex).getRollCommand();
    gp_scalar pitch_prev = aircraftState.at(stepIndex).getPitchCommand();
    gp_scalar throttle_prev = aircraftState.at(stepIndex).getThrottleCommand();

    // now walk next steps of actual path
    while (++stepIndex < aircraftState.size()) {
      auto& stepAircraftState = aircraftState.at(stepIndex);
      int rawPathIndex = stepAircraftState.getThisPathIndex();
      int pathIndex = std::clamp(rawPathIndex, 0, static_cast<int>(path.size()) - 1);
      const Path& currentPathPoint = path.at(pathIndex);
      PathFrame frame = computePathFrame(path, pathIndex);
      int nextIndex = std::min(pathIndex + 1, static_cast<int>(path.size()) - 1);

      gp_vec3 aircraftPosition = stepAircraftState.getPosition();

      // Compute the distance between the aircraft and the goal
      gp_scalar distanceFromGoal = (currentPathPoint.start - aircraftPosition).norm();
      distanceFromGoal = distanceFromGoal * static_cast<gp_scalar>(100.0f) / (static_cast<gp_scalar>(2.0f) * SIM_PATH_RADIUS_LIMIT);

      // Compute vector from me to target using local Frenet frame to avoid gimbal issues
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

      gp_quat craftOrientation = stepAircraftState.getOrientation();
      gp_vec3 aircraftForward = craftOrientation * gp_vec3::UnitX();
      gp_vec3 aircraftUp = craftOrientation * gp_vec3::UnitZ();

      gp_scalar crossTrackPercentForLog = 0.0f;
      gp_scalar oscillationPercentForLog = 0.0f;
      gp_scalar orientationPenaltyForLog = 0.0f;

      // Cross-track error measured in the plane perpendicular to the tangent
      gp_vec3 craftOffset = aircraftPosition - currentPathPoint.start;
      gp_vec3 lateral = craftOffset - frame.tangent * craftOffset.dot(frame.tangent);
      gp_scalar crossTrackSigned = lateral.dot(frame.insideNormal);
      gp_scalar crossTrackMagnitude = lateral.norm();
      if (SIM_PATH_RADIUS_LIMIT > 0.0f) {
        crossTrackPercentForLog = (crossTrackMagnitude * static_cast<gp_scalar>(100.0f)) / SIM_PATH_RADIUS_LIMIT;
        cross_track_error_sum += pow(crossTrackPercentForLog, CROSS_TRACK_WEIGHT);
      }

      gp_scalar oscillationPercent = 0.0f;
      bool oscillationTriggered = false;
      if (has_prev_cross_track) {
        if ((crossTrackSigned > CROSS_TRACK_SIGN_THRESHOLD && prev_cross_track_signed < -CROSS_TRACK_SIGN_THRESHOLD) ||
            (crossTrackSigned < -CROSS_TRACK_SIGN_THRESHOLD && prev_cross_track_signed > CROSS_TRACK_SIGN_THRESHOLD)) {
          gp_scalar denomOsc = std::max(1.0f, static_cast<gp_scalar>(2.0f) * SIM_PATH_RADIUS_LIMIT);
          oscillationPercent = (std::abs(prev_cross_track_signed) + std::abs(crossTrackSigned)) * static_cast<gp_scalar>(100.0f) / denomOsc;
          cross_track_osc_sum += pow(oscillationPercent, CROSS_TRACK_OSC_WEIGHT);
          cross_track_osc_events++;
          oscillationTriggered = true;
        }
      }
      oscillationPercentForLog = oscillationTriggered ? oscillationPercent : 0.0;
      prev_cross_track_signed = crossTrackSigned;
      has_prev_cross_track = true;

      // Quaternion-based orientation alignment (avoid gimbal lock)
      gp_scalar forwardDot = std::clamp(aircraftForward.dot(frame.tangent), -1.0f, 1.0f);
      gp_scalar upDot = std::clamp(aircraftUp.dot(frame.insideNormal), -1.0f, 1.0f);
      gp_scalar forwardAngle = std::acos(forwardDot);
      gp_scalar upAngle = std::acos(upDot);
      gp_scalar forwardPercent = (forwardAngle * static_cast<gp_scalar>(100.0f)) / static_cast<gp_scalar>(M_PI);
      gp_scalar upPercent = (upAngle * static_cast<gp_scalar>(100.0f)) / static_cast<gp_scalar>(M_PI);
      gp_scalar orientationStep = pow(forwardPercent, ORIENTATION_TANGENT_WEIGHT) +
                               pow(upPercent, ORIENTATION_UP_WEIGHT);
      gp_scalar inversionPenalty = (upDot < 0.0f) ? (-upDot) * PATH_UP_INVERSION_PENALTY : 0.0f;
      orientation_alignment_sum += orientationStep + inversionPenalty;
      orientationPenaltyForLog = orientationStep + inversionPenalty;

      // Movement efficiency: compare aircraft movement to optimal path movement
      gp_scalar movement_efficiency = 0.0f;
      if (stepIndex > 1 && pathIndex < path.size() - 1) {
        // Aircraft movement since last step
        gp_vec3 aircraft_movement = stepAircraftState.getPosition() - aircraftState.at(stepIndex-1).getPosition();
        gp_scalar aircraft_distance = aircraft_movement.norm();
        
        // Optimal path movement (direction from current waypoint to next)
        int movementNextIndex = std::min(pathIndex + 1, static_cast<int>(path.size()) - 1);
        gp_vec3 path_direction = path.at(movementNextIndex).start - path.at(pathIndex).start;
        if (path_direction.norm() < static_cast<gp_scalar>(1e-5)) {
          path_direction = frame.tangent;
        } else {
          path_direction.normalize();
        }
        gp_scalar optimal_distance = aircraft_distance; // Same distance, optimal direction
        
        if (aircraft_distance > static_cast<gp_scalar>(0.1f)) { // Avoid division by zero
          // Compare directions: dot product shows alignment (-1 to +1)
          gp_scalar direction_alignment = aircraft_movement.normalized().dot(path_direction);
          
          // Movement efficiency penalty: 0 = perfect alignment, 2 = opposite direction
          movement_efficiency = (static_cast<gp_scalar>(1.0f) - direction_alignment); // 0 to 2 range
          
          // Additional penalty for excessive turning/rolling
          gp_scalar current_roll = abs(stepAircraftState.getRollCommand());
          gp_scalar current_pitch = abs(stepAircraftState.getPitchCommand());
          
          // Control saturation penalty - heavily penalize near ±1.0 commands
          gp_scalar saturation_penalty = 0.0f;
          if (current_roll > CONTROL_SATURATION_THRESHOLD) {
            saturation_penalty += (current_roll - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
          }
          if (current_pitch > CONTROL_SATURATION_THRESHOLD) {
            saturation_penalty += (current_pitch - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
          }
          
          // Control rate penalty - penalize rapid changes (bang-bang behavior)
          gp_scalar roll_change = abs(stepAircraftState.getRollCommand() - roll_prev);
          gp_scalar pitch_change = abs(stepAircraftState.getPitchCommand() - pitch_prev);
          gp_scalar control_rate_penalty = (roll_change + pitch_change) * CONTROL_RATE_PENALTY;
          
          // Estimate "optimal" control needed for this path segment
          // Simple model: straight segments need minimal control, curves need some roll
          gp_vec3 prev_path_dir = (pathIndex > 0) ?
            (path.at(pathIndex).start - path.at(pathIndex-1).start) : path_direction;
          if (prev_path_dir.norm() < static_cast<gp_scalar>(1e-5)) {
            prev_path_dir = frame.tangent;
          } else {
            prev_path_dir.normalize();
          }
          gp_scalar path_curvature = (static_cast<gp_scalar>(1.0f) - prev_path_dir.dot(path_direction)); // 0 to 2
          gp_scalar optimal_roll = path_curvature * static_cast<gp_scalar>(0.3f); // Scale to reasonable roll estimate
          gp_scalar optimal_pitch = static_cast<gp_scalar>(0.1f); // Assume minimal pitch needed
          
          // Penalty for excessive control relative to path requirements
          gp_scalar control_excess = (current_roll - optimal_roll) + (current_pitch - optimal_pitch);
          control_excess = std::max(static_cast<gp_scalar>(0.0f), control_excess); // Only penalize excess, not deficit
          
          movement_efficiency += control_excess * 0.5; // Add control excess penalty
          movement_efficiency += saturation_penalty; // Add saturation penalty
          movement_efficiency += control_rate_penalty; // Add rate change penalty
        }
      }
      // Normalize to [100:0] scale - increased max penalty due to saturation/rate penalties  
      movement_efficiency = movement_efficiency * static_cast<gp_scalar>(100.0f) / static_cast<gp_scalar>(MOVEMENT_EFFICIENCY_MAX_PENALTY);

      // ready for next cycle
      roll_prev = stepAircraftState.getRollCommand();
      pitch_prev = stepAircraftState.getPitchCommand();
      throttle_prev = stepAircraftState.getThrottleCommand();

      // accumulate the error
      distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
      angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
      movement_efficiency_sum += pow(movement_efficiency, MOVEMENT_EFFICIENCY_WEIGHT); // Stronger penalty weight than distance/angle
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

    // tally up the normlized fitness based on steps and progress
    gp_scalar normalized_distance_error = (distance_error_sum / simulation_steps);
    gp_scalar normalized_angle_align = (angle_error_sum / simulation_steps);
    gp_scalar normalized_movement_efficiency = (movement_efficiency_sum / simulation_steps);
    gp_scalar normalized_cross_track = (cross_track_error_sum / simulation_steps);
    gp_scalar normalized_cross_track_osc = cross_track_osc_events > 0
      ? (cross_track_osc_sum / cross_track_osc_events)
      : 0.0f;
    gp_scalar normalized_orientation_penalty = (orientation_alignment_sum / simulation_steps);
    localFitness = normalized_distance_error + normalized_angle_align + normalized_movement_efficiency
      + normalized_cross_track + normalized_cross_track_osc + normalized_orientation_penalty;

    if (isnan(localFitness)) {
      nanDetector++;
    }

    if (crashReason != CrashReason::None) {
      gp_scalar fractional_distance_remaining = static_cast<gp_scalar>(1.0f) - path.at(aircraftState.back().getThisPathIndex()).distanceFromStart / path.back().distanceFromStart;
      localFitness += SIM_CRASH_PENALTY * fractional_distance_remaining;
    }

    // accumulate the local fitness
    stdFitness += localFitness;
  }

  // normalize
  stdFitness /= context.evalResults.pathList.size();
  
  // Mark fitness as valid for the GP library
  fitnessValid = 1;
  
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

  // Init GP system.
  GPInit(1, -1);

  // Initialize ConfigManager (this replaces the old GPConfiguration call)
  ConfigManager::initialize(configFile, *logger.info());

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
  *logger.info() << "WindSeedStride: " << ConfigManager::getExtraConfig().windSeedStride << endl << endl;

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
         << std::defaultfloat << std::endl;
    bout.flush();
    lastThroughputTime = now;
    lastSimRunCount = currentRuns;
  };

  // prime the paths?
  generationPaths = generateSmoothPaths(ConfigManager::getExtraConfig().generatorMethod,
                                        ConfigManager::getExtraConfig().simNumPathsPerGen,
                                        SIM_PATH_BOUNDS, SIM_PATH_BOUNDS,
                                        ConfigManager::getExtraConfig().randomPathSeedB);
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
        sendRPC(*context.socket, evalData);
        
        // Get simulation results
        context.evalResults = receiveRPC<EvalResults>(*context.socket);
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
          
          gp_scalar distance_error_sum = 0;
          gp_scalar angle_error_sum = 0;
          gp_scalar movement_efficiency_sum = 0;
          gp_scalar cross_track_error_sum = 0;
          gp_scalar cross_track_osc_sum = 0;
          int cross_track_osc_events = 0;
          gp_scalar orientation_alignment_sum = 0;
          gp_scalar prev_cross_track_signed = 0;
          bool has_prev_cross_track = false;
          int simulation_steps = 0;
          
          gp_scalar roll_prev = aircraftState.at(stepIndex).getRollCommand();
          gp_scalar pitch_prev = aircraftState.at(stepIndex).getPitchCommand();
          gp_scalar throttle_prev = aircraftState.at(stepIndex).getThrottleCommand();
          
          while (++stepIndex < aircraftState.size()) {
            auto& stepAircraftState = aircraftState.at(stepIndex);
            int rawPathIndex = stepAircraftState.getThisPathIndex();
            int pathIndex = std::clamp(rawPathIndex, 0, static_cast<int>(path.size()) - 1);
            const Path& currentPathPoint = path.at(pathIndex);
            PathFrame frame = computePathFrame(path, pathIndex);
            int nextIndex = std::min(pathIndex + 1, static_cast<int>(path.size()) - 1);
            
            gp_vec3 aircraftPosition = stepAircraftState.getPosition();
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

            gp_quat craftOrientation = stepAircraftState.getOrientation();
            gp_vec3 aircraftForward = craftOrientation * gp_vec3::UnitX();
            gp_vec3 aircraftUp = craftOrientation * gp_vec3::UnitZ();

            gp_scalar crossTrackPercentForLog = 0.0f;
            gp_scalar oscillationPercentForLog = 0.0f;
            gp_scalar orientationPenaltyForLog = 0.0f;

            gp_vec3 craftOffset = aircraftPosition - currentPathPoint.start;
            gp_vec3 lateral = craftOffset - frame.tangent * craftOffset.dot(frame.tangent);
            gp_scalar crossTrackSigned = lateral.dot(frame.insideNormal);
            gp_scalar crossTrackMagnitude = lateral.norm();
            if (SIM_PATH_RADIUS_LIMIT > static_cast<gp_scalar>(0.0f)) {
              crossTrackPercentForLog = (crossTrackMagnitude * static_cast<gp_scalar>(100.0f)) / SIM_PATH_RADIUS_LIMIT;
              cross_track_error_sum += pow(crossTrackPercentForLog, CROSS_TRACK_WEIGHT);
            }

            gp_scalar oscillationPercent = 0.0f;
            bool oscillationTriggered = false;
            if (has_prev_cross_track) {
              if ((crossTrackSigned > CROSS_TRACK_SIGN_THRESHOLD && prev_cross_track_signed < -CROSS_TRACK_SIGN_THRESHOLD) ||
                  (crossTrackSigned < -CROSS_TRACK_SIGN_THRESHOLD && prev_cross_track_signed > CROSS_TRACK_SIGN_THRESHOLD)) {
                gp_scalar denomOsc = std::max(static_cast<gp_scalar>(1.0f), static_cast<gp_scalar>(2.0f) * SIM_PATH_RADIUS_LIMIT);
                oscillationPercent = (std::abs(prev_cross_track_signed) + std::abs(crossTrackSigned)) * static_cast<gp_scalar>(100.0f) / denomOsc;
                cross_track_osc_sum += pow(oscillationPercent, CROSS_TRACK_OSC_WEIGHT);
                cross_track_osc_events++;
                oscillationTriggered = true;
              }
            }
            oscillationPercentForLog = oscillationTriggered ? oscillationPercent : 0.0;
            prev_cross_track_signed = crossTrackSigned;
            has_prev_cross_track = true;

            gp_scalar forwardDot = std::clamp(aircraftForward.dot(frame.tangent), static_cast<gp_scalar>(-1.0f), static_cast<gp_scalar>(1.0f));
            gp_scalar upDot = std::clamp(aircraftUp.dot(frame.insideNormal), static_cast<gp_scalar>(-1.0f), static_cast<gp_scalar>(1.0f));
            gp_scalar forwardAngle = static_cast<gp_scalar>(std::acos(forwardDot));
            gp_scalar upAngle = static_cast<gp_scalar>(std::acos(upDot));
            gp_scalar forwardPercent = (forwardAngle * static_cast<gp_scalar>(100.0f)) / static_cast<gp_scalar>(M_PI);
            gp_scalar upPercent = (upAngle * static_cast<gp_scalar>(100.0f)) / static_cast<gp_scalar>(M_PI);
            gp_scalar orientationStep = pow(forwardPercent, ORIENTATION_TANGENT_WEIGHT) +
                                        pow(upPercent, ORIENTATION_UP_WEIGHT);
            gp_scalar inversionPenalty = (upDot < static_cast<gp_scalar>(0.0f)) ? (-upDot) * PATH_UP_INVERSION_PENALTY : static_cast<gp_scalar>(0.0f);
            orientation_alignment_sum += orientationStep + inversionPenalty;
            orientationPenaltyForLog = orientationStep + inversionPenalty;
            
            // Movement efficiency: compare aircraft movement to optimal path movement
            gp_scalar movement_efficiency = 0.0f;
            if (stepIndex > 1 && pathIndex < path.size() - 1) {
              // Aircraft movement since last step
              gp_vec3 aircraft_movement = stepAircraftState.getPosition() - aircraftState.at(stepIndex-1).getPosition();
              gp_scalar aircraft_distance = aircraft_movement.norm();
              
              // Optimal path movement (direction from current waypoint to next)
              int movementNextIndex = std::min(pathIndex + 1, static_cast<int>(path.size()) - 1);
              gp_vec3 path_direction = path.at(movementNextIndex).start - path.at(pathIndex).start;
              if (path_direction.norm() < static_cast<gp_scalar>(1e-5f)) {
                path_direction = frame.tangent;
              } else {
                path_direction.normalize();
              }
              if (aircraft_distance > static_cast<gp_scalar>(0.1f)) { // Avoid division by zero
                // Compare directions: dot product shows alignment (-1 to +1)
                gp_scalar direction_alignment = aircraft_movement.normalized().dot(path_direction);
                
                // Movement efficiency penalty: 0 = perfect alignment, 2 = opposite direction
                movement_efficiency = (static_cast<gp_scalar>(1.0f) - direction_alignment); // 0 to 2 range
                
                // Additional penalty for excessive turning/rolling
                gp_scalar current_roll = abs(stepAircraftState.getRollCommand());
                gp_scalar current_pitch = abs(stepAircraftState.getPitchCommand());
                
                // Control saturation penalty - heavily penalize near ±1.0 commands
                gp_scalar saturation_penalty = 0.0f;
                if (current_roll > CONTROL_SATURATION_THRESHOLD) {
                  saturation_penalty += (current_roll - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
                }
                if (current_pitch > CONTROL_SATURATION_THRESHOLD) {
                  saturation_penalty += (current_pitch - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
                }
                
                // Control rate penalty - penalize rapid changes (bang-bang behavior)
                gp_scalar roll_change = abs(stepAircraftState.getRollCommand() - roll_prev);
                gp_scalar pitch_change = abs(stepAircraftState.getPitchCommand() - pitch_prev);
                gp_scalar control_rate_penalty = (roll_change + pitch_change) * CONTROL_RATE_PENALTY;
                
                // Estimate "optimal" control needed for this path segment
                // Simple model: straight segments need minimal control, curves need some roll
                gp_vec3 prev_path_dir = (pathIndex > 0) ?
                  (path.at(pathIndex).start - path.at(pathIndex-1).start) : path_direction;
                if (prev_path_dir.norm() < static_cast<gp_scalar>(1e-5f)) {
                  prev_path_dir = frame.tangent;
                } else {
                  prev_path_dir.normalize();
                }
                gp_scalar path_curvature = (static_cast<gp_scalar>(1.0f) - prev_path_dir.dot(path_direction)); // 0 to 2
                gp_scalar optimal_roll = path_curvature * static_cast<gp_scalar>(0.3f); // Scale to reasonable roll estimate
                gp_scalar optimal_pitch = static_cast<gp_scalar>(0.1f); // Assume minimal pitch needed
                
                // Penalty for excessive control relative to path requirements
                gp_scalar control_excess = (current_roll - optimal_roll) + (current_pitch - optimal_pitch);
                control_excess = std::max(static_cast<gp_scalar>(0.0f), control_excess); // Only penalize excess, not deficit
                
                movement_efficiency += control_excess * static_cast<gp_scalar>(0.5f); // Add control excess penalty
                movement_efficiency += saturation_penalty; // Add saturation penalty
                movement_efficiency += control_rate_penalty; // Add rate change penalty
              }
            }
            // Normalize to [100:0] scale - increased max penalty due to saturation/rate penalties  
            movement_efficiency = movement_efficiency * static_cast<gp_scalar>(100.0f) / static_cast<gp_scalar>(MOVEMENT_EFFICIENCY_MAX_PENALTY);
            
            roll_prev = stepAircraftState.getRollCommand();
            pitch_prev = stepAircraftState.getPitchCommand();
            throttle_prev = stepAircraftState.getThrottleCommand();
            
            distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
            angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
            movement_efficiency_sum += pow(movement_efficiency, MOVEMENT_EFFICIENCY_WEIGHT);
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
          
          gp_scalar normalized_distance_error = (distance_error_sum / simulation_steps);
          gp_scalar normalized_angle_align = (angle_error_sum / simulation_steps);
          gp_scalar normalized_movement_efficiency = (movement_efficiency_sum / simulation_steps);
          gp_scalar normalized_cross_track = (cross_track_error_sum / simulation_steps);
          gp_scalar normalized_cross_track_osc = cross_track_osc_events > 0
            ? (cross_track_osc_sum / cross_track_osc_events)
            : static_cast<gp_scalar>(0.0f);
          gp_scalar normalized_orientation_penalty = (orientation_alignment_sum / simulation_steps);
          localFitness = normalized_distance_error + normalized_angle_align + normalized_movement_efficiency
            + normalized_cross_track + normalized_cross_track_osc + normalized_orientation_penalty;
          
          if (isnan(localFitness)) {
            nanDetector++;
          }
          
          if (crashReason != CrashReason::None) {
            gp_scalar fractional_distance_remaining = static_cast<gp_scalar>(1.0f) - path.at(aircraftState.back().getThisPathIndex()).distanceFromStart / path.back().distanceFromStart;
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

      std::ostringstream oss;
      boost::archive::text_oarchive oa(oss);
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
      // For this generation, build a smooth path goal
      generationPaths = generateSmoothPaths(ConfigManager::getExtraConfig().generatorMethod,
                                            ConfigManager::getExtraConfig().simNumPathsPerGen,
                                            SIM_PATH_BOUNDS, SIM_PATH_BOUNDS,
                                            ConfigManager::getExtraConfig().randomPathSeedB);
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

      MyGP* best = pop->NthMyGP(pop->bestOfPopulation);
      best->setScenarioIndex(computeScenarioIndexForIndividual(pop->bestOfPopulation));
      best->evaluate();

      // reverse order names for s3...
      computedKeyName = startTime + "/gen" + std::to_string(10000 - gen) + ".dmp";
      pop->endOfEvaluation();

      printEval = false;

      // Create a report of this generation and how well it is doing
      if (nanDetector > 0) {
        *logger.warn() << "NanDetector count: " << nanDetector << endl;
      }
      pop->createGenerationReport(0, gen, fout, bout, *logger.info());
      logGenerationStats(gen);
    }

    // TODO send exit message to workers

    // go ahead and dump out the best of the best
    // ofstream bestGP("best.dat");
    // pop->NthMyGP(pop->bestOfPopulation)->save(bestGP);

    *logger.info() << "GP complete!" << endl;
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

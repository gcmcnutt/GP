
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
extern std::atomic_ulong nanDetector;

std::vector<std::vector<Path>> generationPaths;
std::vector<ScenarioDescriptor> generationScenarios;
std::atomic_ulong nanDetector = 0;
std::ofstream fout;
EvalResults bestOfEvalResults;
EvalResults aggregatedEvalResults;
EvalResults* activeEvalCollector = nullptr;
std::atomic<uint64_t> evaluationProgress{0};



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
    scenario.windSeed = seedBase;
    scenario.pathVariantIndex = -1;
    scenario.windVariantIndex = 0;
    generationScenarios.push_back(std::move(scenario));
    return;
  }

  for (int windIdx = 0; windIdx < windScenarioCount; ++windIdx) {
    ScenarioDescriptor scenario;
    scenario.pathList = basePaths;
    scenario.windSeed = seedBase + seedStride * static_cast<unsigned int>(windIdx);
    scenario.pathVariantIndex = -1;
    scenario.windVariantIndex = windIdx;
    generationScenarios.push_back(std::move(scenario));
  }

  if (generationScenarios.empty()) {
    ScenarioDescriptor scenario;
    scenario.pathVariantIndex = -1;
    scenario.windVariantIndex = 0;
    scenario.windSeed = seedBase;
    scenario.pathList = basePaths;
    generationScenarios.push_back(std::move(scenario));
  }
}

int computeScenarioIndexForIndividual(int individualIndex) {
  const GPVariables& gpCfg = ConfigManager::getGPConfig();
  int scenarioCount = std::max<int>(generationScenarios.size(), 1);

  if (gpCfg.DemeticGrouping && gpCfg.DemeSize > 0) {
    int demeIndex = individualIndex / gpCfg.DemeSize;
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
      MyGP* best = NthMyGP(bestOfPopulation);

      if (!generationScenarios.empty()) {
        int originalScenarioIndex = best->getScenarioIndex();
        clearEvalResults(aggregatedEvalResults);
        aggregatedEvalResults.scenario.pathVariantIndex = -1;
        aggregatedEvalResults.scenario.windVariantIndex = -1;
        aggregatedEvalResults.scenario.windSeed = 0;

        EvalResults* previousCollector = activeEvalCollector;
        activeEvalCollector = &aggregatedEvalResults;

        for (size_t scenarioIdx = 0; scenarioIdx < generationScenarios.size(); ++scenarioIdx) {
          best->setScenarioIndex(static_cast<int>(scenarioIdx));
          threadPool->enqueue([best](WorkerContext& context) {
            best->evalTask(context);
          });
          threadPool->wait_for_tasks();
        }

        best->setScenarioIndex(originalScenarioIndex);
        activeEvalCollector = previousCollector;
        bestOfEvalResults = aggregatedEvalResults;
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

  EvalData evalData;
  evalData.gp = buffer;
  evalData.pathList = scenario.pathList;
  evalData.scenario.pathVariantIndex = scenario.pathVariantIndex;
  evalData.scenario.windVariantIndex = scenario.windVariantIndex;
  evalData.scenario.windSeed = scenario.windSeed;

  sendRPC(*context.socket, evalData);

  // How did it go?
  context.evalResults = receiveRPC<EvalResults>(*context.socket);
  if (context.evalResults.scenario.windSeed == 0 && evalData.scenario.windSeed != 0) {
    context.evalResults.scenario = evalData.scenario;
  }
  if (context.evalResults.pathList.empty() && !evalData.pathList.empty()) {
    context.evalResults.pathList = evalData.pathList;
  }
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

    // compute this path fitness
    double localFitness = 0;
    int stepIndex = 0; // where are we on the flight path?

    // error accumulators
    double distance_error_sum = 0;
    double angle_error_sum = 0;
    double movement_efficiency_sum = 0;
    int simulation_steps = 0;

    // initial states
    double roll_prev = aircraftState.at(stepIndex).getRollCommand();
    double pitch_prev = aircraftState.at(stepIndex).getPitchCommand();
    double throttle_prev = aircraftState.at(stepIndex).getThrottleCommand();

    // now walk next steps of actual path
    while (++stepIndex < aircraftState.size()) {
      auto& stepAircraftState = aircraftState.at(stepIndex);
      int pathIndex = stepAircraftState.getThisPathIndex();

      // Compute the distance between the aircraft and the goal
      double distanceFromGoal = (path.at(pathIndex).start - stepAircraftState.getPosition()).norm();
      // normalize [100:0]
      distanceFromGoal = distanceFromGoal * 100.0 / (2 * SIM_PATH_RADIUS_LIMIT);

      // Compute vector from me to target
      Eigen::Vector3d target_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start);
      Eigen::Vector3d aircraft_to_target = (path.at(pathIndex).start - stepAircraftState.getPosition());
      double dot_product = target_direction.dot(aircraft_to_target);
      double angle_rad = std::acos(std::clamp(dot_product / (target_direction.norm() * aircraft_to_target.norm()), -1.0, 1.0));
      // normalize [100:0]
      angle_rad = angle_rad * 100.0 / M_PI;

      // Movement efficiency: compare aircraft movement to optimal path movement
      double movement_efficiency = 0.0;
      if (stepIndex > 1 && pathIndex < path.size() - 1) {
        // Aircraft movement since last step
        Eigen::Vector3d aircraft_movement = stepAircraftState.getPosition() - aircraftState.at(stepIndex-1).getPosition();
        double aircraft_distance = aircraft_movement.norm();
        
        // Optimal path movement (direction from current waypoint to next)
        Eigen::Vector3d path_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start).normalized();
        double optimal_distance = aircraft_distance; // Same distance, optimal direction
        
        if (aircraft_distance > 0.1) { // Avoid division by zero
          // Compare directions: dot product shows alignment (-1 to +1)
          double direction_alignment = aircraft_movement.normalized().dot(path_direction);
          
          // Movement efficiency penalty: 0 = perfect alignment, 2 = opposite direction
          movement_efficiency = (1.0 - direction_alignment); // 0 to 2 range
          
          // Additional penalty for excessive turning/rolling
          double current_roll = abs(stepAircraftState.getRollCommand());
          double current_pitch = abs(stepAircraftState.getPitchCommand());
          
          // Control saturation penalty - heavily penalize near ±1.0 commands
          double saturation_penalty = 0.0;
          if (current_roll > CONTROL_SATURATION_THRESHOLD) {
            saturation_penalty += (current_roll - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
          }
          if (current_pitch > CONTROL_SATURATION_THRESHOLD) {
            saturation_penalty += (current_pitch - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
          }
          
          // Control rate penalty - penalize rapid changes (bang-bang behavior)
          double roll_change = abs(stepAircraftState.getRollCommand() - roll_prev);
          double pitch_change = abs(stepAircraftState.getPitchCommand() - pitch_prev);
          double control_rate_penalty = (roll_change + pitch_change) * CONTROL_RATE_PENALTY;
          
          // Estimate "optimal" control needed for this path segment
          // Simple model: straight segments need minimal control, curves need some roll
          Eigen::Vector3d prev_path_dir = (pathIndex > 0) ? 
            (path.at(pathIndex).start - path.at(pathIndex-1).start).normalized() : path_direction;
          double path_curvature = (1.0 - prev_path_dir.dot(path_direction)); // 0 to 2
          double optimal_roll = path_curvature * 0.3; // Scale to reasonable roll estimate
          double optimal_pitch = 0.1; // Assume minimal pitch needed
          
          // Penalty for excessive control relative to path requirements
          double control_excess = (current_roll - optimal_roll) + (current_pitch - optimal_pitch);
          control_excess = std::max(0.0, control_excess); // Only penalize excess, not deficit
          
          movement_efficiency += control_excess * 0.5; // Add control excess penalty
          movement_efficiency += saturation_penalty; // Add saturation penalty
          movement_efficiency += control_rate_penalty; // Add rate change penalty
        }
      }
      // Normalize to [100:0] scale - increased max penalty due to saturation/rate penalties  
      movement_efficiency = movement_efficiency * 100.0 / MOVEMENT_EFFICIENCY_MAX_PENALTY;

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
          fout << "Pth:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP  movEffP\n";
          printHeader = false;
        }

        // convert aircraft_orientaton to euler
        Eigen::Matrix3d rotMatrix = stepAircraftState.getOrientation().toRotationMatrix();

        // Extract Euler angles
        // Note: atan2 returns angle in range [-pi, pi]
        Eigen::Vector3d euler;

        // Handle special case near pitch = ±90° (gimbal lock)
        if (std::abs(rotMatrix(2, 0)) > 0.99999) {
          // Gimbal lock case
          euler[0] = 0; // Roll becomes undefined, set to zero

          // Determine pitch based on r31 sign
          if (rotMatrix(2, 0) > 0) {
            euler[1] = -M_PI / 2; // pitch = -90°
            euler[2] = -atan2(rotMatrix(1, 2), rotMatrix(0, 2)); // yaw
          }
          else {
            euler[1] = M_PI / 2;  // pitch = 90°
            euler[2] = atan2(rotMatrix(1, 2), rotMatrix(0, 2)); // yaw
          }
        }
        else {
          // Normal case
          euler[0] = atan2(rotMatrix(2, 1), rotMatrix(2, 2)); // roll (phi)
          euler[1] = -asin(rotMatrix(2, 0));                 // pitch (theta)
          euler[2] = atan2(rotMatrix(1, 0), rotMatrix(0, 0)); // yaw (psi)
        }

        char outbuf[1000]; // XXX use c++20
        sprintf(outbuf, "%03d:%04d: %06ld %3d % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n",
          i, simulation_steps,
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
          movement_efficiency
        );
        fout << outbuf;
      }
    }

    // tally up the normlized fitness based on steps and progress
    double normalized_distance_error = (distance_error_sum / simulation_steps);
    double normalized_angle_align = (angle_error_sum / simulation_steps);
    double normalized_movement_efficiency = (movement_efficiency_sum / simulation_steps);
    localFitness = normalized_distance_error + normalized_angle_align + normalized_movement_efficiency;

    if (isnan(localFitness)) {
      nanDetector++;
    }

    if (crashReason != CrashReason::None) {
      double fractional_distance_remaining = 1.0 - path.at(aircraftState.back().getThisPathIndex()).distanceFromStart / path.back().distanceFromStart;
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

  std::string startTime = generate_iso8601_timestamp();

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
  ofstream bout(strStatFile.str());

  // prime the paths?
  generationPaths = generateSmoothPaths(ConfigManager::getExtraConfig().generatorMethod, ConfigManager::getExtraConfig().simNumPathsPerGen, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);
  rebuildGenerationScenarios(generationPaths);
  *logger.debug() << "Wind scenarios this generation: paths="
                 << ConfigManager::getExtraConfig().simNumPathsPerGen
                 << " windScenarios=" << ConfigManager::getExtraConfig().windScenarioCount
                 << " total=" << generationScenarios.size() << endl;
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
        evalData.scenario.pathVariantIndex = scenario.pathVariantIndex;
        evalData.scenario.windVariantIndex = scenario.windVariantIndex;
        evalData.scenario.windSeed = scenario.windSeed;
        sendRPC(*context.socket, evalData);
        
        // Get simulation results
        context.evalResults = receiveRPC<EvalResults>(*context.socket);
        if (context.evalResults.pathList.empty() && !evalData.pathList.empty()) {
          context.evalResults.pathList = evalData.pathList;
        }
        if (context.evalResults.scenario.windSeed == 0 && evalData.scenario.windSeed != 0) {
          context.evalResults.scenario = evalData.scenario;
        }

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
          
          double localFitness = 0;
          int stepIndex = 0;
          
          double distance_error_sum = 0;
          double angle_error_sum = 0;
          double movement_efficiency_sum = 0;
          int simulation_steps = 0;
          
          double roll_prev = aircraftState.at(stepIndex).getRollCommand();
          double pitch_prev = aircraftState.at(stepIndex).getPitchCommand();
          double throttle_prev = aircraftState.at(stepIndex).getThrottleCommand();
          
          while (++stepIndex < aircraftState.size()) {
            auto& stepAircraftState = aircraftState.at(stepIndex);
            int pathIndex = stepAircraftState.getThisPathIndex();
            
            double distanceFromGoal = (path.at(pathIndex).start - stepAircraftState.getPosition()).norm();
            distanceFromGoal = distanceFromGoal * 100.0 / (2 * SIM_PATH_RADIUS_LIMIT);
            
            Eigen::Vector3d target_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start);
            Eigen::Vector3d aircraft_to_target = (path.at(pathIndex).start - stepAircraftState.getPosition());
            double dot_product = target_direction.dot(aircraft_to_target);
            double angle_rad = std::acos(std::clamp(dot_product / (target_direction.norm() * aircraft_to_target.norm()), -1.0, 1.0));
            angle_rad = angle_rad * 100.0 / M_PI;
            
            // Movement efficiency: compare aircraft movement to optimal path movement
            double movement_efficiency = 0.0;
            if (stepIndex > 1 && pathIndex < path.size() - 1) {
              // Aircraft movement since last step
              Eigen::Vector3d aircraft_movement = stepAircraftState.getPosition() - aircraftState.at(stepIndex-1).getPosition();
              double aircraft_distance = aircraft_movement.norm();
              
              // Optimal path movement (direction from current waypoint to next)
              Eigen::Vector3d path_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start).normalized();
              double optimal_distance = aircraft_distance; // Same distance, optimal direction
              
              if (aircraft_distance > 0.1) { // Avoid division by zero
                // Compare directions: dot product shows alignment (-1 to +1)
                double direction_alignment = aircraft_movement.normalized().dot(path_direction);
                
                // Movement efficiency penalty: 0 = perfect alignment, 2 = opposite direction
                movement_efficiency = (1.0 - direction_alignment); // 0 to 2 range
                
                // Additional penalty for excessive turning/rolling
                double current_roll = abs(stepAircraftState.getRollCommand());
                double current_pitch = abs(stepAircraftState.getPitchCommand());
                
                // Control saturation penalty - heavily penalize near ±1.0 commands
                double saturation_penalty = 0.0;
                if (current_roll > CONTROL_SATURATION_THRESHOLD) {
                  saturation_penalty += (current_roll - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
                }
                if (current_pitch > CONTROL_SATURATION_THRESHOLD) {
                  saturation_penalty += (current_pitch - CONTROL_SATURATION_THRESHOLD) * CONTROL_SATURATION_PENALTY;
                }
                
                // Control rate penalty - penalize rapid changes (bang-bang behavior)
                double roll_change = abs(stepAircraftState.getRollCommand() - roll_prev);
                double pitch_change = abs(stepAircraftState.getPitchCommand() - pitch_prev);
                double control_rate_penalty = (roll_change + pitch_change) * CONTROL_RATE_PENALTY;
                
                // Estimate "optimal" control needed for this path segment
                // Simple model: straight segments need minimal control, curves need some roll
                Eigen::Vector3d prev_path_dir = (pathIndex > 0) ? 
                  (path.at(pathIndex).start - path.at(pathIndex-1).start).normalized() : path_direction;
                double path_curvature = (1.0 - prev_path_dir.dot(path_direction)); // 0 to 2
                double optimal_roll = path_curvature * 0.3; // Scale to reasonable roll estimate
                double optimal_pitch = 0.1; // Assume minimal pitch needed
                
                // Penalty for excessive control relative to path requirements
                double control_excess = (current_roll - optimal_roll) + (current_pitch - optimal_pitch);
                control_excess = std::max(0.0, control_excess); // Only penalize excess, not deficit
                
                movement_efficiency += control_excess * 0.5; // Add control excess penalty
                movement_efficiency += saturation_penalty; // Add saturation penalty
                movement_efficiency += control_rate_penalty; // Add rate change penalty
              }
            }
            // Normalize to [100:0] scale - increased max penalty due to saturation/rate penalties  
            movement_efficiency = movement_efficiency * 100.0 / MOVEMENT_EFFICIENCY_MAX_PENALTY;
            
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
                fout << "Pth:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP  movEffP\n";
              }
              
              Eigen::Matrix3d rotMatrix = stepAircraftState.getOrientation().toRotationMatrix();
              Eigen::Vector3d euler;
              
              if (std::abs(rotMatrix(2, 0)) > 0.99999) {
                euler[0] = 0;
                if (rotMatrix(2, 0) > 0) {
                  euler[1] = -M_PI / 2;
                  euler[2] = -atan2(rotMatrix(1, 2), rotMatrix(0, 2));
                } else {
                  euler[1] = M_PI / 2;
                  euler[2] = atan2(rotMatrix(1, 2), rotMatrix(0, 2));
                }
              } else {
                euler[0] = atan2(rotMatrix(2, 1), rotMatrix(2, 2));
                euler[1] = -asin(rotMatrix(2, 0));
                euler[2] = atan2(rotMatrix(1, 0), rotMatrix(0, 0));
              }
              
              char outbuf[1000];
              sprintf(outbuf, "%03d:%04d: %06ld %3d % 8.2f% 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f %8.2f %8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f % 8.2f\n",
                i, simulation_steps,
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
                movement_efficiency
              );
              fout << outbuf;
            }
          }
          
          double normalized_distance_error = (distance_error_sum / simulation_steps);
          double normalized_angle_align = (angle_error_sum / simulation_steps);
          double normalized_movement_efficiency = (movement_efficiency_sum / simulation_steps);
          localFitness = normalized_distance_error + normalized_angle_align + normalized_movement_efficiency;
          
          if (isnan(localFitness)) {
            nanDetector++;
          }
          
          if (crashReason != CrashReason::None) {
            double fractional_distance_remaining = 1.0 - path.at(aircraftState.back().getThisPathIndex()).distanceFromStart / path.back().distanceFromStart;
            localFitness += SIM_CRASH_PENALTY * fractional_distance_remaining;
          }
          
          stdFitness += localFitness;
        }
        
        stdFitness /= context.evalResults.pathList.size();
        
        // Mark fitness as valid for the GP library
        fitnessValid = 1;
      }
      
      double getFinalFitness() const { return stdFitness; }
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
    double cumulativeFitness = 0.0;
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
    
    double averageFitness = scenariosEvaluated > 0 ? cumulativeFitness / static_cast<double>(scenariosEvaluated) : 0.0;
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

    // This next for statement is the actual genetic programming system
    // which is in essence just repeated reproduction and crossover loop
    // through all the generations ...
    MyPopulation* newPop = NULL;

    for (int gen = 1; gen <= ConfigManager::getGPConfig().NumberOfGenerations; gen++)
    {
      // For this generation, build a smooth path goal
      generationPaths = generateSmoothPaths(ConfigManager::getExtraConfig().generatorMethod, ConfigManager::getExtraConfig().simNumPathsPerGen, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);
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
    }

    // TODO send exit message to workers

    // go ahead and dump out the best of the best
    // ofstream bestGP("best.dat");
    // pop->NthMyGP(pop->bestOfPopulation)->save(bestGP);

    *logger.info() << "GP complete!" << endl;
  }
}

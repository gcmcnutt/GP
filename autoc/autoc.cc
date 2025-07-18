
// autoc.cc

/* -------------------------------------------------------------------
From skeleton/skeleton.cc
------------------------------------------------------------------- */

#include <boost/asio.hpp>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <new>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <memory>

#include "gp.h"
#include "gp_bytecode.h"
#include "gpconfig.h"
#include "minisim.h"
#include "threadpool.h"
#include "autoc.h"
#include "logger.h"
#include "pathgen.h"
#include "gp_evaluator.h"

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

// Special GP class for evaluation mode that uses generated code
class EvaluationGP : public MyGP
{
public:
  EvaluationGP() : MyGP(1) {}
  
  virtual void evaluate() override {
    // Use the generated evaluateGP function instead of tree evaluation
    tasks.push_back(this);
  }
  
  virtual void evalTask(WorkerContext& context) override {
    stdFitness = 0;
    
    // Send evaluation request to simulator using generated code
    std::vector<char> buffer;
    boost::iostreams::stream<boost::iostreams::back_insert_device<std::vector<char>>> outStream(buffer);
    save(outStream);
    outStream.flush();
    EvalData evalData = { buffer, generationPaths };
    sendRPC(*context.socket, evalData);
    
    // Get simulation results
    context.evalResults = receiveRPC<EvalResults>(*context.socket);
    
    // Instead of normal GP evaluation, use the generated evaluateGP function
    // to compute control commands and measure fitness
    for (int i = 0; i < context.evalResults.pathList.size(); i++) {
      auto& path = context.evalResults.pathList.at(i);
      auto& aircraftState = context.evalResults.aircraftStateList.at(i);
      auto& crashReason = context.evalResults.crashReasonList.at(i);
      
      double localFitness = 0;
      int stepIndex = 0;
      
      double distance_error_sum = 0;
      double angle_error_sum = 0;
      double control_smoothness_sum = 0;
      int simulation_steps = 0;
      
      double roll_prev = aircraftState.at(stepIndex).getRollCommand();
      double pitch_prev = aircraftState.at(stepIndex).getPitchCommand();
      double throttle_prev = aircraftState.at(stepIndex).getThrottleCommand();
      
      while (++stepIndex < aircraftState.size()) {
        auto& stepAircraftState = aircraftState.at(stepIndex);
        int pathIndex = stepAircraftState.getThisPathIndex();
        
        // Apply generated GP control logic
        double gpResult = evaluateGP(const_cast<AircraftState&>(stepAircraftState), path, 0.0);
        
        double distanceFromGoal = (path.at(pathIndex).start - stepAircraftState.getPosition()).norm();
        distanceFromGoal = distanceFromGoal * 100.0 / (2 * SIM_PATH_RADIUS_LIMIT);
        
        Eigen::Vector3d target_direction = (path.at(pathIndex + 1).start - path.at(pathIndex).start);
        Eigen::Vector3d aircraft_to_target = (path.at(pathIndex).start - stepAircraftState.getPosition());
        double dot_product = target_direction.dot(aircraft_to_target);
        double angle_rad = std::acos(std::clamp(dot_product / (target_direction.norm() * aircraft_to_target.norm()), -1.0, 1.0));
        angle_rad = angle_rad * 100.0 / M_PI;
        
        double smoothness = pow(roll_prev - stepAircraftState.getRollCommand(), 2.0);
        smoothness += pow(pitch_prev - stepAircraftState.getPitchCommand(), 2.0);
        smoothness += pow(throttle_prev - stepAircraftState.getThrottleCommand(), 2.0);
        smoothness = sqrt(smoothness);
        smoothness = smoothness * 100.0 / 3.46;
        
        roll_prev = stepAircraftState.getRollCommand();
        pitch_prev = stepAircraftState.getPitchCommand();
        throttle_prev = stepAircraftState.getThrottleCommand();
        
        distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
        angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
        control_smoothness_sum += pow(smoothness, FITNESS_CONTROL_WEIGHT);
        simulation_steps++;
        
        if (printEval) {
          bestOfEvalResults = context.evalResults;
          
          if (stepIndex == 1) {
            fout << "Pth:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP controlP\n";
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
            smoothness
          );
          fout << outbuf;
        }
      }
      
      double normalized_distance_error = (distance_error_sum / simulation_steps);
      double normalized_angle_align = (angle_error_sum / simulation_steps);
      double normalized_control_smoothness = (control_smoothness_sum / simulation_steps);
      localFitness = normalized_distance_error + normalized_angle_align + normalized_control_smoothness;
      
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
  }
};
std::vector<std::vector<Path>> generationPaths;
std::atomic_ulong nanDetector = 0;
std::ofstream fout;
EvalResults bestOfEvalResults;

// Define configuration parameters and the neccessary array to
// read/write the configuration to a file.  If you need more
// variables, just add them below and insert an entry in the
// configArray.
GPVariables cfg;
ExtraConfig extraCfg;
struct GPConfigVarInformation configArray[] =
{
  {"PopulationSize", DATAINT, &cfg.PopulationSize},
  {"NumberOfGenerations", DATAINT, &cfg.NumberOfGenerations},
  {"CreationType", DATAINT, &cfg.CreationType},
  {"CrossoverProbability", DATADOUBLE, &cfg.CrossoverProbability},
  {"CreationProbability", DATADOUBLE, &cfg.CreationProbability},
  {"MaximumDepthForCreation", DATAINT, &cfg.MaximumDepthForCreation},
  {"MaximumDepthForCrossover", DATAINT, &cfg.MaximumDepthForCrossover},
  {"SelectionType", DATAINT, &cfg.SelectionType},
  {"TournamentSize", DATAINT, &cfg.TournamentSize},
  {"DemeticGrouping", DATAINT, &cfg.DemeticGrouping},
  {"DemeSize", DATAINT, &cfg.DemeSize},
  {"DemeticMigProbability", DATADOUBLE, &cfg.DemeticMigProbability},
  {"SwapMutationProbability", DATADOUBLE, &cfg.SwapMutationProbability},
  {"ShrinkMutationProbability", DATADOUBLE, &cfg.ShrinkMutationProbability},
  {"AddBestToNewPopulation", DATAINT, &cfg.AddBestToNewPopulation},
  {"SteadyState", DATAINT, &cfg.SteadyState},
  {"SimNumPathsPerGeneration", DATAINT, &extraCfg.simNumPathsPerGen},
  {"EvalThreads", DATAINT, &extraCfg.evalThreads},
  {"PathGeneratorMethod", DATASTRING, &extraCfg.generatorMethod},
  {"MinisimProgram", DATASTRING, &extraCfg.minisimProgram},
  {"MinisimPortOverride", DATAINT, &extraCfg.minisimPortOverride},
  {"S3Bucket", DATASTRING, &extraCfg.s3Bucket},
  {"S3Profile", DATASTRING, &extraCfg.s3Profile},
  {"EvaluateMode", DATAINT, &extraCfg.evaluateMode},
  {"BytecodeFile", DATASTRING, &extraCfg.bytecodeFile},
  {"", DATAINT, NULL}
};


ThreadPool* threadPool;
std::string computedKeyName;

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
  if (strcmp("default", extraCfg.s3Profile) != 0) {

    Aws::Client::ClientConfiguration clientConfig;
    clientConfig.endpointOverride = "http://localhost:9000"; // MinIO server address
    clientConfig.scheme = Aws::Http::Scheme::HTTP; // Use HTTP instead of HTTPS
    clientConfig.verifySSL = false; // Disable SSL verification for local testing

    auto policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Never;
    auto credentialsProvider = Aws::MakeShared<Aws::Auth::ProfileConfigFileAWSCredentialsProvider>("CredentialsProvider", extraCfg.s3Profile);
    auto s3_client = Aws::MakeShared<Aws::S3::S3Client>("S3Client",
      credentialsProvider,
      clientConfig,
      Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Never,
      false
    );
    return s3_client;
  }
  else {
    return Aws::MakeShared<Aws::S3::S3Client>("S3Client");
  }
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

  virtual void endOfEvaluation() {

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
    if (printEval) {
      // now put the resulting elements into the S3 object
      Aws::S3::Model::PutObjectRequest request;
      request.SetBucket(extraCfg.s3Bucket);

      // path name is $base/RunDate/gen$gen.dmp
      request.SetKey(computedKeyName);

      std::ostringstream oss;
      boost::archive::text_oarchive oa(oss);
      oa << bestOfEvalResults;

      // TODO: dump out fitness
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
  EvalData evalData = { buffer, generationPaths };
  sendRPC(*context.socket, evalData);

  // How did it go?
  context.evalResults = receiveRPC<EvalResults>(*context.socket);
  // context.evalResults.dump(std::cerr);

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
    double control_smoothness_sum = 0;
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

      // control smoothness -- internal vector distance
      double smoothness = pow(roll_prev - stepAircraftState.getRollCommand(), 2.0);
      smoothness += pow(pitch_prev - stepAircraftState.getPitchCommand(), 2.0);
      smoothness += pow(throttle_prev - stepAircraftState.getThrottleCommand(), 2.0);
      smoothness = sqrt(smoothness);
      // normalize [100:0]
      smoothness = smoothness * 100.0 / 3.46;

      // ready for next cycle
      roll_prev = stepAircraftState.getRollCommand();
      pitch_prev = stepAircraftState.getPitchCommand();
      throttle_prev = stepAircraftState.getThrottleCommand();

      // accumulate the error
      distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
      angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
      control_smoothness_sum += pow(smoothness, FITNESS_CONTROL_WEIGHT);
      simulation_steps++;

      // use the ugly global to communicate best of gen
      if (printEval) {

        // TODO: need reference to best task somehow
        bestOfEvalResults = context.evalResults;

        if (printHeader) {
          fout << "Pth:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP controlP\n";
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
          smoothness
        );
        fout << outbuf;
      }
    }

    // tally up the normlized fitness based on steps and progress
    double normalized_distance_error = (distance_error_sum / simulation_steps);
    double normalized_angle_align = (angle_error_sum / simulation_steps);
    double normalized_control_smoothness = (control_smoothness_sum / simulation_steps);
    localFitness = normalized_distance_error + normalized_angle_align + normalized_control_smoothness;

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
}

void newHandler()
{
  cerr << "\nFatal error: Out of memory." << endl;
  exit(1);
}


int main()
{
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

  // Read configuration file.
  GPConfiguration config(*logger.info(), "autoc.ini", configArray);

  // AWS setup
  Aws::SDKOptions options;
  Aws::InitAPI(options);

  // initialize workers
  threadPool = new ThreadPool(extraCfg);

  // Print the configuration
  *logger.info() << cfg << endl;
  *logger.info() << "SimNumPathsPerGen: " << extraCfg.simNumPathsPerGen << endl;
  *logger.info() << "EvalThreads: " << extraCfg.evalThreads << endl;
  *logger.info() << "MinisimProgram: " << extraCfg.minisimProgram << endl;
  *logger.info() << "MinisimPortOverride: " << extraCfg.minisimPortOverride << endl;
  *logger.info() << "EvaluateMode: " << extraCfg.evaluateMode << endl << endl;

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
  generationPaths = generateSmoothPaths(extraCfg.generatorMethod, extraCfg.simNumPathsPerGen, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);

  if (extraCfg.evaluateMode) {
    // Evaluation mode: use bytecode interpreter with external simulators
    *logger.info() << "Running in bytecode evaluation mode with external simulators..." << endl;
    *logger.info() << "Loading bytecode file: " << extraCfg.bytecodeFile << endl;
    
    // Load the bytecode program
    GPBytecodeInterpreter interpreter;
    if (!interpreter.loadProgram(extraCfg.bytecodeFile)) {
      *logger.error() << "Failed to load bytecode file: " << extraCfg.bytecodeFile << endl;
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
        EvalData evalData = { bytecodeBuffer, generationPaths };
        sendRPC(*context.socket, evalData);
        
        // Get simulation results
        context.evalResults = receiveRPC<EvalResults>(*context.socket);
        
        // Use same fitness computation as normal GP evaluation
        for (int i = 0; i < context.evalResults.pathList.size(); i++) {
          auto& path = context.evalResults.pathList.at(i);
          auto& aircraftState = context.evalResults.aircraftStateList.at(i);
          auto& crashReason = context.evalResults.crashReasonList.at(i);
          
          double localFitness = 0;
          int stepIndex = 0;
          
          double distance_error_sum = 0;
          double angle_error_sum = 0;
          double control_smoothness_sum = 0;
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
            
            double smoothness = pow(roll_prev - stepAircraftState.getRollCommand(), 2.0);
            smoothness += pow(pitch_prev - stepAircraftState.getPitchCommand(), 2.0);
            smoothness += pow(throttle_prev - stepAircraftState.getThrottleCommand(), 2.0);
            smoothness = sqrt(smoothness);
            smoothness = smoothness * 100.0 / 3.46;
            
            roll_prev = stepAircraftState.getRollCommand();
            pitch_prev = stepAircraftState.getPitchCommand();
            throttle_prev = stepAircraftState.getThrottleCommand();
            
            distance_error_sum += pow(distanceFromGoal, FITNESS_DISTANCE_WEIGHT);
            angle_error_sum += pow(angle_rad, FITNESS_ALIGNMENT_WEIGHT);
            control_smoothness_sum += pow(smoothness, FITNESS_CONTROL_WEIGHT);
            simulation_steps++;
            
            if (printEval) {
              bestOfEvalResults = context.evalResults;
              
              if (stepIndex == 1) {
                fout << "Pth:Step:   Time Idx  totDist   pathX    pathY    pathZ        X        Y        Z       dr       dp       dy   relVel     roll    pitch    power    distP   angleP controlP\n";
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
                smoothness
              );
              fout << outbuf;
            }
          }
          
          double normalized_distance_error = (distance_error_sum / simulation_steps);
          double normalized_angle_align = (angle_error_sum / simulation_steps);
          double normalized_control_smoothness = (control_smoothness_sum / simulation_steps);
          localFitness = normalized_distance_error + normalized_angle_align + normalized_control_smoothness;
          
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
    
    *logger.info() << "Running bytecode evaluation on " << generationPaths.size() << " paths..." << endl;
    
    // Initialize evaluation
    printEval = true; // Enable detailed output
    
    // Create bytecode evaluation GP
    BytecodeEvaluationGP evalGP(bytecodeBuffer);
    evalGP.evaluate();
    
    // Execute the evaluation using the thread pool (same as MyPopulation::endOfEvaluation)
    for (auto& task : tasks) {
      threadPool->enqueue([task](WorkerContext& context) {
        task->evalTask(context);
        });
    }
    
    // wait for all tasks to finish
    threadPool->wait_for_tasks();
    tasks.clear();
    
    *logger.info() << "Bytecode evaluation complete!" << endl;
    *logger.info() << "Computed fitness: " << evalGP.getFinalFitness() << endl;
    *logger.info() << "Original GP fitness: " << interpreter.getFitness() << endl;
    
    // Set up S3 storage for results (use gen-1.dmp so renderer can interpret it)
    computedKeyName = startTime + "/gen-1.dmp";
    
    // Store results to S3
    if (bestOfEvalResults.pathList.size() > 0) {
      Aws::S3::Model::PutObjectRequest request;
      request.SetBucket(extraCfg.s3Bucket);
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
    MyPopulation* pop = new MyPopulation(cfg, adfNs);
    pop->create();
    *logger.info() << "Ok." << endl;
    pop->createGenerationReport(1, 0, fout, bout, *logger.info());

    // This next for statement is the actual genetic programming system
    // which is in essence just repeated reproduction and crossover loop
    // through all the generations ...
    MyPopulation* newPop = NULL;

    for (int gen = 1; gen <= cfg.NumberOfGenerations; gen++)
    {
      // For this generation, build a smooth path goal
      generationPaths = generateSmoothPaths(extraCfg.generatorMethod, extraCfg.simNumPathsPerGen, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);

      // Create a new generation from the old one by applying the genetic operators
      if (!cfg.SteadyState)
        newPop = new MyPopulation(cfg, adfNs);
      pop->generate(*newPop);

      // TODO fix this pattern to use a dynamic logger
      printEval = true;
      MyGP* best = pop->NthMyGP(pop->bestOfPopulation);
      best->evaluate();

      // reverse order names for s3...
      computedKeyName = startTime + "/gen" + std::to_string(10000 - gen) + ".dmp";
      pop->endOfEvaluation();

      printEval = false;

      // Delete the old generation and make the new the old one
      if (!cfg.SteadyState)
      {
        delete pop;
        pop = newPop;
      }

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

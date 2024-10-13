
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

#include "gp.h"
#include "gpconfig.h"
#include "minisim.h"
#include "threadpool.h"
#include "autoc.h"
#include "logger.h"
#include "pathgen.h"

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/PutObjectRequest.h>

#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/support/date_time.hpp>

using namespace std;
namespace logging = boost::log;

Logger logger;

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
  {"MinisimProgram", DATASTRING, &extraCfg.minisimProgram},
  {"MinisimPortOverride", DATAINT, &extraCfg.minisimPortOverride},
  {"S3Bucket", DATASTRING, &extraCfg.s3Bucket},
  {"", DATAINT, NULL}
};


ThreadPool* threadPool;
std::ofstream fout;
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


class MyPopulation : public GPPopulation
{
public:
  // Constructor (mandatory)
  MyPopulation(GPVariables& GPVar_, GPAdfNodeSet& adfNs_) :
    GPPopulation(GPVar_, adfNs_) {}

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
    Aws::S3::S3Client s3_client;

    // dispatch all the GPs now (TODO this may still work inline with evaluate)
    for (auto& task : tasks) {
      threadPool->enqueue([task](WorkerContext& context) {
        task->evalTask(context);
        });
    }

    // wait for all tasks to finish
    threadPool->wait_for_tasks();
    tasks.clear();

    if (printEval) {
      // now put the resulting elements into the S3 object
      Aws::S3::Model::PutObjectRequest request;
      request.SetBucket(extraCfg.s3Bucket);

      // path name is $base/RunDate/gen$gen.dmp
      request.SetKey(computedKeyName);

      std::ostringstream oss;
      boost::archive::text_oarchive oa(oss);
      oa << evalResults;

      // TODO: dump the selected GP
      // TODO: dump out fitness
      std::shared_ptr<Aws::StringStream> ss = Aws::MakeShared<Aws::StringStream>("");
      *ss << oss.str();
      request.SetBody(ss);

      auto outcome = s3_client.PutObject(request);
      if (!outcome.IsSuccess()) {
        *logger.error() << "Error: " << outcome.GetError().GetMessage() << std::endl;
      }

      // clear out elements for next pass
      evalResults.crashReasonList.clear();
      evalResults.pathList.clear();
      evalResults.aircraftStateList.clear();
    }
  }

  // Print (not mandatory) 
  // virtual void printOn (ostream& os);

  // Access genetic programs (not mandatory)
  MyGP* NthMyGP(int n) {
    return (MyGP*)GPContainer::Nth(n);
  }
};



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
  *logger.info() << "MinisimPortOverride: " << extraCfg.minisimPortOverride << endl << endl;

  // Create the adf function/terminal set and print it out.
  GPAdfNodeSet adfNs;
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
  generationPaths = generateSmoothPaths(extraCfg.simNumPathsPerGen, NUM_SEGMENTS_PER_PATH, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);

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
    generationPaths = generateSmoothPaths(extraCfg.simNumPathsPerGen, NUM_SEGMENTS_PER_PATH, SIM_PATH_BOUNDS, SIM_PATH_BOUNDS);

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

    if (printEval) {
      best->printAircraftState(fout, evalResults);

      // TODO // now update points for Renderer
      // planPath.push_back(path.at(pathIndex).start); // XXX push the full path
      // actualPath.push_back(stepAircraftState.getPosition());

      // // negative here as world up is negative Z
      // actualOrientation.push_back(stepAircraftState.getOrientation() * -Eigen::Vector3d::UnitZ());
    }

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

  // go ahead and dump out the best of the best
  // ofstream bestGP("best.dat");
  // pop->NthMyGP(pop->bestOfPopulation)->save(bestGP);

  *logger.info() << "GP complete!" << endl;
}

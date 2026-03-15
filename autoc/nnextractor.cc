#include <regex>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iomanip>

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/s3/model/ListObjectsV2Request.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/client/ClientConfiguration.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "../include/gp.h"
#include "minisim.h"
#include "nn_serialization.h"
#include "config_manager.h"
#include "autoc.h"
#include "threadpool.h"

// Dummy implementations required by linker (real ones are in autoc.cc)
void MyGP::evaluate() {}
void MyGP::evalTask(WorkerContext& context) {}

void printUsage(const char* progName) {
  std::cout << "Usage: " << progName << " [OPTIONS]\n";
  std::cout << "Options:\n";
  std::cout << "  -k, --keyname KEYNAME    Specify run key name (autoc-timestamp)\n";
  std::cout << "  -g, --generation GEN     Specify generation number (default: latest)\n";
  std::cout << "  -o, --output FILE        Output NN weight file (default: nn_weights.dat)\n";
  std::cout << "  -i, --config FILE        Use specified config file (default: autoc.ini)\n";
  std::cout << "  -h, --help               Show this help message\n";
  std::cout << "\n";
  std::cout << "Extracts the best NN genome from an S3 evolution archive.\n";
  std::cout << "The archive contains Boost-serialized EvalResults with the NN genome\n";
  std::cout << "embedded in the .gp field (NN01 binary format).\n";
  std::cout << "\n";
  std::cout << "Examples:\n";
  std::cout << "  " << progName << "                          # Extract latest run, latest generation\n";
  std::cout << "  " << progName << " -k autoc-20250101-12345  # Extract specific run\n";
  std::cout << "  " << progName << " -g 50                    # Extract generation 50\n";
  std::cout << "  " << progName << " -o my_weights.dat        # Custom output file\n";
}

// Extract the generation number from the S3 key
int extractGenNumber(const std::string& input) {
  std::regex pattern("autoc-.*\\/gen(\\d+)\\.dmp");
  std::smatch matches;

  if (std::regex_search(input, matches, pattern) && matches.size() > 1) {
    return std::stoi(matches[1].str());
  }

  return -1;
}

int main(int argc, char** argv) {
  static struct option long_options[] = {
    {"keyname", required_argument, 0, 'k'},
    {"generation", required_argument, 0, 'g'},
    {"output", required_argument, 0, 'o'},
    {"config", required_argument, 0, 'i'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

  std::string computedKeyName = "";
  int specifiedGeneration = -1;
  std::string outputFile = "nn_weights.dat";
  std::string configFile = "autoc.ini";
  int option_index = 0;
  int c;

  while ((c = getopt_long(argc, argv, "k:g:o:i:h", long_options, &option_index)) != -1) {
    switch (c) {
      case 'k': computedKeyName = optarg; break;
      case 'g': specifiedGeneration = std::stoi(optarg); break;
      case 'o': outputFile = optarg; break;
      case 'i': configFile = optarg; break;
      case 'h': printUsage(argv[0]); return 0;
      case '?': printUsage(argv[0]); return 1;
      default: break;
    }
  }

  // Positional argument for backward compat
  if (computedKeyName.empty() && optind < argc) {
    computedKeyName = argv[optind];
  }

  // Initialize configuration and AWS
  ConfigManager::initialize(configFile);

  Aws::SDKOptions options;
  Aws::InitAPI(options);

  auto s3_client = ConfigManager::getS3Client();
  std::string bucket = ConfigManager::getConfig().s3Bucket;

  // Find latest run if not specified
  if (computedKeyName.empty()) {
    Aws::S3::Model::ListObjectsV2Request listFolders;
    listFolders.SetBucket(bucket);
    listFolders.SetPrefix("autoc-");
    listFolders.SetDelimiter("/");

    bool isTruncated = false;
    do {
      auto outcome = s3_client->ListObjectsV2(listFolders);
      if (!outcome.IsSuccess()) {
        std::cerr << "Error listing S3 objects: " << outcome.GetError().GetMessage() << std::endl;
        return 1;
      }
      const auto& result = outcome.GetResult();
      for (const auto& commonPrefix : result.GetCommonPrefixes()) {
        computedKeyName = commonPrefix.GetPrefix();
        break;
      }
      isTruncated = result.GetIsTruncated();
      if (isTruncated) {
        listFolders.SetContinuationToken(result.GetNextContinuationToken());
      }
    } while (isTruncated);
  }

  if (computedKeyName.empty()) {
    std::cerr << "No autoc runs found in S3 bucket: " << bucket << std::endl;
    return 1;
  }

  // Find target generation
  std::string keyName;
  if (specifiedGeneration >= 0) {
    keyName = computedKeyName + "gen" + std::to_string(specifiedGeneration) + ".dmp";
  } else {
    Aws::S3::Model::ListObjectsV2Request listItem;
    listItem.SetBucket(bucket);
    listItem.SetPrefix(computedKeyName + "gen");
    bool isTruncated = false;
    do {
      auto outcome = s3_client->ListObjectsV2(listItem);
      if (!outcome.IsSuccess()) {
        std::cerr << "Error listing generations: " << outcome.GetError().GetMessage() << std::endl;
        return 1;
      }
      for (const auto& object : outcome.GetResult().GetContents()) {
        keyName = object.GetKey();
        break;
      }
      isTruncated = outcome.GetResult().GetIsTruncated();
      if (isTruncated) {
        listItem.SetContinuationToken(outcome.GetResult().GetNextContinuationToken());
      }
    } while (isTruncated);
  }

  if (keyName.empty()) {
    std::cerr << "No generation data found for run: " << computedKeyName << std::endl;
    return 1;
  }

  std::cout << "Fetching: " << keyName << " from " << bucket << std::endl;

  // Fetch from S3
  Aws::S3::Model::GetObjectRequest request;
  request.SetBucket(bucket);
  request.SetKey(keyName);
  auto outcome = s3_client->GetObject(request);

  if (!outcome.IsSuccess()) {
    std::cerr << "Error retrieving " << keyName << ": " << outcome.GetError().GetMessage() << std::endl;
    return 1;
  }

  // Deserialize EvalResults
  EvalResults evalResults;
  {
    std::ostringstream oss;
    oss << outcome.GetResult().GetBody().rdbuf();
    std::string retrievedData = oss.str();

    try {
      std::istringstream iss(retrievedData, std::ios::binary);
      boost::archive::binary_iarchive ia(iss);
      ia >> evalResults;
    }
    catch (const std::exception& e) {
      std::cerr << "Error deserializing EvalResults: " << e.what() << std::endl;
      return 1;
    }
  }

  // Extract NN genome from .gp field
  if (evalResults.gp.empty()) {
    std::cerr << "No controller data found in EvalResults" << std::endl;
    return 1;
  }

  // Verify it's NN format
  if (!nn_detect_format(reinterpret_cast<const uint8_t*>(evalResults.gp.data()),
                        evalResults.gp.size())) {
    std::cerr << "Error: S3 archive does not contain NN data (no NN01 magic bytes)." << std::endl;
    std::cerr << "This may be a GP archive. Use gpextractor instead." << std::endl;
    return 1;
  }

  // Deserialize to get metadata for display
  NNGenome genome;
  if (!nn_deserialize(reinterpret_cast<const uint8_t*>(evalResults.gp.data()),
                      evalResults.gp.size(), genome)) {
    std::cerr << "Error deserializing NN genome from archive" << std::endl;
    return 1;
  }

  // Write the raw NN01 binary directly (already serialized)
  std::ofstream file(outputFile, std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot create output file: " << outputFile << std::endl;
    return 1;
  }
  file.write(evalResults.gp.data(), evalResults.gp.size());
  file.close();

  // Print summary
  int generation = extractGenNumber(keyName);
  std::cout << "\nExtracted NN genome to: " << outputFile << std::endl;
  std::cout << "  S3 Key:    " << keyName << std::endl;
  std::cout << "  Generation: " << generation << std::endl;
  std::cout << "  Topology:   ";
  for (size_t i = 0; i < genome.topology.size(); i++) {
    if (i > 0) std::cout << " -> ";
    std::cout << genome.topology[i];
  }
  std::cout << std::endl;
  std::cout << "  Weights:    " << genome.weights.size() << std::endl;
  std::cout << "  Fitness:    " << std::fixed << std::setprecision(6) << genome.fitness << std::endl;
  std::cout << "  Sigma:      " << std::fixed << std::setprecision(6) << genome.mutation_sigma << std::endl;
  std::cout << "  File size:  " << evalResults.gp.size() << " bytes" << std::endl;

  Aws::ShutdownAPI(options);
  return 0;
}

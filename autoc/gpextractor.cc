#include <regex>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

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
#include "gp_bytecode.h"
#include "config_manager.h"
#include "autoc.h"
#include "threadpool.h"

// Node set definitions and createNodeSet function are now in autoc-eval.cc

// MyGene and MyGP classes are now in autoc.h

// Provide dummy implementations for gpextractor (real implementations are in autoc.cc)
void MyGP::evaluate() {
  // Dummy implementation for gpextractor - not used for loading/extraction
}

void MyGP::evalTask(WorkerContext& context) {
  // Dummy implementation for gpextractor - not used for loading/extraction
}

// Utility function from minisim.cc  
boost::iostreams::stream<boost::iostreams::array_source> charArrayToIstream(const std::vector<char>& charArray) {
  return boost::iostreams::stream<boost::iostreams::array_source>(
    boost::iostreams::array_source(charArray.data(), charArray.size())
  );
}

std::shared_ptr<Aws::S3::S3Client> getS3Client() {
  return ConfigManager::getS3Client();
}

void printUsage(const char* progName) {
  std::cout << "Usage: " << progName << " [OPTIONS]\n";
  std::cout << "Options:\n";
  std::cout << "  -k, --keyname KEYNAME    Specify GP log key name\n";
  std::cout << "  -g, --generation GEN     Specify generation number\n";
  std::cout << "  -b, --bytecode           Generate bytecode for interpreter\n";
  std::cout << "  -o, --output FILE        Output file for bytecode (default: gp_program.dat)\n";
  std::cout << "  -i, --config FILE        Use specified config file (default: autoc.ini)\n";
  std::cout << "  -h, --help               Show this help message\n";
  std::cout << "\n";
  std::cout << "Examples:\n";
  std::cout << "  " << progName << "                          # Extract latest run, latest generation\n";
  std::cout << "  " << progName << " -k autoc-20250101-12345  # Extract specific run, latest generation\n";
  std::cout << "  " << progName << " -g 50                    # Extract latest run, generation 50\n";
  std::cout << "  " << progName << " -b -o gp_program.dat     # Generate bytecode file\n";
}

// Extract the generation number from the key
int extractGenNumber(const std::string& input) {
  std::regex pattern("autoc-.*\\/gen(\\d+)\\.dmp");
  std::smatch matches;

  if (std::regex_search(input, matches, pattern) && matches.size() > 1) {
    return std::stoi(matches[1].str());
  }

  return -1;
}

// Generate bytecode from GP tree preserving left-to-right evaluation order
// This is critical for operations with side effects like SETROLL, SETPITCH, etc.
void generateBytecode(MyGene* gene, std::vector<GPBytecode>& program) {
  if (!gene) {
    return;
  }
  
  int nodeValue = gene->geneNode().value();
  if (nodeValue == -1) {
    return;
  }
  
  // Handle different node types to preserve evaluation order
  switch (nodeValue) {
    // Operations that need left-to-right evaluation (like GP tree traversal)
    case ADD:
    case SUB:
    case MUL:
    case DIV:
    case EQ:
    case GT:
      // Generate left child first, then right child, then operation
      if (gene->containerSize() >= 2) {
        generateBytecode(gene->NthMyChild(0), program);
        generateBytecode(gene->NthMyChild(1), program);
      }
      break;
      
    case IF:
      // IF needs condition, then true branch, then false branch
      if (gene->containerSize() >= 3) {
        generateBytecode(gene->NthMyChild(0), program); // condition
        generateBytecode(gene->NthMyChild(1), program); // true branch
        generateBytecode(gene->NthMyChild(2), program); // false branch
      }
      break;
      
    case PROGN:
      // PROGN executes children in order, returns last
      for (int i = 0; i < gene->containerSize(); i++) {
        generateBytecode(gene->NthMyChild(i), program);
      }
      break;
      
    // Single-argument operations (unary)
    case SIN:
    case COS:
    case SETPITCH:
    case SETROLL:
    case SETTHROTTLE:
    case ABS:
    case SQRT:
      // Generate child first
      if (gene->containerSize() >= 1) {
        generateBytecode(gene->NthMyChild(0), program);
      }
      break;
      
    // Two-argument operations (binary) - evaluate left-to-right
    case ATAN2:
    case MIN:
    case MAX:
      if (gene->containerSize() >= 2) {
        generateBytecode(gene->NthMyChild(0), program);
        generateBytecode(gene->NthMyChild(1), program);
      }
      break;
      
    // Three-argument operations (ternary)
    case CLAMP:
      if (gene->containerSize() >= 3) {
        generateBytecode(gene->NthMyChild(0), program); // value
        generateBytecode(gene->NthMyChild(1), program); // min
        generateBytecode(gene->NthMyChild(2), program); // max
      }
      break;
      
    // Sensor functions with optional arguments
    case GETDPHI:
    case GETDTHETA:
    case GETDTARGET:
      if (gene->containerSize() >= 1) {
        generateBytecode(gene->NthMyChild(0), program);
      }
      break;
      
    // Terminals and constants - no children to process
    case OP_PI:
    case ZERO:
    case ONE:
    case TWO:
    case GETPITCH:
    case GETROLL:
    case GETTHROTTLE:
    case GETVEL:
    case GETDHOME:
    case GETALPHA:
    case GETBETA:
    case GETVELX:
    case GETVELY:
    case GETVELZ:
    case GETROLL_RAD:
    case GETPITCH_RAD:
      // No children to process
      break;
      
    default:
      // For unknown operations, fall back to processing all children
      for (int i = 0; i < gene->containerSize(); i++) {
        generateBytecode(gene->NthMyChild(i), program);
      }
      break;
  }
  
  // Add the current node instruction
  GPBytecode instruction;
  instruction.opcode = nodeValue;
  instruction.argc = gene->containerSize();
  instruction.constant = 0.0f;
  
  // For constants, store the value
  switch (nodeValue) {
    case OP_PI:
      instruction.constant = M_PI;
      break;
    case ZERO:
      instruction.constant = 0.0f;
      break;
    case ONE:
      instruction.constant = 1.0f;
      break;
    case TWO:
      instruction.constant = 2.0f;
      break;
  }
  
  program.push_back(instruction);
}

// Export GP as bytecode file
bool exportBytecode(MyGP& gp, const std::string& outputFile, const std::string& s3Key, int generation) {
  std::ofstream file(outputFile, std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot create bytecode file: " << outputFile << std::endl;
    return false;
  }
  
  // Generate bytecode program
  std::vector<GPBytecode> program;
  MyGene* rootGene = (MyGene*)gp.NthGene(0);
  generateBytecode(rootGene, program);
  
  // Create header
  GPBytecodeHeader header;
  header.magic = GPBytecodeHeader::MAGIC;
  header.version = GPBytecodeHeader::VERSION;
  header.instruction_count = program.size();
  header.fitness_int = (uint32_t)(gp.getFitness() * 1000000.0); // Fixed point
  header.length = gp.length();
  header.depth = gp.depth();
  header.generation = generation;
  
  // Copy S3 key (truncate if too long)
  strncpy(header.s3_key, s3Key.c_str(), sizeof(header.s3_key) - 1);
  header.s3_key[sizeof(header.s3_key) - 1] = '\0';
  
  // Write header
  file.write(reinterpret_cast<const char*>(&header), sizeof(header));
  if (file.fail()) {
    std::cerr << "Error writing bytecode header" << std::endl;
    return false;
  }
  
  // Write instructions
  for (const auto& instruction : program) {
    file.write(reinterpret_cast<const char*>(&instruction), sizeof(instruction));
    if (file.fail()) {
      std::cerr << "Error writing bytecode instruction" << std::endl;
      return false;
    }
  }
  
  file.close();
  
  std::cout << "Generated bytecode file: " << outputFile << std::endl;
  std::cout << "  Instructions: " << program.size() << std::endl;
  std::cout << "  Original GP Length: " << gp.length() << ", Depth: " << gp.depth() << std::endl;
  std::cout << "  Fitness: " << gp.getFitness() << std::endl;
  
  return true;
}

int main(int argc, char** argv) {
  // Parse command line arguments
  static struct option long_options[] = {
    {"keyname", required_argument, 0, 'k'},
    {"generation", required_argument, 0, 'g'},
    {"bytecode", no_argument, 0, 'b'},
    {"output", required_argument, 0, 'o'},
    {"config", required_argument, 0, 'i'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  std::string computedKeyName = "";
  int specifiedGeneration = -1;
  bool generateBytecode = false;
  std::string outputFile = "gp_program.dat";
  std::string configFile = "autoc.ini";
  int option_index = 0;
  int c;
  
  while ((c = getopt_long(argc, argv, "k:g:bo:i:h", long_options, &option_index)) != -1) {
    switch (c) {
      case 'k':
        computedKeyName = optarg;
        break;
      case 'g':
        specifiedGeneration = std::stoi(optarg);
        break;
      case 'b':
        generateBytecode = true;
        break;
      case 'o':
        outputFile = optarg;
        break;
      case 'i':
        configFile = optarg;
        break;
      case 'h':
        printUsage(argv[0]);
        return 0;
      case '?':
        printUsage(argv[0]);
        return 1;
      default:
        break;
    }
  }
  
  // Handle positional arguments for backward compatibility
  if (computedKeyName.empty() && optind < argc) {
    computedKeyName = argv[optind];
  }
  
  std::string keyName = "";
  
  // Initialize configuration
  ConfigManager::initialize(configFile);
  
  // AWS setup
  Aws::SDKOptions options;
  Aws::InitAPI(options);

  std::shared_ptr<Aws::S3::S3Client> s3_client = getS3Client();
  
  // should we look up the latest run?
  if (computedKeyName.empty()) {
    Aws::S3::Model::ListObjectsV2Request listFolders;
    listFolders.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
    listFolders.SetPrefix("autoc-");
    listFolders.SetDelimiter("/");

    bool isTruncated = false;
    do {
      auto outcome = s3_client->ListObjectsV2(listFolders);
      if (outcome.IsSuccess()) {
        const auto& result = outcome.GetResult();

        // Process common prefixes (these are our 'folders')
        for (const auto& commonPrefix : result.GetCommonPrefixes()) {
          computedKeyName = commonPrefix.GetPrefix();
          break;
        }

        // Check if there are more results to retrieve
        isTruncated = result.GetIsTruncated();
        if (isTruncated) {
          listFolders.SetContinuationToken(result.GetNextContinuationToken());
        }
      }
      else {
        std::cerr << "Error listing objects: " << outcome.GetError().GetMessage() << std::endl;
        return 1;
      }
    } while (isTruncated);
  }

  // ok for this run, look for the last generation (or specific generation)
  if (specifiedGeneration >= 0) {
    // Use specified generation
    keyName = computedKeyName + "gen" + std::to_string(specifiedGeneration) + ".dmp";
  } else {
    // Find the latest generation
    Aws::S3::Model::ListObjectsV2Request listItem;
    listItem.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
    listItem.SetPrefix(computedKeyName + "gen");
    bool isTruncated = false;
    do {
      auto outcome = s3_client->ListObjectsV2(listItem);
      if (outcome.IsSuccess()) {
        // Objects are already in reverse lexicographical order
        for (const auto& object : outcome.GetResult().GetContents()) {
          keyName = object.GetKey();
          break;
        }

        // Check if the response is truncated
        isTruncated = outcome.GetResult().GetIsTruncated();
        if (isTruncated) {
          // Set the continuation token for the next request
          listItem.SetContinuationToken(outcome.GetResult().GetNextContinuationToken());
        }
      }
      else {
        std::cerr << "Error listing objects: " << outcome.GetError().GetMessage() << std::endl;
        return 1;
      }
    } while (isTruncated);
  }

  if (keyName.empty()) {
    std::cerr << "No generation data found for run: " << computedKeyName << std::endl;
    return 1;
  }

  // Fetch the object from S3
  EvalResults evalResults;
  Aws::S3::Model::GetObjectRequest request;
  request.SetBucket(ConfigManager::getExtraConfig().s3Bucket);
  request.SetKey(keyName);
  auto outcome = s3_client->GetObject(request);
  
  if (!outcome.IsSuccess()) {
    std::cerr << "Error retrieving object " << keyName << " from S3: " << outcome.GetError().GetMessage() << std::endl;
    return 1;
  }

  // Deserialize the data
  std::ostringstream oss;
  oss << outcome.GetResult().GetBody().rdbuf();
  std::string retrievedData = oss.str();

  try {
    std::istringstream iss(retrievedData, std::ios::binary);
    boost::archive::binary_iarchive ia(iss);
    ia >> evalResults;
  }
  catch (const std::exception& e) {
    std::cerr << "Error during deserialization: " << e.what() << std::endl;
    return 1;
  }

  // Initialize GP system
  GPInit(0, 1);  // Don't print copyright, use seed 1

  // Register kernel classes for loading 
  GPRegisterKernelClasses();
  
  // Register our custom classes
  GPRegisterClass(new MyGene());
  GPRegisterClass(new MyGP());

  // Set up node sets (required for loading GPs)
  GPAdfNodeSet adfNs;
  createNodeSet(adfNs);

  // Create a MyGP object and load the serialized data
  MyGP gp;
  if (evalResults.gp.empty()) {
    std::cerr << "No GP data found in results" << std::endl;
    return 1;
  }

  // Use the same loading method as minisim.cc
  boost::iostreams::stream<boost::iostreams::array_source> gpStream = charArrayToIstream(evalResults.gp);
  
  char* loadResult = gp.load(gpStream);
  if (loadResult != nullptr) {
    std::cerr << "Error loading GP: " << loadResult << std::endl;
    return 1;
  }
  
  // Resolve node values (required for printOn to work)
  gp.resolveNodeValues(adfNs);

  if (generateBytecode) {
    // Generate bytecode file
    int generation = extractGenNumber(keyName);
    exportBytecode(gp, outputFile, keyName, generation);
  } else {
    // Print the GP in human readable format
    std::cout << "Extracted GP from: " << keyName << std::endl;
    std::cout << "Generation: " << extractGenNumber(keyName) << std::endl;
    std::cout << "GP Length: " << gp.length() << std::endl;
    std::cout << "GP Depth: " << gp.depth() << std::endl;
    std::cout << "Fitness: " << gp.getFitness() << std::endl;
    std::cout << std::endl;
    
    // Print the GP tree structure
    gp.printOn(std::cout);
  }

  // Cleanup AWS
  Aws::ShutdownAPI(options);

  return 0;
}
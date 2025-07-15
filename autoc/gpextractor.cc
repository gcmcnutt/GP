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

#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "../include/gp.h"
#include "minisim.h"

// Node set definitions (from autoc.h)
enum Operators {
  ADD = 0, SUB, MUL, DIV,
  IF, EQ, GT,
  SIN, COS,
  GETDPHI, GETDTHETA, GETDTARGET, GETDHOME, GETVEL,
  GETPITCH, GETROLL, GETTHROTTLE,
  SETPITCH, SETROLL, SETTHROTTLE,
  PI, ZERO, ONE, TWO, PROGN, _END
};
const int OPERATORS_NR_ITEM = _END;

// Simple createNodeSet function for GP loading
void createNodeSet(GPAdfNodeSet& adfNs)
{
  // Reserve space for the node sets
  adfNs.reserveSpace(1);

  // Now define the function and terminal set for each ADF and place
  // function/terminal sets into overall ADF container
  GPNodeSet& ns = *new GPNodeSet(OPERATORS_NR_ITEM);

  adfNs.put(0, ns);

  // Define functions/terminals and place them into the appropriate
  // sets.  Terminals take two arguments, functions three (the third
  // parameter is the number of arguments the function has)
  ns.putNode(*new GPNode(ADD, "ADD", 2));
  ns.putNode(*new GPNode(SUB, "SUB", 2));
  ns.putNode(*new GPNode(MUL, "MUL", 2));
  ns.putNode(*new GPNode(DIV, "DIV", 2));
  ns.putNode(*new GPNode(IF, "IF", 3));
  ns.putNode(*new GPNode(EQ, "EQ", 2));
  ns.putNode(*new GPNode(GT, "GT", 2));
  ns.putNode(*new GPNode(SETPITCH, "SETPITCH", 1));
  ns.putNode(*new GPNode(SETROLL, "SETROLL", 1));
  ns.putNode(*new GPNode(SETTHROTTLE, "SETTHROTTLE", 1));
  ns.putNode(*new GPNode(GETPITCH, "GETPITCH"));
  ns.putNode(*new GPNode(GETROLL, "GETROLL"));
  ns.putNode(*new GPNode(GETTHROTTLE, "GETTHROTTLE"));
  ns.putNode(*new GPNode(SIN, "SIN", 1));
  ns.putNode(*new GPNode(COS, "COS", 1));
  ns.putNode(*new GPNode(PI, "PI"));
  ns.putNode(*new GPNode(ZERO, "0"));
  ns.putNode(*new GPNode(ONE, "1"));
  ns.putNode(*new GPNode(TWO, "2"));
  ns.putNode(*new GPNode(PROGN, "PROGN", 2));
  ns.putNode(*new GPNode(GETDPHI, "GETDPHI", 1));
  ns.putNode(*new GPNode(GETDTHETA, "GETDTHETA", 1));
  ns.putNode(*new GPNode(GETDTARGET, "GETDTARGET", 1));
  ns.putNode(*new GPNode(GETVEL, "GETVEL"));
  ns.putNode(*new GPNode(GETDHOME, "GETDHOME"));
}

// Minimal MyGene class for loading (from autoc.h)
const int MyGeneID = GPUserID;
const int MyGPID = GPUserID + 1;

class MyGene : public GPGene
{
public:
  MyGene(const MyGene& gpo) : GPGene(gpo) { }
  virtual GPObject& duplicate() { return *(new MyGene(*this)); }
  MyGene(GPNode& gpo) : GPGene(gpo) {}
  MyGene() {}
  virtual int isA() { return MyGeneID; }
  virtual GPObject* createObject() { return new MyGene; }
  virtual GPGene* createChild(GPNode& gpo) {
    return new MyGene(gpo);
  }
  
  // Helper method for code generation
  int getNodeValue() const { return node ? node->value() : -1; }
  MyGene* getNthChild(int n) { return (MyGene*)NthChild(n); }
};

// Minimal MyGP class for loading (from autoc.h)
class MyGP : public GP
{
public:
  MyGP(MyGP& gpo) : GP(gpo) { }
  virtual GPObject& duplicate() { return *(new MyGP(*this)); }
  MyGP(int genes) : GP(genes) {}
  MyGP() {}
  virtual GPGene* createGene(GPNode& gpo) {
    return new MyGene(gpo);
  }
  virtual void evaluate() {} // Dummy implementation
  virtual int isA() { return MyGPID; }
  virtual GPObject* createObject() { return new MyGP; }
};

// Utility function from minisim.cc  
boost::iostreams::stream<boost::iostreams::array_source> charArrayToIstream(const std::vector<char>& charArray) {
  return boost::iostreams::stream<boost::iostreams::array_source>(
    boost::iostreams::array_source(charArray.data(), charArray.size())
  );
}

std::shared_ptr<Aws::S3::S3Client> getS3Client() {
  // real S3 or local minio?
  Aws::Client::ClientConfiguration clientConfig;
  Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::RequestDependent;
  if (strcmp("default", "minio" /*extraCfg.s3Profile*/) != 0) {
    clientConfig.endpointOverride = "http://localhost:9000"; // MinIO server address
    clientConfig.scheme = Aws::Http::Scheme::HTTP; // Use HTTP instead of HTTPS
    clientConfig.verifySSL = false; // Disable SSL verification for local testing

    policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Never;
  }

  auto credentialsProvider = Aws::MakeShared<Aws::Auth::ProfileConfigFileAWSCredentialsProvider>("CredentialsProvider", "minio" /*extraCfg.s3Profile*/);

  return Aws::MakeShared<Aws::S3::S3Client>("S3Client",
    credentialsProvider,
    clientConfig,
    Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Always,
    false
  );
}

void printUsage(const char* progName) {
  std::cout << "Usage: " << progName << " [OPTIONS]\n";
  std::cout << "Options:\n";
  std::cout << "  -k, --keyname KEYNAME    Specify GP log key name\n";
  std::cout << "  -g, --generation GEN     Specify generation number\n";
  std::cout << "  -c, --codegen            Generate C++ evaluator code\n";
  std::cout << "  -o, --output FILE        Output file for generated code\n";
  std::cout << "  -h, --help               Show this help message\n";
  std::cout << "\n";
  std::cout << "Examples:\n";
  std::cout << "  " << progName << "                          # Extract latest run, latest generation\n";
  std::cout << "  " << progName << " -k autoc-20250101-12345  # Extract specific run, latest generation\n";
  std::cout << "  " << progName << " -g 50                    # Extract latest run, generation 50\n";
  std::cout << "  " << progName << " -c -o gp_evaluator.cpp   # Generate C++ evaluator code\n";
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

// Generate C++ code for a GP tree node
std::string generateNodeCode(MyGene* gene, int& tempVarCounter) {
  if (!gene) {
    return "0.0";
  }
  
  int nodeValue = gene->getNodeValue();
  if (nodeValue == -1) {
    return "0.0";
  }
  std::string result;
  
  switch (nodeValue) {
    case ADD: {
      std::string left = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string right = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      result = "(" + left + " + " + right + ")";
      break;
    }
    case SUB: {
      std::string left = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string right = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      result = "(-" + left + " - " + right + ")";
      break;
    }
    case MUL: {
      std::string left = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string right = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      result = "(" + left + " * " + right + ")";
      break;
    }
    case DIV: {
      std::string dividend = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string divisor = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      std::string tempVar = "temp" + std::to_string(tempVarCounter++);
      result = "([&]() { double " + tempVar + " = " + divisor + "; return (" + tempVar + " == 0.0) ? 0.0 : (" + dividend + " / " + tempVar + "); })()";
      break;
    }
    case IF: {
      std::string condition = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string trueVal = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      std::string falseVal = generateNodeCode(gene->getNthChild(2), tempVarCounter);
      result = "((" + condition + ") ? (" + trueVal + ") : (" + falseVal + "))";
      break;
    }
    case EQ: {
      std::string left = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string right = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      result = "((" + left + ") == (" + right + ") ? 1.0 : 0.0)";
      break;
    }
    case GT: {
      std::string left = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string right = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      result = "((" + left + ") > (" + right + ") ? 1.0 : 0.0)";
      break;
    }
    case SETPITCH: {
      std::string value = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      result = "aircraftState.setPitchCommand(" + value + ")";
      break;
    }
    case SETROLL: {
      std::string value = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      result = "aircraftState.setRollCommand(" + value + ")";
      break;
    }
    case SETTHROTTLE: {
      std::string value = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      result = "aircraftState.setThrottleCommand(" + value + ")";
      break;
    }
    case GETPITCH:
      result = "aircraftState.getPitchCommand()";
      break;
    case GETROLL:
      result = "aircraftState.getRollCommand()";
      break;
    case GETTHROTTLE:
      result = "aircraftState.getThrottleCommand()";
      break;
    case SIN: {
      std::string value = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      result = "sin(" + value + ")";
      break;
    }
    case COS: {
      std::string value = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      result = "cos(" + value + ")";
      break;
    }
    case PI:
      result = "M_PI";
      break;
    case ZERO:
      result = "0.0";
      break;
    case ONE:
      result = "1.0";
      break;
    case TWO:
      result = "2.0";
      break;
    case PROGN: {
      std::string first = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string second = generateNodeCode(gene->getNthChild(1), tempVarCounter);
      result = "([&]() { " + first + "; return " + second + "; })()";
      break;
    }
    case GETVEL:
      result = "aircraftState.getRelVel()";
      break;
    case GETDPHI: {
      std::string argValue = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string tempVar = "temp" + std::to_string(tempVarCounter++);
      result = "([&]() { "
               "int " + tempVar + " = getIndex(path, " + argValue + "); "
               "Eigen::Vector3d craftToTarget = path.at(" + tempVar + ").start - aircraftState.getPosition(); "
               "Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget; "
               "Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z()); "
               "return std::atan2(projectedVector.y(), -projectedVector.z()); "
               "})()";
      break;
    }
    case GETDTHETA: {
      std::string argValue = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string tempVar = "temp" + std::to_string(tempVarCounter++);
      result = "([&]() { "
               "int " + tempVar + " = getIndex(path, " + argValue + "); "
               "Eigen::Vector3d craftToTarget = path.at(" + tempVar + ").start - aircraftState.getPosition(); "
               "Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget; "
               "Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z()); "
               "double rollEstimate = std::atan2(projectedVector.y(), -projectedVector.z()); "
               "Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX())); "
               "Eigen::Quaterniond virtualOrientation = aircraftState.getOrientation() * rollRotation; "
               "Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget; "
               "return std::atan2(-newLocalTargetVector.z(), newLocalTargetVector.x()); "
               "})()";
      break;
    }
    case GETDTARGET: {
      std::string argValue = generateNodeCode(gene->getNthChild(0), tempVarCounter);
      std::string tempVar = "temp" + std::to_string(tempVarCounter++);
      result = "([&]() { "
               "int " + tempVar + " = getIndex(path, " + argValue + "); "
               "double distance = (path.at(" + tempVar + ").start - aircraftState.getPosition()).norm(); "
               "return CLAMP_DEF((distance - 10.0) / aircraftState.getRelVel(), -1.0, 1.0); "
               "})()";
      break;
    }
    case GETDHOME:
      result = "(Eigen::Vector3d(0, 0, -10.0) - aircraftState.getPosition()).norm()";
      break;
    default:
      result = "0.0";
      break;
  }
  
  return result;
}

// Generate complete C++ evaluator file
void generateCppEvaluator(MyGP& gp, const std::string& outputFile, const std::string& s3Key, int generation) {
  std::ofstream out(outputFile);
  if (!out.is_open()) {
    std::cerr << "Error: Cannot open output file " << outputFile << std::endl;
    return;
  }
  
  out << "// Generated GP Evaluator - DO NOT EDIT MANUALLY\n";
  out << "// S3 Key: " << s3Key << "\n";
  out << "// Generation: " << generation << "\n";
  out << "// Generated from GP with fitness: " << gp.getFitness() << "\n";
  out << "// GP Length: " << gp.length() << ", Depth: " << gp.depth() << "\n";
  out << "//\n";
  out << "// Human readable GP structure:\n";
  
  // Capture printOn output in a string and add as comments
  std::ostringstream gpStructure;
  gp.printOn(gpStructure);
  std::string gpStr = gpStructure.str();
  
  // Add each line of the GP structure as a comment
  std::istringstream iss(gpStr);
  std::string line;
  while (std::getline(iss, line)) {
    out << "// " << line << "\n";
  }
  out << "\n";
  
  out << "#include <cmath>\n";
  out << "#include <algorithm>\n";
  out << "#include <vector>\n";
  out << "#ifndef GP_BUILD\n";
  out << "#include <ArduinoEigenDense.h>\n";
  out << "#else\n";
  out << "#include <Eigen/Dense>\n";
  out << "#include <Eigen/Geometry>\n";
  out << "#endif\n";
  out << "#include \"GP/autoc/aircraft_state.h\"\n\n";
  
  // Add getIndex helper function
  out << "static int getIndex(const std::vector<Path>& path, double arg) {\n";
  out << "  extern AircraftState aircraftState;\n";
  out << "  if (std::isnan(arg)) {\n";
  out << "    return aircraftState.getThisPathIndex();\n";
  out << "  }\n";
  out << "  int steps = CLAMP_DEF((int)arg, -5, 5);\n";
  out << "  double distanceSoFar = path.at(aircraftState.getThisPathIndex()).distanceFromStart;\n";
  out << "  double distanceGoal = distanceSoFar + steps * 22.0 * (200.0 / 1000.0);\n";
  out << "  int currentStep = aircraftState.getThisPathIndex();\n";
  out << "  if (steps > 0) {\n";
  out << "    while (path.at(currentStep).distanceFromStart < distanceGoal && currentStep < path.size() - 1) {\n";
  out << "      currentStep++;\n";
  out << "    }\n";
  out << "  } else {\n";
  out << "    while (path.at(currentStep).distanceFromStart > distanceGoal && currentStep > 0) {\n";
  out << "      currentStep--;\n";
  out << "    }\n";
  out << "  }\n";
  out << "  return currentStep;\n";
  out << "}\n\n";
  
  out << "double evaluateGP(AircraftState& aircraftState, const std::vector<Path>& path, double arg) {\n";
  out << "  double returnValue;\n\n";
  
  // Generate the main evaluation code
  MyGene* rootGene = (MyGene*)gp.NthGene(0);
  int tempVarCounter = 0;
  std::string evalCode = generateNodeCode(rootGene, tempVarCounter);
  
  out << "  returnValue = " << evalCode << ";\n\n";
  
  // Add range limiting
  out << "  const double RANGELIMIT = 1000000.0;\n";
  out << "  if (returnValue < -RANGELIMIT) return -RANGELIMIT;\n";
  out << "  if (returnValue > RANGELIMIT) return RANGELIMIT;\n";
  out << "  if (std::abs(returnValue) < 0.000001) return 0.0;\n";
  out << "  return returnValue;\n";
  out << "}\n";
  
  out.close();
  std::cout << "Generated C++ evaluator: " << outputFile << std::endl;
}

int main(int argc, char** argv) {
  // Parse command line arguments
  static struct option long_options[] = {
    {"keyname", required_argument, 0, 'k'},
    {"generation", required_argument, 0, 'g'},
    {"codegen", no_argument, 0, 'c'},
    {"output", required_argument, 0, 'o'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };
  
  std::string computedKeyName = "";
  int specifiedGeneration = -1;
  bool generateCode = false;
  std::string outputFile = "gp_evaluator.cpp";
  int option_index = 0;
  int c;
  
  while ((c = getopt_long(argc, argv, "k:g:co:h", long_options, &option_index)) != -1) {
    switch (c) {
      case 'k':
        computedKeyName = optarg;
        break;
      case 'g':
        specifiedGeneration = std::stoi(optarg);
        break;
      case 'c':
        generateCode = true;
        break;
      case 'o':
        outputFile = optarg;
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
  
  // AWS setup
  Aws::SDKOptions options;
  Aws::InitAPI(options);

  std::shared_ptr<Aws::S3::S3Client> s3_client = getS3Client();
  
  // should we look up the latest run?
  if (computedKeyName.empty()) {
    Aws::S3::Model::ListObjectsV2Request listFolders;
    listFolders.SetBucket("autoc-storage"); // TODO extraCfg
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
    listItem.SetBucket("autoc-storage"); // TODO extraCfg
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
  request.SetBucket("autoc-storage"); // TODO extraCfg
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
    std::istringstream iss(retrievedData);
    boost::archive::text_iarchive ia(iss);
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

  if (generateCode) {
    // Generate C++ evaluator code
    int generation = extractGenNumber(keyName);
    generateCppEvaluator(gp, outputFile, keyName, generation);
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
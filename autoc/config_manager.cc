#include "config_manager.h"
#include <vector>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <cstring>
#include "aircraft_state.h"
#include "pathgen.h"
#include "autoc.h"

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/client/ClientConfiguration.h>

// Static member definitions
GPVariables* ConfigManager::gpConfig = nullptr;
ExtraConfig* ConfigManager::extraConfig = nullptr;
GPConfiguration* ConfigManager::gpConfiguration = nullptr;
bool ConfigManager::initialized = false;

void ConfigManager::initialize(const std::string& filename, std::ostream& out) {
    if (initialized) {
        return; // Already initialized
    }
    
    // Check if config file exists
    std::ifstream configFile(filename);
    if (!configFile.good()) {
        out << "FATAL ERROR: Configuration file '" << filename << "' not found!" << std::endl;
        out << "Please ensure autoc.ini exists in the current directory." << std::endl;
        exit(1);
    }
    configFile.close();
    
    // Create config instances
    gpConfig = new GPVariables();
    extraConfig = new ExtraConfig();
    
    // Create the configuration array (moved from autoc.cc)
    GPConfigVarInformation* configArray = createConfigArray();
    
    // Use GPConfiguration to parse the file
    gpConfiguration = new GPConfiguration(out, const_cast<char*>(filename.c_str()), configArray);
    
    // Print S3 configuration once during initialization
    if (strcmp(extraConfig->s3Profile, "default") != 0) {
        out << "S3 Configuration: Using MinIO S3 (profile: " << extraConfig->s3Profile << ", bucket: " << extraConfig->s3Bucket << ")" << std::endl;
    } else {
        out << "S3 Configuration: Using AWS S3 (profile: " << extraConfig->s3Profile << ", bucket: " << extraConfig->s3Bucket << ")" << std::endl;
    }
    
    initialized = true;
}

GPVariables& ConfigManager::getGPConfig() {
    if (!initialized) {
        std::cerr << "Error: ConfigManager not initialized. Call ConfigManager::initialize() first." << std::endl;
        static GPVariables defaultConfig;
        return defaultConfig;
    }
    return *gpConfig;
}

ExtraConfig& ConfigManager::getExtraConfig() {
    if (!initialized) {
        std::cerr << "Error: ConfigManager not initialized. Call ConfigManager::initialize() first." << std::endl;
        static ExtraConfig defaultConfig;
        return defaultConfig;
    }
    return *extraConfig;
}

bool ConfigManager::isInitialized() {
    return initialized;
}

// Create the configuration array (moved from autoc.cc lines 220-248)
GPConfigVarInformation* ConfigManager::createConfigArray() {
    static GPConfigVarInformation configArray[] = {
        {"PopulationSize", DATAINT, &gpConfig->PopulationSize},
        {"NumberOfGenerations", DATAINT, &gpConfig->NumberOfGenerations},
        {"CreationType", DATAINT, &gpConfig->CreationType},
        {"CrossoverProbability", DATADOUBLE, &gpConfig->CrossoverProbability},
        {"CreationProbability", DATADOUBLE, &gpConfig->CreationProbability},
        {"MaximumDepthForCreation", DATAINT, &gpConfig->MaximumDepthForCreation},
        {"MaximumDepthForCrossover", DATAINT, &gpConfig->MaximumDepthForCrossover},
        {"SelectionType", DATAINT, &gpConfig->SelectionType},
        {"TournamentSize", DATAINT, &gpConfig->TournamentSize},
        {"DemeticGrouping", DATAINT, &gpConfig->DemeticGrouping},
        {"DemeSize", DATAINT, &gpConfig->DemeSize},
        {"DemeticMigProbability", DATADOUBLE, &gpConfig->DemeticMigProbability},
        {"SwapMutationProbability", DATADOUBLE, &gpConfig->SwapMutationProbability},
        {"ShrinkMutationProbability", DATADOUBLE, &gpConfig->ShrinkMutationProbability},
        {"AddBestToNewPopulation", DATAINT, &gpConfig->AddBestToNewPopulation},
        {"SteadyState", DATAINT, &gpConfig->SteadyState},
        {"SimNumPathsPerGeneration", DATAINT, &extraConfig->simNumPathsPerGen},
        {"EvalThreads", DATAINT, &extraConfig->evalThreads},
        {"PathGeneratorMethod", DATASTRING, &extraConfig->generatorMethod},
        {"MinisimProgram", DATASTRING, &extraConfig->minisimProgram},
        {"MinisimPortOverride", DATAINT, &extraConfig->minisimPortOverride},
        {"S3Bucket", DATASTRING, &extraConfig->s3Bucket},
        {"S3Profile", DATASTRING, &extraConfig->s3Profile},
        {"EvaluateMode", DATAINT, &extraConfig->evaluateMode},
        {"BytecodeFile", DATASTRING, &extraConfig->bytecodeFile},
        {"WindScenarios", DATAINT, &extraConfig->windScenarioCount},
        {"WindSeedBase", DATAINT, &extraConfig->windSeedBase},
        {"WindSeedStride", DATAINT, &extraConfig->windSeedStride},
        {"", DATAINT, NULL}
    };
    return configArray;
}

std::shared_ptr<Aws::S3::S3Client> ConfigManager::getS3Client() {
    if (!initialized) {
        std::cerr << "Error: ConfigManager not initialized. Call ConfigManager::initialize() first." << std::endl;
        return nullptr;
    }
    
    const ExtraConfig& extraCfg = getExtraConfig();
    
    // real S3 or local minio?
    Aws::Client::ClientConfiguration clientConfig;
    Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::RequestDependent;
    
    if (strcmp(extraCfg.s3Profile, "default") != 0) {
        clientConfig.endpointOverride = "http://localhost:9000"; // MinIO server address
        clientConfig.scheme = Aws::Http::Scheme::HTTP; // Use HTTP instead of HTTPS
        clientConfig.verifySSL = false; // Disable SSL verification for local testing
        
        policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Never;
    }
    
    auto credentialsProvider = Aws::MakeShared<Aws::Auth::ProfileConfigFileAWSCredentialsProvider>("CredentialsProvider", extraCfg.s3Profile);
    
    return Aws::MakeShared<Aws::S3::S3Client>("S3Client",
        credentialsProvider,
        clientConfig,
        policy,
        false
    );
}
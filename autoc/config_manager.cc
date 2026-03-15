#include "config_manager.h"
#include <fstream>
#include <cstdlib>
#include <INIReader.h>

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/core/client/ClientConfiguration.h>

// Static member definitions
AutocConfig* ConfigManager::config = nullptr;
bool ConfigManager::initialized = false;

void ConfigManager::initialize(const std::string& filename, std::ostream& out) {
    if (initialized) {
        return;
    }

    // Check if config file exists
    std::ifstream configFile(filename);
    if (!configFile.good()) {
        out << "FATAL ERROR: Configuration file '" << filename << "' not found!" << std::endl;
        out << "Please ensure autoc.ini exists in the current directory." << std::endl;
        exit(1);
    }
    configFile.close();

    INIReader reader(filename);
    if (reader.ParseError() != 0) {
        out << "FATAL ERROR: Cannot parse configuration file '" << filename << "'" << std::endl;
        exit(1);
    }

    config = new AutocConfig();

    // Evolution core (formerly GPVariables)
    config->populationSize = reader.GetInteger("", "PopulationSize", config->populationSize);
    config->numberOfGenerations = reader.GetInteger("", "NumberOfGenerations", config->numberOfGenerations);
    config->creationType = reader.GetInteger("", "CreationType", config->creationType);
    config->crossoverProbability = reader.GetReal("", "CrossoverProbability", config->crossoverProbability);
    config->creationProbability = reader.GetReal("", "CreationProbability", config->creationProbability);
    config->maximumDepthForCreation = reader.GetInteger("", "MaximumDepthForCreation", config->maximumDepthForCreation);
    config->maximumDepthForCrossover = reader.GetInteger("", "MaximumDepthForCrossover", config->maximumDepthForCrossover);
    config->selectionType = reader.GetInteger("", "SelectionType", config->selectionType);
    config->tournamentSize = reader.GetInteger("", "TournamentSize", config->tournamentSize);
    config->demeticGrouping = reader.GetInteger("", "DemeticGrouping", config->demeticGrouping);
    config->demeSize = reader.GetInteger("", "DemeSize", config->demeSize);
    config->demeticMigProbability = reader.GetReal("", "DemeticMigProbability", config->demeticMigProbability);
    config->swapMutationProbability = reader.GetReal("", "SwapMutationProbability", config->swapMutationProbability);
    config->shrinkMutationProbability = reader.GetReal("", "ShrinkMutationProbability", config->shrinkMutationProbability);
    config->addBestToNewPopulation = reader.GetInteger("", "AddBestToNewPopulation", config->addBestToNewPopulation);
    config->steadyState = reader.GetInteger("", "SteadyState", config->steadyState);

    // Simulation
    config->simNumPathsPerGen = reader.GetInteger("", "SimNumPathsPerGeneration", config->simNumPathsPerGen);
    config->evalThreads = reader.GetInteger("", "EvalThreads", config->evalThreads);
    config->generatorMethod = reader.Get("", "PathGeneratorMethod", config->generatorMethod);
    config->minisimProgram = reader.Get("", "MinisimProgram", config->minisimProgram);
    config->minisimPortOverride = static_cast<unsigned short>(reader.GetInteger("", "MinisimPortOverride", config->minisimPortOverride));

    // S3
    config->s3Bucket = reader.Get("", "S3Bucket", config->s3Bucket);
    config->s3Profile = reader.Get("", "S3Profile", config->s3Profile);

    // Eval mode
    config->evaluateMode = reader.GetInteger("", "EvaluateMode", config->evaluateMode);
    config->bytecodeFile = reader.Get("", "BytecodeFile", config->bytecodeFile);

    // Scenarios
    config->windScenarioCount = reader.GetInteger("", "WindScenarios", config->windScenarioCount);
    config->randomPathSeedB = reader.GetInteger("", "RandomPathSeedB", config->randomPathSeedB);
    config->gpSeed = reader.GetInteger("", "GPSeed", config->gpSeed);
    config->trainingNodes = reader.Get("", "TrainingNodes", config->trainingNodes);

    // Entry and wind direction variations
    config->enableEntryVariations = reader.GetInteger("", "EnableEntryVariations", config->enableEntryVariations);
    config->enableWindVariations = reader.GetInteger("", "EnableWindVariations", config->enableWindVariations);
    config->entryHeadingSigma = reader.GetReal("", "EntryHeadingSigma", config->entryHeadingSigma);
    config->entryRollSigma = reader.GetReal("", "EntryRollSigma", config->entryRollSigma);
    config->entryPitchSigma = reader.GetReal("", "EntryPitchSigma", config->entryPitchSigma);
    config->entrySpeedSigma = reader.GetReal("", "EntrySpeedSigma", config->entrySpeedSigma);
    config->windDirectionSigma = reader.GetReal("", "WindDirectionSigma", config->windDirectionSigma);

    // Entry position variations
    config->entryPositionRadiusSigma = reader.GetReal("", "EntryPositionRadiusSigma", config->entryPositionRadiusSigma);
    config->entryPositionAltSigma = reader.GetReal("", "EntryPositionAltSigma", config->entryPositionAltSigma);

    // Variation landscape ramp
    config->variationRampStep = reader.GetInteger("", "VariationRampStep", config->variationRampStep);

    // Variable rabbit speed
    config->rabbitSpeedNominal = reader.GetReal("", "RabbitSpeedNominal", config->rabbitSpeedNominal);
    config->rabbitSpeedSigma = reader.GetReal("", "RabbitSpeedSigma", config->rabbitSpeedSigma);
    config->rabbitSpeedMin = reader.GetReal("", "RabbitSpeedMin", config->rabbitSpeedMin);
    config->rabbitSpeedMax = reader.GetReal("", "RabbitSpeedMax", config->rabbitSpeedMax);
    config->rabbitSpeedCycleMin = reader.GetReal("", "RabbitSpeedCycleMin", config->rabbitSpeedCycleMin);
    config->rabbitSpeedCycleMax = reader.GetReal("", "RabbitSpeedCycleMax", config->rabbitSpeedCycleMax);

    // Neural network evolution
    config->controllerType = reader.Get("", "ControllerType", config->controllerType);
    config->nnMutationSigma = reader.GetReal("", "NNMutationSigma", config->nnMutationSigma);
    config->nnCrossoverAlpha = reader.GetReal("", "NNCrossoverAlpha", config->nnCrossoverAlpha);
    config->nnWeightFile = reader.Get("", "NNWeightFile", config->nnWeightFile);
    config->nnInitMethod = reader.Get("", "NNInitMethod", config->nnInitMethod);

    // Print S3 configuration
    if (config->s3Profile != "default") {
        out << "S3 Configuration: Using MinIO S3 (profile: " << config->s3Profile << ", bucket: " << config->s3Bucket << ")" << std::endl;
    } else {
        out << "S3 Configuration: Using AWS S3 (profile: " << config->s3Profile << ", bucket: " << config->s3Bucket << ")" << std::endl;
    }

    initialized = true;
}

AutocConfig& ConfigManager::getConfig() {
    if (!initialized) {
        std::cerr << "Error: ConfigManager not initialized. Call ConfigManager::initialize() first." << std::endl;
        static AutocConfig defaultConfig;
        return defaultConfig;
    }
    return *config;
}

bool ConfigManager::isInitialized() {
    return initialized;
}

std::shared_ptr<Aws::S3::S3Client> ConfigManager::getS3Client() {
    if (!initialized) {
        std::cerr << "Error: ConfigManager not initialized. Call ConfigManager::initialize() first." << std::endl;
        return nullptr;
    }

    const AutocConfig& cfg = getConfig();

    Aws::Client::ClientConfiguration clientConfig;
    Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::RequestDependent;

    if (cfg.s3Profile != "default") {
        clientConfig.endpointOverride = "http://localhost:9000";
        clientConfig.scheme = Aws::Http::Scheme::HTTP;
        clientConfig.verifySSL = false;
        policy = Aws::Client::AWSAuthV4Signer::PayloadSigningPolicy::Never;
    }

    auto credentialsProvider = Aws::MakeShared<Aws::Auth::ProfileConfigFileAWSCredentialsProvider>(
        "CredentialsProvider", cfg.s3Profile.c_str());

    return Aws::MakeShared<Aws::S3::S3Client>("S3Client",
        credentialsProvider,
        clientConfig,
        policy,
        false
    );
}

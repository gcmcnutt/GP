#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <string>
#include <iostream>
#include <memory>

// Forward declarations
namespace Aws { namespace S3 { class S3Client; } }

// Unified configuration — replaces GPVariables + ExtraConfig.
// All keys parsed from autoc.ini via inih INIReader.
struct AutocConfig {
    // --- Evolution core (formerly in GPVariables) ---
    int populationSize = 500;
    int numberOfGenerations = 50;
    int creationType = 0;
    double crossoverProbability = 90.0;
    double creationProbability = 5.0;
    int maximumDepthForCreation = 6;
    int maximumDepthForCrossover = 17;
    int selectionType = 0;
    int tournamentSize = 7;
    int demeticGrouping = 0;
    int demeSize = 100;
    double demeticMigProbability = 0.0;
    double swapMutationProbability = 5.0;
    double shrinkMutationProbability = 0.0;
    int addBestToNewPopulation = 1;
    int steadyState = 0;

    // --- Simulation ---
    int simNumPathsPerGen = 1;
    std::string generatorMethod = "classic";
    int evalThreads = 1;
    std::string minisimProgram = "../build/minisim";
    unsigned short minisimPortOverride = 0;

    // --- S3 ---
    std::string s3Bucket = "autoc-storage";
    std::string s3Profile = "default";

    // --- Eval mode ---
    int evaluateMode = 0;
    std::string bytecodeFile = "gp_program.dat";

    // --- Scenarios ---
    int windScenarioCount = 1;
    int randomPathSeedB = 67890;
    int gpSeed = -1;
    std::string trainingNodes = "";

    // --- Entry and wind direction variations ---
    int enableEntryVariations = 0;
    int enableWindVariations = 0;
    double entryHeadingSigma = 45.0;
    double entryRollSigma = 22.5;
    double entryPitchSigma = 7.5;
    double entrySpeedSigma = 0.1;
    double windDirectionSigma = 45.0;

    // --- Entry position variations ---
    double entryPositionRadiusSigma = 0.0;
    double entryPositionAltSigma = 0.0;

    // --- Variation landscape ramp ---
    int variationRampStep = 0;

    // --- Variable rabbit speed ---
    double rabbitSpeedNominal = 16.0;
    double rabbitSpeedSigma = 0.0;
    double rabbitSpeedMin = 8.0;
    double rabbitSpeedMax = 25.0;
    double rabbitSpeedCycleMin = 0.5;
    double rabbitSpeedCycleMax = 5.0;

    // --- Neural network evolution ---
    std::string controllerType = "GP";
    double nnMutationSigma = 0.1;
    double nnCrossoverAlpha = -1.0;
    std::string nnWeightFile = "nn_weights.dat";
    std::string nnInitMethod = "xavier";
};

class ConfigManager {
public:
    // Initialize configuration from INI file using inih
    static void initialize(const std::string& filename = "autoc.ini", std::ostream& out = std::cout);

    // Get unified configuration
    static AutocConfig& getConfig();

    // Check if initialized
    static bool isInitialized();

    // Get S3 client instance with proper configuration
    static std::shared_ptr<Aws::S3::S3Client> getS3Client();

private:
    static AutocConfig* config;
    static bool initialized;
};

#endif // CONFIG_MANAGER_H

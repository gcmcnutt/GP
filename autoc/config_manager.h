#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <string>
#include <iostream>
#include <memory>

// Forward declarations
namespace Aws { namespace S3 { class S3Client; } }

// Unified configuration parsed from autoc.ini via inih INIReader.
struct AutocConfig {
    // --- Evolution ---
    int populationSize = 500;
    int numberOfGenerations = 50;
    double crossoverProbability = 90.0;
    double creationProbability = 5.0;
    int tournamentSize = 7;
    double swapMutationProbability = 5.0;
    int addBestToNewPopulation = 1;

    // --- NN-specific ---
    double nnMutationSigma = 0.1;
    double nnCrossoverAlpha = -1.0;
    std::string nnWeightFile = "nn_weights.dat";
    std::string nnInitMethod = "xavier";

    // --- Simulation ---
    int simNumPathsPerGen = 1;
    std::string generatorMethod = "classic";
    int evalThreads = 1;
    std::string minisimProgram = "./build/minisim";
    unsigned short minisimPortOverride = 0;

    // --- S3 ---
    std::string s3Bucket = "autoc-storage";
    std::string s3Profile = "default";

    // --- Eval mode ---
    int evaluateMode = 0;

    // --- Scenarios ---
    int windScenarioCount = 1;
    int randomPathSeedB = 67890;
    int seed = -1;

    // --- Demetic grouping (scenario assignment) ---
    int demeticGrouping = 0;
    int demeSize = 100;
    double demeticMigProbability = 0.0;

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
};

class ConfigManager {
public:
    static void initialize(const std::string& filename = "autoc.ini", std::ostream& out = std::cout);
    static AutocConfig& getConfig();
    static bool isInitialized();
    static std::shared_ptr<Aws::S3::S3Client> getS3Client();

private:
    static AutocConfig* config;
    static bool initialized;
};

#endif // CONFIG_MANAGER_H

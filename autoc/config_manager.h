#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <string>
#include <iostream>
#include <memory>
#include "gp.h"
#include "gpconfig.h"

// Forward declarations
class ExtraConfig;
namespace Aws { namespace S3 { class S3Client; } }

class ConfigManager {
public:
    // Initialize configuration using GPConfiguration system
    static void initialize(const std::string& filename = "autoc.ini", std::ostream& out = std::cout);
    
    // Get global configuration instances
    static GPVariables& getGPConfig();
    static ExtraConfig& getExtraConfig();
    
    // Check if initialized
    static bool isInitialized();
    
    // Get S3 client instance with proper configuration
    static std::shared_ptr<Aws::S3::S3Client> getS3Client();
    
private:
    static GPVariables* gpConfig;
    static ExtraConfig* extraConfig;  
    static GPConfiguration* gpConfiguration;
    static bool initialized;
    
    // Create the configuration array (moved from autoc.cc lines 220-248)
    static GPConfigVarInformation* createConfigArray();
};

#endif // CONFIG_MANAGER_H
#include "blackbox_extractor_simple.h"
#include <iostream>
#include <fstream>

// Create a sample blackbox file for testing
void createSampleBlackboxFile(const std::string& filename) {
    std::ofstream file(filename);
    
    // Write header
    file << "# Blackbox data from iNav flight controller\n";
    file << "loopIteration,time,axisP[0],axisP[1],GPS_coord[0],GPS_coord[1],GPS_altitude,attitude[0],attitude[1],attitude[2]\n";
    
    // Write sample data (simplified format)
    // Format: timestamp, field1, field2, field3, lat, lon, alt, roll, pitch, yaw
    file << "1000,0,0,0,470000000,80000000,1000,100,50,1800\n";
    file << "1200,0,0,0,470000010,80000010,1020,120,60,1810\n";
    file << "1400,0,0,0,470000020,80000020,1040,140,70,1820\n";
    file << "1600,0,0,0,470000030,80000030,1060,160,80,1830\n";
    file << "1800,0,0,0,470000040,80000040,1080,180,90,1840\n";
    
    file.close();
}

int main() {
    std::cout << "Testing Simple Blackbox Extractor..." << std::endl;
    
    // Create sample blackbox file
    std::string testFile = "/tmp/test_blackbox_simple.txt";
    createSampleBlackboxFile(testFile);
    
    // Test the extractor
    BlackboxExtractor extractor;
    
    if (!extractor.loadBlackboxFile(testFile)) {
        std::cerr << "Failed to load blackbox file" << std::endl;
        return 1;
    }
    
    std::cout << "Successfully loaded blackbox file" << std::endl;
    std::cout << "Number of records: " << extractor.getRecordCount() << std::endl;
    std::cout << "Flight duration: " << extractor.getFlightDuration() << " ms" << std::endl;
    
    // Extract aircraft states
    std::vector<SimpleAircraftState> states = extractor.extractAircraftStates();
    
    std::cout << "Extracted " << states.size() << " aircraft states" << std::endl;
    
    // Print first few states
    for (size_t i = 0; i < std::min(states.size(), size_t(3)); ++i) {
        const auto& state = states[i];
        std::cout << "State " << i << ": " 
                  << "Time=" << state.timeMsec << "ms "
                  << "Pos=[" << state.position.transpose() << "] "
                  << "Vel=" << state.relativeVelocity << " "
                  << "Pitch=" << state.pitchCommand << " "
                  << "Roll=" << state.rollCommand << " "
                  << "Throttle=" << state.throttleCommand << std::endl;
    }
    
    // Clean up
    std::remove(testFile.c_str());
    
    std::cout << "Test completed successfully!" << std::endl;
    return 0;
}
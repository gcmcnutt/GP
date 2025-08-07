#ifndef GP_EVALUATOR_EMBEDDED_H
#define GP_EVALUATOR_EMBEDDED_H

#include "gp_evaluator_portable.h"

// Embedded/Arduino-specific extensions
class GPEvaluatorEmbedded {
public:
    // Simple evaluation function for Arduino use - single path, single GP program
    static double evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, 
                                  const GPBytecode* program, int program_size, double arg = 0.0);
    
    // Convert GP control outputs (-1 to 1) to MSP channel values (1000-2000)
    static int convertToMSPChannel(double gp_command);
    
    // Helper to create single path provider
    static SinglePathProvider createPathProvider(const Path& currentPath, int currentIndex);
};

// Global GP program management for embedded systems
void setEmbeddedGPProgram(const GPBytecode* program, int size);

// C-style function for easy Arduino integration
extern "C" {
    double evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, double arg);
}

#endif
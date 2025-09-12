#ifndef GP_EVALUATOR_EMBEDDED_H
#define GP_EVALUATOR_EMBEDDED_H

#include "gp_evaluator_portable.h"

// Embedded/Arduino-specific extensions
// Note: With static bytecode arrays, most embedded-specific functionality
// is no longer needed. The portable evaluator handles everything.

class GPEvaluatorEmbedded {
public:
    // Simple evaluation function for embedded use - delegates to portable evaluator
    static double evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, 
                                  const GPBytecode* program, int program_size, double arg = 0.0);
    
    // Helper to create single path provider
    static SinglePathProvider createPathProvider(const Path& currentPath, int currentIndex);
};

#endif
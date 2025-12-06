#ifndef GP_EVALUATOR_DESKTOP_H
#define GP_EVALUATOR_DESKTOP_H

#include "gp_evaluator_portable.h"
#include "gp_bytecode.h"

// Desktop-specific extensions for GP build
class GPEvaluatorDesktop {
public:
    // Convert GP build GPBytecode to portable format (no conversion needed in GP build)
    static void convertBytecode(const std::vector<struct GPBytecode>& gpBytecode, 
                               std::vector<struct GPBytecode>& portableBytecode);
    
    // Evaluate using GP build bytecode format
    static double evaluateGPBytecode(const std::vector<struct GPBytecode>& program, 
                                    std::vector<Path>& path, AircraftState& aircraftState, 
                                    double contextArg = 0.0);

    static double evaluateGPBytecodeReference(const std::vector<struct GPBytecode>& program,
                                              std::vector<Path>& path, AircraftState& aircraftState,
                                              double contextArg = 0.0);
};

#endif

#include "gp_evaluator_desktop.h"

void GPEvaluatorDesktop::convertBytecode(const std::vector<struct GPBytecode>& gpBytecode, 
                                        std::vector<struct GPBytecode>& portableBytecode) {
    // For GP build, both formats are the same - just copy
    portableBytecode = gpBytecode;
}

double GPEvaluatorDesktop::evaluateGPBytecode(const std::vector<struct GPBytecode>& program, 
                                             std::vector<Path>& path, AircraftState& aircraftState, 
                                             double contextArg) {
    // Create path provider and evaluate using desktop implementation
    VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
    return evaluateBytecodePortable(program.data(), program.size(), 
                                   pathProvider, aircraftState, contextArg);
}

// Desktop evaluator is now just a wrapper around the portable implementation
// The actual implementation is in gp_evaluator_portable.cc
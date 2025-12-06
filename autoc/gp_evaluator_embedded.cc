#include "gp_evaluator_embedded.h"

double GPEvaluatorEmbedded::evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, 
                                           const GPBytecode* program, int program_size, double arg) {
    // Create single-path provider for embedded use
    SinglePathProvider pathProvider(currentPath, aircraftState.getThisPathIndex());
    
    // Use portable bytecode evaluator - this is the core evaluation path
    return evaluateBytecodePortable(program, program_size, pathProvider, aircraftState, fastmath::GPScalar::fromDouble(arg)).toDouble();
}

SinglePathProvider GPEvaluatorEmbedded::createPathProvider(const Path& currentPath, int currentIndex) {
    return SinglePathProvider(currentPath, currentIndex);
}

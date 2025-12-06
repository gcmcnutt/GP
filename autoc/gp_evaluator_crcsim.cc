#include "gp_evaluator_crcsim.h"

double GPEvaluatorCRRCSim::evaluateBytecodeForCRRCSim(const GPBytecode* program, int program_size,
                                                     std::vector<Path>& path, AircraftState& aircraftState, 
                                                     double arg) {
    // Create vector path provider for CRRCSim
    VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
    
    // Use portable bytecode evaluator
    return evaluateBytecodePortable(program, program_size, pathProvider, aircraftState, fastmath::GPScalar::fromDouble(arg)).toDouble();
}

void GPEvaluatorCRRCSim::convertFromBoostBytecode(const void* boostBytecode, int boost_size,
                                                 std::vector<GPBytecode>& portableBytecode) {
    // This would handle conversion from boost serialized bytecode
    // For now, placeholder implementation
    // TODO: Implement boost bytecode deserialization for CRRCSim
    portableBytecode.clear();
}

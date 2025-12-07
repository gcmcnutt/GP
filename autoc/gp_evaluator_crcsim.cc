#include "gp_evaluator_crcsim.h"

gp_scalar GPEvaluatorCRRCSim::evaluateBytecodeForCRRCSim(const GPBytecode* program, int program_size,
                                                     std::vector<Path>& path, AircraftState& aircraftState, 
                                                     gp_scalar arg) {
    // Create vector path provider for CRRCSim
    VectorPathProvider pathProvider(path, aircraftState.getThisPathIndex());
    
    // Use portable bytecode evaluator
    return evaluateBytecodePortable(program, program_size, pathProvider, aircraftState, arg);
}

void GPEvaluatorCRRCSim::convertFromBoostBytecode(const void* boostBytecode, int boost_size,
                                                 std::vector<GPBytecode>& portableBytecode) {
    // This would handle conversion from boost serialized bytecode
    // For now, placeholder implementation
    // TODO: Implement boost bytecode deserialization for CRRCSim
    portableBytecode.clear();
}

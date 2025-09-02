#include "gp_evaluator_embedded.h"

double GPEvaluatorEmbedded::evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, 
                                           const GPBytecode* program, int program_size, double arg) {
    // Create single-path provider for embedded use
    SinglePathProvider pathProvider(currentPath, aircraftState.getThisPathIndex());
    
    // Use portable bytecode evaluator
    return evaluateBytecodePortable(program, program_size, pathProvider, aircraftState, arg);
}

SinglePathProvider GPEvaluatorEmbedded::createPathProvider(const Path& currentPath, int currentIndex) {
    return SinglePathProvider(currentPath, currentIndex);
}

// Global GP program storage for embedded systems
static GPBytecode* embedded_gp_program = nullptr;
static int embedded_gp_program_size = 0;

// Set the global GP program for embedded evaluation
void setEmbeddedGPProgram(const GPBytecode* program, int size) {
    embedded_gp_program = const_cast<GPBytecode*>(program);
    embedded_gp_program_size = size;
}

// C-style function for easy Arduino integration
extern "C" double evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, double arg) {
    if (embedded_gp_program && embedded_gp_program_size > 0) {
        return GPEvaluatorEmbedded::evaluateGPSimple(aircraftState, currentPath, 
                                                    embedded_gp_program, embedded_gp_program_size, arg);
    }
    
    // Placeholder: simple proportional control until GP program is loaded
    // This gives basic flight capability for testing
    Eigen::Vector3d toTarget = currentPath.start - aircraftState.getPosition();
    double distance = toTarget.norm();
    
    if (distance > 1.0) {
        // Simple proportional control toward target
        Eigen::Vector3d targetLocal = aircraftState.getOrientation().inverse() * toTarget.normalized();
        aircraftState.setRollCommand(CLAMP_DEF(targetLocal.y() * 2.0, -1.0, 1.0));
        aircraftState.setPitchCommand(CLAMP_DEF(-targetLocal.z() * 1.0, -1.0, 1.0));
        aircraftState.setThrottleCommand(0.5); // Constant throttle
    }
    
    return distance;
}
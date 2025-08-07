#include "gp_evaluator_embedded.h"

double GPEvaluatorEmbedded::evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, 
                                           const GPBytecode* program, int program_size, double arg) {
    // Create single-path provider for embedded use
    SinglePathProvider pathProvider(currentPath, aircraftState.getThisPathIndex());
    
    // Use portable bytecode evaluator
    return evaluateBytecodePortable(program, program_size, pathProvider, aircraftState, arg);
}

// Embedded implementation of bytecode evaluation (only compiled when not GP_BUILD)
#ifndef GP_BUILD
double evaluateBytecodePortable(const GPBytecode* program, int program_size, 
                               PathProvider& pathProvider, AircraftState& aircraftState, 
                               double contextArg) {
    const int MAX_STACK_SIZE = 64;  // Smaller stack for embedded
    double stack[MAX_STACK_SIZE];
    int stack_ptr = 0;
    
    // Execute bytecode instructions
    for (int i = 0; i < program_size; i++) {
        const GPBytecode& instruction = program[i];
        
        switch (instruction.opcode) {
            // Binary operations - pop two, push result
            case ADD:
            case SUB:
            case MUL:
            case DIV:
            case EQ:
            case GT:
            case ATAN2:
            case MIN:
            case MAX: {
                if (stack_ptr < 2) return 0.0;
                double args[2] = {stack[stack_ptr-2], stack[stack_ptr-1]};
                stack_ptr -= 2;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 2);
                break;
            }
            
            // Unary operations - pop one, push result
            case SIN:
            case COS:
            case SETPITCH:
            case SETROLL:
            case SETTHROTTLE:
            case ABS:
            case SQRT:
            case GETDPHI:
            case GETDTHETA:
            case GETDTARGET: {
                if (stack_ptr < 1) return 0.0;
                double args[1] = {stack[stack_ptr-1]};
                stack_ptr -= 1;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 1);
                break;
            }
            
            // Ternary operations - pop three, push result
            case IF:
            case CLAMP: {
                if (stack_ptr < 3) return 0.0;
                double args[3] = {stack[stack_ptr-3], stack[stack_ptr-2], stack[stack_ptr-1]};
                stack_ptr -= 3;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 3);
                break;
            }
            
            // PROGN - pop two, discard first, keep second
            case PROGN: {
                if (stack_ptr < 2) return 0.0;
                double args[2] = {stack[stack_ptr-2], stack[stack_ptr-1]};
                stack_ptr -= 2;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 2);
                break;
            }
            
            // Zero-argument operations (terminals and sensors)
            default: {
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, nullptr, 0);
                break;
            }
        }
        
        // Stack overflow protection
        if (stack_ptr >= MAX_STACK_SIZE) {
            return 0.0;
        }
    }
    
    // Return final result
    if (stack_ptr != 1) {
        return 0.0;
    }
    
    return applyRangeLimit(stack[0]);
}
#endif

int GPEvaluatorEmbedded::convertToMSPChannel(double gp_command) {
    // Clamp GP command to [-1, 1] and map to MSP range [1000, 2000]
    double clamped = CLAMP_DEF(gp_command, -1.0, 1.0);
    return (int)(1500.0 + clamped * 500.0);
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
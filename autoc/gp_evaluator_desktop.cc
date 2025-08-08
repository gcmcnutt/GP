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

// Desktop implementation of portable bytecode evaluation
double evaluateBytecodePortable(const struct GPBytecode* program, int program_size, 
                               PathProvider& pathProvider, AircraftState& aircraftState, 
                               double contextArg) {
    const int MAX_STACK_SIZE = 256;
    double stack[MAX_STACK_SIZE];
    int stack_ptr = 0;
    
    // Execute bytecode instructions
    for (int i = 0; i < program_size; i++) {
        const struct GPBytecode& instruction = program[i];
        
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
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 2, contextArg);
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
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 1, contextArg);
                break;
            }
            
            // Ternary operations - pop three, push result
            case IF:
            case CLAMP: {
                if (stack_ptr < 3) return 0.0;
                double args[3] = {stack[stack_ptr-3], stack[stack_ptr-2], stack[stack_ptr-1]};
                stack_ptr -= 3;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 3, contextArg);
                break;
            }
            
            // PROGN - pop two, discard first, keep second
            case PROGN: {
                if (stack_ptr < 2) return 0.0;
                double args[2] = {stack[stack_ptr-2], stack[stack_ptr-1]};
                stack_ptr -= 2;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 2, contextArg);
                break;
            }
            
            // Zero-argument operations (terminals and sensors)
            default: {
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, nullptr, 0, contextArg);
                break;
            }
        }
        
        // Stack overflow protection
        if (stack_ptr >= MAX_STACK_SIZE) {
            std::cerr << "Error: Stack overflow in bytecode execution" << std::endl;
            return 0.0;
        }
    }
    
    // Return final result
    if (stack_ptr != 1) {
        std::cerr << "Error: Invalid stack state after bytecode execution (stack_ptr=" << stack_ptr << ")" << std::endl;
        return 0.0;
    }
    
    return applyRangeLimit(stack[0]);
}
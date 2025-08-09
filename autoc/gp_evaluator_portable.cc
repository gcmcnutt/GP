#include "gp_evaluator_portable.h"
#include <cmath>

#ifdef GP_BUILD
#include "gp_bytecode.h"  // Full GPBytecode definition for GP builds
#endif

#ifdef GP_BUILD
#include <iostream>
#endif

double evaluateGPOperator(int opcode, PathProvider& pathProvider, 
                         AircraftState& aircraftState,
                         const double* args, int argc, double contextArg) {
    double result = 0.0;
    
    switch (opcode) {
        // Math operators - identical on all platforms
        case ADD: 
            result = args[0] + args[1]; 
            break;
        case SUB: 
            result = -args[0] - args[1]; // autoc convention
            break;
        case MUL: 
            result = args[0] * args[1]; 
            break;
        case DIV: 
            result = (args[1] == 0) ? 0 : args[0] / args[1]; 
            break;
        
        // Control operators - use CLAMP_DEF macro (platform-specific)
        case SETPITCH: 
            result = aircraftState.setPitchCommand(args[0]); 
            break;
        case SETROLL: 
            result = aircraftState.setRollCommand(args[0]); 
            break;
        case SETTHROTTLE: 
            result = aircraftState.setThrottleCommand(args[0]); 
            break;
        
        // State queries - identical logic
        case GETPITCH: 
            result = aircraftState.getPitchCommand(); 
            break;
        case GETROLL: 
            result = aircraftState.getRollCommand(); 
            break;
        case GETTHROTTLE: 
            result = aircraftState.getThrottleCommand(); 
            break;
        case GETVEL: 
            result = aircraftState.getRelVel(); 
            break;
        
        // Navigation - use PathProvider abstraction
        case GETDPHI: 
            result = executeGetDPhi(pathProvider, aircraftState, args ? args[0] : contextArg); 
            break;
        case GETDTHETA: 
            result = executeGetDTheta(pathProvider, aircraftState, args ? args[0] : contextArg); 
            break;
        case GETDTARGET: 
            result = executeGetDTarget(pathProvider, aircraftState, args ? args[0] : contextArg); 
            break;
        case GETDHOME: 
            result = executeGetDHome(aircraftState); 
            break;
        
        // Trigonometry - use C math library (available on all platforms)
        case SIN: 
            result = sin(args[0]); 
            break;
        case COS: 
            result = cos(args[0]); 
            break;
        case ATAN2: 
            result = ATAN2_DEF(args[0], args[1]); 
            break;
        
        // Math helpers - use platform macros  
        case CLAMP: 
            result = CLAMP_DEF(args[0], args[1], args[2]); 
            break;
        case ABS: 
            result = ABS_DEF(args[0]); 
            break;
        case SQRT: 
            result = (args[0] >= 0) ? SQRT_DEF(args[0]) : 0.0; 
            break;
        case MIN: 
            result = MIN_DEF(args[0], args[1]); 
            break;
        case MAX: 
            result = MAX_DEF(args[0], args[1]); 
            break;
        
        // Logical operators
        case IF: 
            result = args[0] ? args[1] : args[2]; 
            break;
        case EQ: 
            result = (args[0] == args[1]) ? 1.0 : 0.0; 
            break;
        case GT: 
            result = (args[0] > args[1]) ? 1.0 : 0.0; 
            break;
        case PROGN: 
            result = args[1]; // Return second arg, ignore first
            break;
        
        // Constants
        case OP_PI:
            result = M_PI; 
            break;
        case ZERO: 
            result = 0.0; 
            break;
        case ONE: 
            result = 1.0; 
            break;
        case TWO: 
            result = 2.0; 
            break;
        
        // Velocity/attitude sensors
        case GETVELX: 
            result = aircraftState.getVelocity().x(); 
            break;
        case GETVELY: 
            result = aircraftState.getVelocity().y(); 
            break;
        case GETVELZ: 
            result = aircraftState.getVelocity().z(); 
            break;
        
        case GETALPHA: {
            // Transform actual velocity vector to body frame
            Eigen::Vector3d velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
            // Angle of attack is angle between velocity and body X axis (forward)
            result = ATAN2_DEF(-velocity_body.z(), velocity_body.x());
            break;
        }
        
        case GETBETA: {
            // Transform actual velocity vector to body frame
            Eigen::Vector3d velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
            // Sideslip is angle between velocity and body XZ plane
            result = ATAN2_DEF(velocity_body.y(), velocity_body.x());
            break;
        }
        
        case GETROLL_RAD: {
            // Convert quaternion to Euler angles using standard aerospace convention
            Eigen::Vector3d euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
            result = euler[2]; // Roll angle (rotation around X-axis)
            break;
        }
        
        case GETPITCH_RAD: {
            // Convert quaternion to Euler angles using standard aerospace convention
            Eigen::Vector3d euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
            result = euler[1]; // Pitch angle (rotation around Y-axis)
            break;
        }
        
        default:
#ifdef GP_BUILD
            std::cerr << "Unknown operator: " << opcode << std::endl;
#endif
            result = 0.0;
            break;
    }
    
    return applyRangeLimit(result);
}

double executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, double arg) {
    // Calculate the vector from craft to target in world frame
    int idx = getPathIndex(pathProvider, aircraftState, arg);
    Eigen::Vector3d craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    
    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    
    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());
    
    // Calculate the angle between the projected vector and the body Z-axis
    return ATAN2_DEF(projectedVector.y(), -projectedVector.z());
}

double executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, double arg) {
    // Calculate the vector from craft to target in world frame
    int idx = getPathIndex(pathProvider, aircraftState, arg);
    Eigen::Vector3d craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    
    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    
    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());
    
    // Calculate the angle between the projected vector and the body Z-axis
    double rollEstimate = ATAN2_DEF(projectedVector.y(), -projectedVector.z());
    
    // *** PITCH: Calculate the vector from craft to target in world frame if it did rotate
    Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
    Eigen::Quaterniond virtualOrientation = aircraftState.getOrientation() * rollRotation;
    
    // Transform target vector to new virtual orientation
    Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;
    
    // Calculate pitch angle
    return ATAN2_DEF(-newLocalTargetVector.z(), newLocalTargetVector.x());
}

double executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, double arg) {
    int idx = getPathIndex(pathProvider, aircraftState, arg);
    double distance = (pathProvider.getPath(idx).start - aircraftState.getPosition()).norm();
    return CLAMP_DEF((distance - 10) / aircraftState.getRelVel(), -1.0, 1.0);
}

double executeGetDHome(AircraftState& aircraftState) {
    return (Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE) - aircraftState.getPosition()).norm();
}

// Portable bytecode evaluation implementation - works on all platforms
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
#ifdef GP_BUILD
            std::cerr << "Error: Stack overflow in bytecode execution" << std::endl;
#endif
            return 0.0;
        }
    }
    
    // Return final result
    if (stack_ptr != 1) {
#ifdef GP_BUILD
        std::cerr << "Error: Invalid stack state after bytecode execution (stack_ptr=" << stack_ptr << ")" << std::endl;
#endif
        return 0.0;
    }
    
    return applyRangeLimit(stack[0]);
}
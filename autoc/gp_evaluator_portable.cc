#include "gp_evaluator_portable.h"
#include <array>
#include <cmath>

#if defined(GP_BUILD) && !defined(GP_TEST)
#include "gp_bytecode.h"  // Full GPBytecode definition for GP builds
#endif

#ifdef GP_BUILD
#include <iostream>
#endif

namespace {
constexpr gp_scalar GP_PI = static_cast<gp_scalar>(M_PI);
constexpr gp_scalar GP_TWO_PI = GP_PI * static_cast<gp_scalar>(2.0f);
constexpr gp_scalar GP_HALF_PI = GP_PI * static_cast<gp_scalar>(0.5f);

// Sine/cosine lookup
constexpr int SIN_COS_LUT_SIZE = 512;
std::array<gp_scalar, SIN_COS_LUT_SIZE + 1> SIN_LUT{};
bool SIN_LUT_INIT = false;

inline gp_scalar wrapAngle(gp_scalar angle) {
    gp_scalar wrapped = std::fmod(angle, GP_TWO_PI);
    if (wrapped < 0) wrapped += GP_TWO_PI;
    return wrapped;
}

inline void initTrigLut() {
    if (SIN_LUT_INIT) return;
    for (int i = 0; i <= SIN_COS_LUT_SIZE; ++i) {
        gp_scalar angle = GP_TWO_PI * static_cast<gp_scalar>(i) / static_cast<gp_scalar>(SIN_COS_LUT_SIZE);
        SIN_LUT[i] = std::sin(angle);
    }
    SIN_LUT_INIT = true;
}

inline gp_scalar fastSin(gp_scalar angle) {
    initTrigLut();
    gp_scalar wrapped = wrapAngle(angle);
    gp_scalar scaled = wrapped * (static_cast<gp_scalar>(SIN_COS_LUT_SIZE) / GP_TWO_PI);
    int i0 = static_cast<int>(scaled);
    gp_scalar frac = scaled - static_cast<gp_scalar>(i0);
    int i1 = (i0 + 1 <= SIN_COS_LUT_SIZE) ? i0 + 1 : 0;
    return SIN_LUT[i0] + (SIN_LUT[i1] - SIN_LUT[i0]) * frac;
}

inline gp_scalar fastCos(gp_scalar angle) {
    return fastSin(angle + GP_HALF_PI);
}

// atan lookup for ratio in [0,1], mirrored for other quadrants
constexpr int ATAN_LUT_SIZE = 512;
std::array<gp_scalar, ATAN_LUT_SIZE + 1> ATAN_LUT{};
bool ATAN_LUT_INIT = false;

inline void initAtanLut() {
    if (ATAN_LUT_INIT) return;
    for (int i = 0; i <= ATAN_LUT_SIZE; ++i) {
        gp_scalar t = static_cast<gp_scalar>(i) / static_cast<gp_scalar>(ATAN_LUT_SIZE);
        ATAN_LUT[i] = std::atan(t);
    }
    ATAN_LUT_INIT = true;
}

inline gp_scalar fastAtan(gp_scalar r) {
    initAtanLut();
    gp_scalar t = ABS_DEF(r);
    gp_scalar scaled = t * static_cast<gp_scalar>(ATAN_LUT_SIZE);
    if (scaled >= ATAN_LUT_SIZE) {
        return (r < 0.0f ? -ATAN_LUT[ATAN_LUT_SIZE] : ATAN_LUT[ATAN_LUT_SIZE]);
    }
    int i0 = static_cast<int>(scaled);
    gp_scalar frac = scaled - static_cast<gp_scalar>(i0);
    gp_scalar a0 = ATAN_LUT[i0];
    gp_scalar a1 = ATAN_LUT[i0 + 1];
    gp_scalar angle = a0 + (a1 - a0) * frac;
    return (r < 0.0f) ? -angle : angle;
}

inline gp_scalar fastAtan2(gp_scalar y, gp_scalar x) {
    if (ABS_DEF(x) < static_cast<gp_scalar>(1e-6f) && ABS_DEF(y) < static_cast<gp_scalar>(1e-6f)) {
        return 0.0f;
    }

    if (ABS_DEF(x) > ABS_DEF(y)) {
        gp_scalar angle = fastAtan(y / x);
        if (x < 0.0f) {
            angle += (y >= 0.0f ? GP_PI : -GP_PI);
        }
        return angle;
    } else {
        gp_scalar angle = fastAtan(x / y);
        gp_scalar base = (y > 0.0f) ? GP_HALF_PI : -GP_HALF_PI;
        return base - angle;
    }
}
}  // namespace

gp_scalar evaluateGPOperator(int opcode, PathProvider& pathProvider, 
                         AircraftState& aircraftState,
                         const gp_scalar* args, int argc, gp_scalar contextArg) {
    gp_scalar result = 0.0f;
    
    switch (opcode) {
        // Math operators - identical on all platforms
        case ADD: 
            result = args[0] + args[1]; 
            break;
        case SUB: 
            result = args[0] - args[1];
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
            result = fastSin(args[0]); 
            break;
        case COS: 
            result = fastCos(args[0]); 
            break;
        case ATAN2: 
            result = fastAtan2(args[0], args[1]); 
            break;
        
        // Math helpers - use platform macros  
        case CLAMP: 
            result = CLAMP_DEF(args[0], args[1], args[2]); 
            break;
        case ABS: 
            result = ABS_DEF(args[0]); 
            break;
        case SQRT: 
            result = (args[0] >= 0) ? SQRT_DEF(args[0]) : 0.0f; 
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
            result = (args[0] == args[1]) ? 1.0f : 0.0f; 
            break;
        case GT: 
            result = (args[0] > args[1]) ? 1.0f : 0.0f; 
            break;
        case PROGN: 
            result = args[1]; // Return second arg, ignore first
            break;
        
        // Constants
        case OP_PI:
            result = GP_PI; 
            break;
        case ZERO: 
            result = 0.0f; 
            break;
        case ONE: 
            result = 1.0f; 
            break;
        case TWO: 
            result = 2.0f; 
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
            gp_vec3 velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
            // Angle of attack is angle between velocity and body X axis (forward)
            result = fastAtan2(-velocity_body.z(), velocity_body.x());
            break;
        }
        
        case GETBETA: {
            // Transform actual velocity vector to body frame
            gp_vec3 velocity_body = aircraftState.getOrientation().inverse() * aircraftState.getVelocity();
            // Sideslip is angle between velocity and body XZ plane
            result = fastAtan2(velocity_body.y(), velocity_body.x());
            break;
        }
        
        case GETROLL_RAD: {
            // Convert quaternion to Euler angles using standard aerospace convention
            gp_vec3 euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
            result = euler[2]; // Roll angle (rotation around X-axis)
            break;
        }
        
        case GETPITCH_RAD: {
            // Convert quaternion to Euler angles using standard aerospace convention
            gp_vec3 euler = aircraftState.getOrientation().toRotationMatrix().eulerAngles(2, 1, 0);
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

gp_scalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {
    // Calculate the vector from craft to target in world frame
    int idx = getPathIndex(pathProvider, aircraftState, arg);
    gp_vec3 craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    
    // Transform the craft-to-target vector to body frame
    gp_vec3 target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    
    // Project the craft-to-target vector onto the body YZ plane
    gp_vec3 projectedVector(0.0f, target_local.y(), target_local.z());
    
    // Calculate the angle between the projected vector and the body Z-axis
    return fastAtan2(projectedVector.y(), -projectedVector.z());
}

gp_scalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {
    // Calculate the vector from craft to target in world frame
    int idx = getPathIndex(pathProvider, aircraftState, arg);
    gp_vec3 craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    
    // Transform the craft-to-target vector to body frame
    gp_vec3 target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    
    // Project the craft-to-target vector onto the body YZ plane
    gp_vec3 projectedVector(0.0f, target_local.y(), target_local.z());
    
    // Calculate the angle between the projected vector and the body Z-axis
    gp_scalar rollEstimate = fastAtan2(projectedVector.y(), -projectedVector.z());
    
    // *** PITCH: Calculate the vector from craft to target in world frame if it did rotate
    gp_quat rollRotation(Eigen::AngleAxis<gp_scalar>(rollEstimate, gp_vec3::UnitX()));
    gp_quat virtualOrientation = aircraftState.getOrientation() * rollRotation;
    
    // Transform target vector to new virtual orientation
    gp_vec3 newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;
    
    // Calculate pitch angle
    return fastAtan2(-newLocalTargetVector.z(), newLocalTargetVector.x());
}

gp_scalar executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg) {
    int idx = getPathIndex(pathProvider, aircraftState, arg);
    gp_scalar distance = (pathProvider.getPath(idx).start - aircraftState.getPosition()).norm();
    return CLAMP_DEF((distance - static_cast<gp_scalar>(10.0f)) / aircraftState.getRelVel(), -1.0f, 1.0f);
}

gp_scalar executeGetDHome(AircraftState& aircraftState) {
    return (gp_vec3(0.0f, 0.0f, SIM_INITIAL_ALTITUDE) - aircraftState.getPosition()).norm();
}

// Portable bytecode evaluation implementation - works on all platforms
gp_scalar evaluateBytecodePortable(const struct GPBytecode* program, int program_size, 
                               PathProvider& pathProvider, AircraftState& aircraftState, 
                               gp_scalar contextArg) {
    const int MAX_STACK_SIZE = 256;
    gp_scalar stack[MAX_STACK_SIZE];
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
                if (stack_ptr < 2) return 0.0f;
                gp_scalar args[2] = {stack[stack_ptr-2], stack[stack_ptr-1]};
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
                if (stack_ptr < 1) return 0.0f;
                gp_scalar args[1] = {stack[stack_ptr-1]};
                stack_ptr -= 1;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 1, contextArg);
                break;
            }
            
            // Ternary operations - pop three, push result
            case IF:
            case CLAMP: {
                if (stack_ptr < 3) return 0.0f;
                gp_scalar args[3] = {stack[stack_ptr-3], stack[stack_ptr-2], stack[stack_ptr-1]};
                stack_ptr -= 3;
                stack[stack_ptr++] = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, args, 3, contextArg);
                break;
            }
            
            // PROGN - pop two, discard first, keep second
            case PROGN: {
                if (stack_ptr < 2) return 0.0f;
                gp_scalar args[2] = {stack[stack_ptr-2], stack[stack_ptr-1]};
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
            return 0.0f;
        }
    }
    
    // Return final result
    if (stack_ptr != 1) {
#ifdef GP_BUILD
            std::cerr << "Error: Invalid stack state after bytecode execution (stack_ptr=" << stack_ptr << ")" << std::endl;
#endif
        return 0.0f;
    }
    
    return applyRangeLimit(stack[0]);
}

#include "gp_evaluator_portable.h"
#include <cmath>
#include "fastmath/fixed_math.h"
#include "fastmath/orientation_math.h"

#ifdef GP_BUILD
#include "gp_bytecode.h"  // Full GPBytecode definition for GP builds
#endif

#ifdef GP_BUILD
#include <iostream>
#endif

fastmath::GPScalar evaluateGPOperator(int opcode, PathProvider& pathProvider, 
                         AircraftState& aircraftState,
                         const fastmath::GPScalar* args, int argc, fastmath::GPScalar contextArg) {
    using fastmath::GPScalar;
    GPScalar result = GPScalar::zero();
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
            result = args[1].isZero() ? GPScalar::zero() : args[0] / args[1]; 
            break;
        
        // Control operators - use CLAMP_DEF macro (platform-specific)
        case SETPITCH: 
            result = GPScalar::fromDouble(aircraftState.setPitchCommand(args[0].toDouble())); 
            break;
        case SETROLL: 
            result = GPScalar::fromDouble(aircraftState.setRollCommand(args[0].toDouble())); 
            break;
        case SETTHROTTLE: 
            result = GPScalar::fromDouble(aircraftState.setThrottleCommand(args[0].toDouble())); 
            break;
        
        // State queries - identical logic
        case GETPITCH: 
            result = aircraftState.getPitchCommandScalar(); 
            break;
        case GETROLL: 
            result = aircraftState.getRollCommandScalar(); 
            break;
        case GETTHROTTLE: 
            result = aircraftState.getThrottleCommandScalar(); 
            break;
        case GETVEL: 
            result = aircraftState.getRelVelScalar(); 
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
            result = fastmath::sinApprox(args[0]); 
            break;
        case COS: 
            result = fastmath::cosApprox(args[0]); 
            break;
        case ATAN2: 
            result = fastmath::atan2Approx(args[0], args[1]); 
            break;
        
        // Math helpers - use platform macros  
        case CLAMP: 
            result = fastmath::clamp(args[0], args[1], args[2]); 
            break;
        case ABS: 
            result = args[0].abs(); 
            break;
        case SQRT: 
            result = fastmath::sqrtApprox(args[0]); 
            break;
        case MIN: 
            result = fastmath::min(args[0], args[1]); 
            break;
        case MAX: 
            result = fastmath::max(args[0], args[1]); 
            break;
        
        // Logical operators
        case IF: 
            result = !args[0].isZero() ? args[1] : args[2]; 
            break;
        case EQ: 
            result = (args[0] == args[1]) ? GPScalar::fromInt(1) : GPScalar::zero(); 
            break;
        case GT: 
            result = (args[0] > args[1]) ? GPScalar::fromInt(1) : GPScalar::zero(); 
            break;
        case PROGN: 
            result = args[1]; // Return second arg, ignore first
            break;
        
        // Constants
        case OP_PI:
            result = GPScalar::fromDouble(M_PI); 
            break;
        case ZERO: 
            result = GPScalar::zero(); 
            break;
        case ONE: 
            result = GPScalar::fromInt(1); 
            break;
        case TWO: 
            result = GPScalar::fromInt(2); 
            break;
        
        // Velocity/attitude sensors
        case GETVELX: 
            result = GPScalar::fromDouble(aircraftState.getVelocity().x()); 
            break;
        case GETVELY: 
            result = GPScalar::fromDouble(aircraftState.getVelocity().y()); 
            break;
        case GETVELZ: 
            result = GPScalar::fromDouble(aircraftState.getVelocity().z()); 
            break;
        
        case GETALPHA: {
            Eigen::Vector3d velocity_body = fastmath::rotateWorldToBody(
                aircraftState.getOrientation(), aircraftState.getVelocity());
            result = fastmath::atan2Approx(
                GPScalar::fromDouble(-velocity_body.z()),
                GPScalar::fromDouble(velocity_body.x()));
            break;
        }
        
        case GETBETA: {
            Eigen::Vector3d velocity_body = fastmath::rotateWorldToBody(
                aircraftState.getOrientation(), aircraftState.getVelocity());
            result = fastmath::atan2Approx(
                GPScalar::fromDouble(velocity_body.y()),
                GPScalar::fromDouble(velocity_body.x()));
            break;
        }
        
        case GETROLL_RAD: {
            double roll = fastmath::rollFromQuaternion(aircraftState.getOrientation());
            result = GPScalar::fromDouble(roll);
            break;
        }
        
        case GETPITCH_RAD: {
            double pitch = fastmath::pitchFromQuaternion(aircraftState.getOrientation());
            result = GPScalar::fromDouble(pitch);
            break;
        }
        
        default:
#ifdef GP_BUILD
            std::cerr << "Unknown operator: " << opcode << std::endl;
#endif
            result = GPScalar::zero();
            break;
    }
    
    return applyRangeLimit(result);
}

fastmath::GPScalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg) {
    // Calculate the vector from craft to target in world frame
    int idx = getPathIndex(pathProvider, aircraftState, arg.toDouble());
    Eigen::Vector3d craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    
    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    
    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());
    
    // Calculate the angle between the projected vector and the body Z-axis
    return fastmath::atan2Approx(fastmath::GPScalar::fromDouble(projectedVector.y()),
                                 fastmath::GPScalar::fromDouble(-projectedVector.z()));
}

fastmath::GPScalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg) {
    // Calculate the vector from craft to target in world frame
    int idx = getPathIndex(pathProvider, aircraftState, arg.toDouble());
    Eigen::Vector3d craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    
    // Transform the craft-to-target vector to body frame
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    
    // Project the craft-to-target vector onto the body YZ plane
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());
    
    // Calculate the angle between the projected vector and the body Z-axis
    double rollEstimate = fastmath::atan2Approx(
        fastmath::GPScalar::fromDouble(projectedVector.y()),
        fastmath::GPScalar::fromDouble(-projectedVector.z())).toDouble();
    
    // *** PITCH: Calculate the vector from craft to target in world frame if it did rotate
    Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
    Eigen::Quaterniond virtualOrientation = aircraftState.getOrientation() * rollRotation;
    
    // Transform target vector to new virtual orientation
    Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;
    
    // Calculate pitch angle
    return fastmath::atan2Approx(
        fastmath::GPScalar::fromDouble(-newLocalTargetVector.z()),
        fastmath::GPScalar::fromDouble(newLocalTargetVector.x()));
}

fastmath::GPScalar executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg) {
    int idx = getPathIndex(pathProvider, aircraftState, arg.toDouble());
    double distance = (pathProvider.getPath(idx).start - aircraftState.getPosition()).norm();
    double relVel = aircraftState.getRelVelScalar().toDouble();
    double clamped = CLAMP_DEF((distance - 10) / relVel, -1.0, 1.0);
    return fastmath::GPScalar::fromDouble(clamped);
}

fastmath::GPScalar executeGetDHome(AircraftState& aircraftState) {
    double distance = (Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE) - aircraftState.getPosition()).norm();
    return fastmath::GPScalar::fromDouble(distance);
}

// Portable bytecode evaluation implementation - works on all platforms
fastmath::GPScalar evaluateBytecodePortable(const struct GPBytecode* program, int program_size, 
                               PathProvider& pathProvider, AircraftState& aircraftState, 
                               fastmath::GPScalar contextArg) {
    const int MAX_STACK_SIZE = 256;
    fastmath::GPScalar stack[MAX_STACK_SIZE];
    int stack_ptr = 0;

    auto pushResult = [&](fastmath::GPScalar value) {
        stack[stack_ptr++] = applyRangeLimit(value);
    };

    auto popScalar = [&]() -> fastmath::GPScalar {
        return stack[--stack_ptr];
    };
    
    // Execute bytecode instructions
    for (int i = 0; i < program_size; i++) {
        const struct GPBytecode& instruction = program[i];
        switch (instruction.opcode) {
            case ADD:
            case SUB:
            case MUL:
            case DIV:
            case EQ:
            case GT:
            case ATAN2:
            case MIN:
            case MAX: {
                if (stack_ptr < 2) return fastmath::GPScalar::zero();
                fastmath::GPScalar rhs = popScalar();
                fastmath::GPScalar lhs = popScalar();
                fastmath::GPScalar value = fastmath::GPScalar::zero();
                switch (instruction.opcode) {
                    case ADD: value = lhs + rhs; break;
                    case SUB: value = lhs - rhs; break;
                    case MUL: value = lhs * rhs; break;
                    case DIV: value = rhs.isZero() ? fastmath::GPScalar::zero() : lhs / rhs; break;
                    case EQ:  value = (lhs == rhs) ? fastmath::GPScalar::fromInt(1) : fastmath::GPScalar::zero(); break;
                    case GT:  value = (lhs > rhs) ? fastmath::GPScalar::fromInt(1) : fastmath::GPScalar::zero(); break;
                    case ATAN2: value = fastmath::atan2Approx(lhs, rhs); break;
                    case MIN: value = fastmath::min(lhs, rhs); break;
                    case MAX: value = fastmath::max(lhs, rhs); break;
                    default: break;
                }
                pushResult(value);
                break;
            }

            case SIN:
            case COS:
            case ABS:
            case SQRT:
            case SETPITCH:
            case SETROLL:
            case SETTHROTTLE:
            case GETDPHI:
            case GETDTHETA:
            case GETDTARGET: {
                if (stack_ptr < 1) return fastmath::GPScalar::zero();
                fastmath::GPScalar arg0 = popScalar();
                fastmath::GPScalar value = fastmath::GPScalar::zero();
                switch (instruction.opcode) {
                    case SIN: value = fastmath::sinApprox(arg0); break;
                    case COS: value = fastmath::cosApprox(arg0); break;
                    case ABS: value = arg0.abs(); break;
                    case SQRT: value = fastmath::sqrtApprox(arg0); break;
                    case SETPITCH:
                        value = fastmath::GPScalar::fromDouble(aircraftState.setPitchCommand(arg0.toDouble()));
                        break;
                    case SETROLL:
                        value = fastmath::GPScalar::fromDouble(aircraftState.setRollCommand(arg0.toDouble()));
                        break;
                    case SETTHROTTLE:
                        value = fastmath::GPScalar::fromDouble(aircraftState.setThrottleCommand(arg0.toDouble()));
                        break;
                    case GETDPHI:
                        value = executeGetDPhi(pathProvider, aircraftState, arg0);
                        break;
                    case GETDTHETA:
                        value = executeGetDTheta(pathProvider, aircraftState, arg0);
                        break;
                    case GETDTARGET:
                        value = executeGetDTarget(pathProvider, aircraftState, arg0);
                        break;
                    default:
                        break;
                }
                pushResult(value);
                break;
            }

            case IF:
            case CLAMP: {
                if (stack_ptr < 3) return fastmath::GPScalar::zero();
                fastmath::GPScalar c = popScalar();
                fastmath::GPScalar b = popScalar();
                fastmath::GPScalar a = popScalar();
                fastmath::GPScalar value = fastmath::GPScalar::zero();
                if (instruction.opcode == IF) {
                    value = (!a.isZero()) ? b : c;
                } else {
                    value = fastmath::clamp(a, b, c);
                }
                pushResult(value);
                break;
            }
            case PROGN: {
                if (stack_ptr < 2) return fastmath::GPScalar::zero();
                fastmath::GPScalar last = popScalar();
                popScalar(); // discard first
                pushResult(last);
                break;
            }

            case GETPITCH:
                pushResult(aircraftState.getPitchCommandScalar());
                break;
            case GETROLL:
                pushResult(aircraftState.getRollCommandScalar());
                break;
            case GETTHROTTLE:
                pushResult(aircraftState.getThrottleCommandScalar());
                break;
            case GETVEL:
                pushResult(aircraftState.getRelVelScalar());
                break;
            case GETVELX:
                pushResult(fastmath::GPScalar::fromDouble(aircraftState.getVelocity().x()));
                break;
            case GETVELY:
                pushResult(fastmath::GPScalar::fromDouble(aircraftState.getVelocity().y()));
                break;
            case GETVELZ:
                pushResult(fastmath::GPScalar::fromDouble(aircraftState.getVelocity().z()));
                break;
            case GETDHOME:
                pushResult(executeGetDHome(aircraftState));
                break;
            case GETALPHA: {
                Eigen::Vector3d velocity_body = fastmath::rotateWorldToBody(
                    aircraftState.getOrientation(), aircraftState.getVelocity());
                pushResult(fastmath::atan2Approx(
                    fastmath::GPScalar::fromDouble(-velocity_body.z()),
                    fastmath::GPScalar::fromDouble(velocity_body.x())));
                break;
            }
            case GETBETA: {
                Eigen::Vector3d velocity_body = fastmath::rotateWorldToBody(
                    aircraftState.getOrientation(), aircraftState.getVelocity());
                pushResult(fastmath::atan2Approx(
                    fastmath::GPScalar::fromDouble(velocity_body.y()),
                    fastmath::GPScalar::fromDouble(velocity_body.x())));
                break;
            }
            case GETROLL_RAD: {
                double roll = fastmath::rollFromQuaternion(aircraftState.getOrientation());
                pushResult(fastmath::GPScalar::fromDouble(roll));
                break;
            }
            case GETPITCH_RAD: {
                double pitch = fastmath::pitchFromQuaternion(aircraftState.getOrientation());
                pushResult(fastmath::GPScalar::fromDouble(pitch));
                break;
            }
            case OP_PI:
                pushResult(fastmath::GPScalar::fromDouble(M_PI));
                break;
            case ZERO:
                pushResult(fastmath::GPScalar::zero());
                break;
            case ONE:
                pushResult(fastmath::GPScalar::fromInt(1));
                break;
            case TWO:
                pushResult(fastmath::GPScalar::fromInt(2));
                break;

            default: {
                fastmath::GPScalar value = evaluateGPOperator(instruction.opcode, pathProvider, aircraftState, nullptr, 0, contextArg);
                pushResult(value);
                break;
            }
        }
        
        if (stack_ptr >= MAX_STACK_SIZE) {
#ifdef GP_BUILD
            std::cerr << "Error: Stack overflow in bytecode execution" << std::endl;
#endif
            return fastmath::GPScalar::zero();
        }
    }

    if (stack_ptr != 1) {
#ifdef GP_BUILD
        std::cerr << "Error: Invalid stack state after bytecode execution (stack_ptr=" << stack_ptr << ")" << std::endl;
#endif
        return fastmath::GPScalar::zero();
    }

    return applyRangeLimit(stack[0]);
}

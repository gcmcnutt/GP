#include "gp_evaluator_portable.h"
#include <array>
#include <cmath>
#include <type_traits>
#include "fastmath/fixed_math.h"

#ifdef GP_BUILD
#include "gp_bytecode.h"  // Full GPBytecode definition for GP builds
#endif

#ifdef GP_BUILD
#include <iostream>
#endif

namespace {

constexpr double kRangeLimit = 1000000.0;

template <typename Scalar>
struct ScalarOps;

template <>
struct ScalarOps<fastmath::GPScalar> {
    static fastmath::GPScalar zero() { return fastmath::GPScalar::zero(); }
    static fastmath::GPScalar fromDouble(double value) { return fastmath::GPScalar::fromDouble(value); }
    static double toDouble(const fastmath::GPScalar& value) { return value.toDouble(); }
    static bool isZero(const fastmath::GPScalar& value) { return value.isZero(); }
    static fastmath::GPScalar sin(const fastmath::GPScalar& value) { return fastmath::sinApprox(value); }
    static fastmath::GPScalar cos(const fastmath::GPScalar& value) { return fastmath::cosApprox(value); }
    static fastmath::GPScalar atan2(const fastmath::GPScalar& y, const fastmath::GPScalar& x) { return fastmath::atan2Approx(y, x); }
    static fastmath::GPScalar sqrt(const fastmath::GPScalar& value) { return fastmath::sqrtApprox(value); }
    static fastmath::GPScalar clamp(const fastmath::GPScalar& value,
                                    const fastmath::GPScalar& min,
                                    const fastmath::GPScalar& max) {
        return fastmath::clamp(value, min, max);
    }
    static fastmath::GPScalar abs(const fastmath::GPScalar& value) { return value.abs(); }
    static fastmath::GPScalar min(const fastmath::GPScalar& a, const fastmath::GPScalar& b) {
        return fastmath::min(a, b);
    }
    static fastmath::GPScalar max(const fastmath::GPScalar& a, const fastmath::GPScalar& b) {
        return fastmath::max(a, b);
    }
    static fastmath::GPScalar applyRangeLimit(const fastmath::GPScalar& value) {
        return ::applyRangeLimit(value);
    }
};

template <>
struct ScalarOps<double> {
    static double zero() { return 0.0; }
    static double fromDouble(double value) { return value; }
    static double toDouble(double value) { return value; }
    static bool isZero(double value) { return std::abs(value) < 1e-12; }
    static double sin(double value) { return std::sin(value); }
    static double cos(double value) { return std::cos(value); }
    static double atan2(double y, double x) { return std::atan2(y, x); }
    static double sqrt(double value) { return value <= 0.0 ? 0.0 : std::sqrt(value); }
    static double clamp(double value, double min, double max) { return CLAMP_DEF(value, min, max); }
    static double abs(double value) { return std::abs(value); }
    static double min(double a, double b) { return std::min(a, b); }
    static double max(double a, double b) { return std::max(a, b); }
    static double applyRangeLimit(double value) {
        if (value < -kRangeLimit) {
            return -kRangeLimit;
        }
        if (value > kRangeLimit) {
            return kRangeLimit;
        }
        if (std::abs(value) < 1e-6) {
            return 0.0;
        }
        return value;
    }
};

template <typename Scalar>
inline Scalar getPitchCommandValue(AircraftState& state) {
    if constexpr (std::is_same_v<Scalar, fastmath::GPScalar>) {
        return state.getPitchCommandScalar();
    } else {
        return ScalarOps<Scalar>::fromDouble(state.getPitchCommand());
    }
}

template <typename Scalar>
inline Scalar getRollCommandValue(AircraftState& state) {
    if constexpr (std::is_same_v<Scalar, fastmath::GPScalar>) {
        return state.getRollCommandScalar();
    } else {
        return ScalarOps<Scalar>::fromDouble(state.getRollCommand());
    }
}

template <typename Scalar>
inline Scalar getThrottleCommandValue(AircraftState& state) {
    if constexpr (std::is_same_v<Scalar, fastmath::GPScalar>) {
        return state.getThrottleCommandScalar();
    } else {
        return ScalarOps<Scalar>::fromDouble(state.getThrottleCommand());
    }
}

template <typename Scalar>
inline Scalar getRelVelValue(AircraftState& state) {
    if constexpr (std::is_same_v<Scalar, fastmath::GPScalar>) {
        return state.getRelVelScalar();
    } else {
        return ScalarOps<Scalar>::fromDouble(state.getRelVel());
    }
}

template <typename Scalar>
Scalar executeGetDPhiImpl(PathProvider& pathProvider, AircraftState& aircraftState, Scalar arg) {
    int idx = getPathIndex(pathProvider, aircraftState, ScalarOps<Scalar>::toDouble(arg));
    Eigen::Vector3d craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());
    return ScalarOps<Scalar>::atan2(
        ScalarOps<Scalar>::fromDouble(projectedVector.y()),
        ScalarOps<Scalar>::fromDouble(-projectedVector.z()));
}

template <typename Scalar>
Scalar executeGetDThetaImpl(PathProvider& pathProvider, AircraftState& aircraftState, Scalar arg) {
    int idx = getPathIndex(pathProvider, aircraftState, ScalarOps<Scalar>::toDouble(arg));
    Eigen::Vector3d craftToTarget = pathProvider.getPath(idx).start - aircraftState.getPosition();
    Eigen::Vector3d target_local = aircraftState.getOrientation().inverse() * craftToTarget;
    Eigen::Vector3d projectedVector(0, target_local.y(), target_local.z());

    double rollEstimate = ScalarOps<Scalar>::toDouble(
        ScalarOps<Scalar>::atan2(
            ScalarOps<Scalar>::fromDouble(projectedVector.y()),
            ScalarOps<Scalar>::fromDouble(-projectedVector.z())));

    Eigen::Quaterniond rollRotation(Eigen::AngleAxisd(rollEstimate, Eigen::Vector3d::UnitX()));
    Eigen::Quaterniond virtualOrientation = aircraftState.getOrientation() * rollRotation;
    Eigen::Vector3d newLocalTargetVector = virtualOrientation.inverse() * craftToTarget;

    return ScalarOps<Scalar>::atan2(
        ScalarOps<Scalar>::fromDouble(-newLocalTargetVector.z()),
        ScalarOps<Scalar>::fromDouble(newLocalTargetVector.x()));
}

template <typename Scalar>
Scalar executeGetDTargetImpl(PathProvider& pathProvider, AircraftState& aircraftState, Scalar arg) {
    int idx = getPathIndex(pathProvider, aircraftState, ScalarOps<Scalar>::toDouble(arg));
    double distance = (pathProvider.getPath(idx).start - aircraftState.getPosition()).norm();
    double relVel = aircraftState.getRelVel();
    double clamped = CLAMP_DEF((distance - 10.0) / relVel, -1.0, 1.0);
    return ScalarOps<Scalar>::fromDouble(clamped);
}

template <typename Scalar>
Scalar executeGetDHomeImpl(AircraftState& aircraftState) {
    double distance = (Eigen::Vector3d(0, 0, SIM_INITIAL_ALTITUDE) - aircraftState.getPosition()).norm();
    return ScalarOps<Scalar>::fromDouble(distance);
}

template <typename Scalar>
Scalar evaluateGPOperatorImpl(int opcode, PathProvider& pathProvider,
                              AircraftState& aircraftState,
                              const Scalar* args, int argc, Scalar contextArg) {
    Scalar result = ScalarOps<Scalar>::zero();
    switch (opcode) {
        case ADD:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = args[0] + args[1];
            break;
        case SUB:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = args[0] - args[1];
            break;
        case MUL:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = args[0] * args[1];
            break;
        case DIV:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::isZero(args[1]) ? ScalarOps<Scalar>::zero() : args[0] / args[1];
            break;
        case SETPITCH:
            if (!args || argc < 1) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::fromDouble(
                aircraftState.setPitchCommand(ScalarOps<Scalar>::toDouble(args[0])));
            break;
        case SETROLL:
            if (!args || argc < 1) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::fromDouble(
                aircraftState.setRollCommand(ScalarOps<Scalar>::toDouble(args[0])));
            break;
        case SETTHROTTLE:
            if (!args || argc < 1) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::fromDouble(
                aircraftState.setThrottleCommand(ScalarOps<Scalar>::toDouble(args[0])));
            break;
        case GETPITCH:
            result = getPitchCommandValue<Scalar>(aircraftState);
            break;
        case GETROLL:
            result = getRollCommandValue<Scalar>(aircraftState);
            break;
        case GETTHROTTLE:
            result = getThrottleCommandValue<Scalar>(aircraftState);
            break;
        case GETVEL:
            result = getRelVelValue<Scalar>(aircraftState);
            break;
        case GETDPHI:
            result = executeGetDPhiImpl(
                pathProvider, aircraftState,
                (args && argc > 0) ? args[0] : contextArg);
            break;
        case GETDTHETA:
            result = executeGetDThetaImpl(
                pathProvider, aircraftState,
                (args && argc > 0) ? args[0] : contextArg);
            break;
        case GETDTARGET:
            result = executeGetDTargetImpl(
                pathProvider, aircraftState,
                (args && argc > 0) ? args[0] : contextArg);
            break;
        case GETDHOME:
            result = executeGetDHomeImpl<Scalar>(aircraftState);
            break;
        case SIN:
            if (!args || argc < 1) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::sin(args[0]);
            break;
        case COS:
            if (!args || argc < 1) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::cos(args[0]);
            break;
        case ATAN2:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::atan2(args[0], args[1]);
            break;
        case CLAMP:
            if (!args || argc < 3) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::clamp(args[0], args[1], args[2]);
            break;
        case ABS:
            if (!args || argc < 1) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::abs(args[0]);
            break;
        case SQRT:
            if (!args || argc < 1) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::sqrt(args[0]);
            break;
        case MIN:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::min(args[0], args[1]);
            break;
        case MAX:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::max(args[0], args[1]);
            break;
        case IF:
            if (!args || argc < 3) return ScalarOps<Scalar>::zero();
            result = ScalarOps<Scalar>::isZero(args[0]) ? args[2] : args[1];
            break;
        case EQ:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = (args[0] == args[1]) ? ScalarOps<Scalar>::fromDouble(1.0) : ScalarOps<Scalar>::zero();
            break;
        case GT:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = (args[0] > args[1]) ? ScalarOps<Scalar>::fromDouble(1.0) : ScalarOps<Scalar>::zero();
            break;
        case PROGN:
            if (!args || argc < 2) return ScalarOps<Scalar>::zero();
            result = args[1];
            break;
        case OP_PI:
            result = ScalarOps<Scalar>::fromDouble(M_PI);
            break;
        case ZERO:
            result = ScalarOps<Scalar>::zero();
            break;
        case ONE:
            result = ScalarOps<Scalar>::fromDouble(1.0);
            break;
        case TWO:
            result = ScalarOps<Scalar>::fromDouble(2.0);
            break;
        case GETVELX:
            result = ScalarOps<Scalar>::fromDouble(aircraftState.getVelocity().x());
            break;
        case GETVELY:
            result = ScalarOps<Scalar>::fromDouble(aircraftState.getVelocity().y());
            break;
        case GETVELZ:
            result = ScalarOps<Scalar>::fromDouble(aircraftState.getVelocity().z());
            break;
        case GETALPHA: {
            Eigen::Vector3d velocity_body = fastmath::rotateWorldToBody(
                aircraftState.getOrientation(), aircraftState.getVelocity());
            result = ScalarOps<Scalar>::atan2(
                ScalarOps<Scalar>::fromDouble(-velocity_body.z()),
                ScalarOps<Scalar>::fromDouble(velocity_body.x()));
            break;
        }
        case GETBETA: {
            Eigen::Vector3d velocity_body = fastmath::rotateWorldToBody(
                aircraftState.getOrientation(), aircraftState.getVelocity());
            result = ScalarOps<Scalar>::atan2(
                ScalarOps<Scalar>::fromDouble(velocity_body.y()),
                ScalarOps<Scalar>::fromDouble(velocity_body.x()));
            break;
        }
        case GETROLL_RAD: {
            double roll = fastmath::rollFromQuaternion(aircraftState.getOrientation());
            result = ScalarOps<Scalar>::fromDouble(roll);
            break;
        }
        case GETPITCH_RAD: {
            double pitch = fastmath::pitchFromQuaternion(aircraftState.getOrientation());
            result = ScalarOps<Scalar>::fromDouble(pitch);
            break;
        }
        default:
#ifdef GP_BUILD
            std::cerr << "Unknown operator: " << opcode << std::endl;
#endif
            result = ScalarOps<Scalar>::zero();
            break;
    }
    return ScalarOps<Scalar>::applyRangeLimit(result);
}

template <typename Scalar>
Scalar evaluateBytecodeImpl(const GPBytecode* program, int program_size,
                            PathProvider& pathProvider, AircraftState& aircraftState,
                            Scalar contextArg) {
    constexpr int MAX_STACK_SIZE = 256;
    std::array<Scalar, MAX_STACK_SIZE> stack{};
    int stack_ptr = 0;

    auto pushResult = [&](Scalar value) {
        stack[stack_ptr++] = ScalarOps<Scalar>::applyRangeLimit(value);
    };

    auto popScalar = [&]() -> Scalar {
        return stack[--stack_ptr];
    };

    for (int i = 0; i < program_size; ++i) {
        const GPBytecode& instruction = program[i];
        switch (instruction.opcode) {
            case ADD:
            case SUB:
            case MUL:
            case DIV:
            case EQ:
            case GT: {
                if (stack_ptr < 2) {
                    return ScalarOps<Scalar>::zero();
                }
                Scalar rhs = popScalar();
                Scalar lhs = popScalar();
                Scalar value = ScalarOps<Scalar>::zero();
                switch (instruction.opcode) {
                    case ADD: value = lhs + rhs; break;
                    case SUB: value = lhs - rhs; break;
                    case MUL: value = lhs * rhs; break;
                    case DIV: value = ScalarOps<Scalar>::isZero(rhs) ? ScalarOps<Scalar>::zero() : lhs / rhs; break;
                    case EQ: value = (lhs == rhs) ? ScalarOps<Scalar>::fromDouble(1.0) : ScalarOps<Scalar>::zero(); break;
                    case GT: value = (lhs > rhs) ? ScalarOps<Scalar>::fromDouble(1.0) : ScalarOps<Scalar>::zero(); break;
                }
                pushResult(value);
                break;
            }
            case SIN:
            case COS:
            case ABS:
        case SQRT:
        case GETPITCH:
        case GETROLL:
        case GETTHROTTLE:
        case GETVEL:
        case SETPITCH:
        case SETROLL:
        case SETTHROTTLE:
        case GETDPHI:
        case GETDTHETA:
        case GETDTARGET:
        case GETDHOME:
        case GETVELX:
        case GETVELY:
        case GETVELZ:
        case GETALPHA:
        case GETBETA:
        case GETROLL_RAD:
        case GETPITCH_RAD: {
                Scalar localArgs[3];
                Scalar* argPtr = nullptr;
                if (instruction.argc > 0) {
                    if (stack_ptr < instruction.argc) {
                        return ScalarOps<Scalar>::zero();
                    }
                    for (int argIndex = instruction.argc - 1; argIndex >= 0; --argIndex) {
                        localArgs[argIndex] = popScalar();
                    }
                    argPtr = localArgs;
                }
                Scalar value = evaluateGPOperatorImpl<Scalar>(
                    instruction.opcode, pathProvider, aircraftState, argPtr, instruction.argc, contextArg);
                pushResult(value);
                break;
            }
            case CLAMP:
            case MIN:
            case MAX:
            case ATAN2: {
                if (stack_ptr < 2) {
                    return ScalarOps<Scalar>::zero();
                }
                Scalar arg2 = popScalar();
                Scalar arg1 = popScalar();
                Scalar value = ScalarOps<Scalar>::zero();
                switch (instruction.opcode) {
                    case CLAMP: {
                        if (stack_ptr < 1) {
                            return ScalarOps<Scalar>::zero();
                        }
                        Scalar arg0 = popScalar();
                        value = ScalarOps<Scalar>::clamp(arg0, arg1, arg2);
                        break;
                    }
                    case MIN:
                        value = ScalarOps<Scalar>::min(arg1, arg2);
                        break;
                    case MAX:
                        value = ScalarOps<Scalar>::max(arg1, arg2);
                        break;
                    case ATAN2:
                        value = ScalarOps<Scalar>::atan2(arg1, arg2);
                        break;
                }
                pushResult(value);
                break;
            }
            case IF: {
                if (stack_ptr < 3) {
                    return ScalarOps<Scalar>::zero();
                }
                Scalar falseBranch = popScalar();
                Scalar trueBranch = popScalar();
                Scalar condition = popScalar();
                pushResult(ScalarOps<Scalar>::isZero(condition) ? falseBranch : trueBranch);
                break;
            }
            case PROGN: {
                if (stack_ptr < instruction.argc) {
                    return ScalarOps<Scalar>::zero();
                }
                Scalar value = popScalar();
                stack_ptr -= instruction.argc - 1;
                pushResult(value);
                break;
            }
            case OP_PI:
            case ZERO:
            case ONE:
            case TWO: {
                double literal = 0.0;
                switch (instruction.opcode) {
                    case OP_PI: literal = M_PI; break;
                    case ZERO: literal = 0.0; break;
                    case ONE: literal = 1.0; break;
                    case TWO: literal = 2.0; break;
                }
                pushResult(ScalarOps<Scalar>::fromDouble(literal));
                break;
            }
            default: {
                Scalar value = evaluateGPOperatorImpl<Scalar>(
                    instruction.opcode, pathProvider, aircraftState, nullptr, 0, contextArg);
                pushResult(value);
                break;
            }
        }
    }

    if (stack_ptr == 0) {
        return ScalarOps<Scalar>::zero();
    }
    return stack[stack_ptr - 1];
}

} // namespace


fastmath::GPScalar evaluateGPOperator(int opcode, PathProvider& pathProvider,
                         AircraftState& aircraftState,
                         const fastmath::GPScalar* args, int argc, fastmath::GPScalar contextArg) {
    return evaluateGPOperatorImpl<fastmath::GPScalar>(opcode, pathProvider, aircraftState, args, argc, contextArg);
}

double evaluateGPOperatorReference(int opcode, PathProvider& pathProvider,
                                   AircraftState& aircraftState,
                                   const double* args, int argc, double contextArg) {
    return evaluateGPOperatorImpl<double>(opcode, pathProvider, aircraftState, args, argc, contextArg);
}

fastmath::GPScalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg) {
    return executeGetDPhiImpl(pathProvider, aircraftState, arg);
}

fastmath::GPScalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg) {
    return executeGetDThetaImpl(pathProvider, aircraftState, arg);
}

fastmath::GPScalar executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg) {
    return executeGetDTargetImpl(pathProvider, aircraftState, arg);
}

fastmath::GPScalar executeGetDHome(AircraftState& aircraftState) {
    return executeGetDHomeImpl<fastmath::GPScalar>(aircraftState);
}

fastmath::GPScalar evaluateBytecodePortable(const struct GPBytecode* program, int program_size,
                               PathProvider& pathProvider, AircraftState& aircraftState,
                               fastmath::GPScalar contextArg) {
    return evaluateBytecodeImpl<fastmath::GPScalar>(program, program_size, pathProvider, aircraftState, contextArg);
}

double evaluateBytecodeReference(const struct GPBytecode* program, int program_size,
                                 PathProvider& pathProvider, AircraftState& aircraftState,
                                 double contextArg) {
    return evaluateBytecodeImpl<double>(program, program_size, pathProvider, aircraftState, contextArg);
}

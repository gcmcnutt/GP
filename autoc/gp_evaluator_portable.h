#ifndef GP_EVALUATOR_PORTABLE_H
#define GP_EVALUATOR_PORTABLE_H

#include "aircraft_state.h"
#include "fastmath/gp_scalar.h"

#ifdef GP_BUILD
#include "autoc.h"
// Forward declaration - actual definition depends on build context
struct GPBytecode;
#else
// Embedded build: define operators enum locally
enum Operators {
  ADD = 0, SUB, MUL, DIV,
  IF, EQ, GT,
  SIN, COS,
  GETDPHI, GETDTHETA, GETDTARGET, GETDHOME, GETVEL,
  GETPITCH, GETROLL, GETTHROTTLE,
  SETPITCH, SETROLL, SETTHROTTLE,
  GETALPHA, GETBETA, GETVELX, GETVELY, GETVELZ,
  GETROLL_RAD, GETPITCH_RAD,
  CLAMP, ATAN2, ABS, SQRT, MIN, MAX,
  OP_PI, ZERO, ONE, TWO, PROGN, _END  // Renamed PI to OP_PI to avoid Arduino macro conflict
};

// Simple bytecode structure for embedded
struct GPBytecode {
    uint8_t opcode;     // Operation code (maps to Operators enum)
    uint8_t argc;       // Number of arguments
    float constant;     // For literal values (PI, 0, 1, 2)
    
    GPBytecode(uint8_t op = 0, uint8_t args = 0, float val = 0.0f) 
        : opcode(op), argc(args), constant(val) {}
};
#endif

// Single portable evaluation function that works everywhere
fastmath::GPScalar evaluateGPOperator(int opcode, PathProvider& pathProvider, 
                         AircraftState& aircraftState, 
                         const fastmath::GPScalar* args, int argc, fastmath::GPScalar contextArg = fastmath::GPScalar::zero());

// Navigation helpers - same logic, different path access
fastmath::GPScalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg);
fastmath::GPScalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg);
fastmath::GPScalar executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, fastmath::GPScalar arg);
fastmath::GPScalar executeGetDHome(AircraftState& aircraftState);

// Range limiting - identical across platforms
inline fastmath::GPScalar applyRangeLimit(fastmath::GPScalar value) {
    const double RANGELIMIT = 1000000.0;
    double asDouble = value.toDouble();
    if (asDouble < -RANGELIMIT) {
        fastmath::recordRangeClamp();
        return fastmath::GPScalar::fromDouble(-RANGELIMIT);
    }
    if (asDouble > RANGELIMIT) {
        fastmath::recordRangeClamp();
        return fastmath::GPScalar::fromDouble(RANGELIMIT);
    }
    if (value.abs() < fastmath::GPScalar::fromDouble(0.000001)) {
        fastmath::recordZeroSnap();
        return fastmath::GPScalar::zero();
    }
    return value;
}

// Stack-based bytecode evaluation using portable operators
fastmath::GPScalar evaluateBytecodePortable(const GPBytecode* program, int program_size, 
                               PathProvider& pathProvider, AircraftState& aircraftState, 
                               fastmath::GPScalar contextArg = fastmath::GPScalar::zero());

#endif

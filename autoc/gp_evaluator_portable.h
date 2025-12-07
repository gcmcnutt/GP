#ifndef GP_EVALUATOR_PORTABLE_H
#define GP_EVALUATOR_PORTABLE_H

#include <cstdint>
#include "aircraft_state.h"

#if defined(GP_BUILD) && !defined(GP_TEST)
#include "autoc.h"
// Forward declaration - actual definition depends on build context
struct GPBytecode;
#endif

// Local operators/bytecode definitions for test/embedded builds
#if !defined(GP_BUILD) || defined(GP_TEST)
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

struct GPBytecode {
    uint8_t opcode;
    uint8_t argc;
    float constant;
    
    GPBytecode(uint8_t op = 0, uint8_t args = 0, float val = 0.0f) 
        : opcode(op), argc(args), constant(val) {}
};
#endif

// Single portable evaluation function that works everywhere
gp_scalar evaluateGPOperator(int opcode, PathProvider& pathProvider, 
                         AircraftState& aircraftState, 
                         const gp_scalar* args, int argc, gp_scalar contextArg = 0.0f);

// Navigation helpers - same logic, different path access
gp_scalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDHome(AircraftState& aircraftState);

// Range limiting - identical across platforms
inline gp_scalar applyRangeLimit(gp_scalar value) {
    const gp_scalar RANGELIMIT = static_cast<gp_scalar>(1000000.0f);
    if (value < -RANGELIMIT) return -RANGELIMIT;
    if (value > RANGELIMIT) return RANGELIMIT;
    if (ABS_DEF(value) < static_cast<gp_scalar>(0.000001f)) return 0.0f;
    return value;
}

// Stack-based bytecode evaluation using portable operators
gp_scalar evaluateBytecodePortable(const GPBytecode* program, int program_size, 
                               PathProvider& pathProvider, AircraftState& aircraftState, 
                               gp_scalar contextArg = 0.0f);

#endif

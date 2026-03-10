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
// NOTE: Must stay in sync with autoc.h - new operators go at END before _END
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
  OP_PI, ZERO, ONE, TWO, PROGN,
  // Temporal state terminals (added 2026-02)
  GETDPHI_PREV, GETDTHETA_PREV,
  GETDPHI_RATE, GETDTHETA_RATE,
  _END  // Renamed PI to OP_PI to avoid Arduino macro conflict
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

// Path interpolation - returns position at goal time (binary search + linear lerp)
// Uses int32_t timestamps for deterministic binary search (no float precision issues)
gp_vec3 getInterpolatedTargetPosition(PathProvider& pathProvider,
                                       int32_t currentTimeMsec,
                                       gp_scalar offsetSteps);

// Navigation helpers - same logic, different path access
gp_scalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDHome(AircraftState& aircraftState);

// Temporal terminals - PREV takes history index, RATE is nullary
gp_scalar executeGetDPhiPrev(AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDThetaPrev(AircraftState& aircraftState, gp_scalar arg);
gp_scalar executeGetDPhiRate(AircraftState& aircraftState);
gp_scalar executeGetDThetaRate(AircraftState& aircraftState);

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

// Expose LUT functions for bit-accurate testing
#ifdef GP_TEST
gp_scalar testFastSin(gp_scalar angle);
gp_scalar testFastCos(gp_scalar angle);
gp_scalar testFastAtan2(gp_scalar y, gp_scalar x);
#endif

#endif

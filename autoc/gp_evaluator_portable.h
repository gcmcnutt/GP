#ifndef GP_EVALUATOR_PORTABLE_H
#define GP_EVALUATOR_PORTABLE_H

#include "aircraft_state.h"

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

// Create compatibility alias for PI
#ifndef GP_BUILD
#define PI OP_PI
#endif

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
double evaluateGPOperator(int opcode, PathProvider& pathProvider, 
                         AircraftState& aircraftState, 
                         const double* args, int argc, double contextArg = 0.0);

// Navigation helpers - same logic, different path access
double executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, double arg);
double executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, double arg);
double executeGetDTarget(PathProvider& pathProvider, AircraftState& aircraftState, double arg);
double executeGetDHome(AircraftState& aircraftState);

// Range limiting - identical across platforms
inline double applyRangeLimit(double value) {
    const double RANGELIMIT = 1000000.0;
    if (value < -RANGELIMIT) return -RANGELIMIT;
    if (value > RANGELIMIT) return RANGELIMIT;
    if (ABS_DEF(value) < 0.000001) return 0.0;
    return value;
}

// Stack-based bytecode evaluation using portable operators
double evaluateBytecodePortable(const GPBytecode* program, int program_size, 
                               PathProvider& pathProvider, AircraftState& aircraftState, 
                               double contextArg = 0.0);

#endif
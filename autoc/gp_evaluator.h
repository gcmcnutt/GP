// GP Evaluator Interface for xiao-gp
#ifndef GP_EVALUATOR_H
#define GP_EVALUATOR_H

#ifdef GP_BUILD
#include <vector>
#include "aircraft_state.h"
#else
#include "GP/autoc/aircraft_state.h"
#endif

#ifdef GP_BUILD
// Main GP evaluation function with vector - for full GP build
double evaluateGP(AircraftState& aircraftState, const std::vector<Path>& path, double arg);
#endif

// Simplified GP evaluation function for embedded use - single path segment
double evaluateGPSimple(AircraftState& aircraftState, const Path& currentPath, double arg);

#endif
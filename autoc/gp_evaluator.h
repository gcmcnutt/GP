// GP Evaluator Interface for xiao-gp
#ifndef GP_EVALUATOR_H
#define GP_EVALUATOR_H

#include <vector>
#ifdef GP_BUILD
#include "aircraft_state.h"
#else
#include "GP/autoc/aircraft_state.h"
#endif

// Main GP evaluation function - implemented in generated code
double evaluateGP(AircraftState& aircraftState, const std::vector<Path>& path, double arg);

#endif
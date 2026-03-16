// sensor_math.h — Navigation sensor functions for temporal history recording
// Shared between minisim (history capture) and NN evaluator (reads history buffers)
// Extracted from gp_evaluator_portable during GP removal (014-nn-training-signal)

#ifndef SENSOR_MATH_H
#define SENSOR_MATH_H

#include <cstdint>
#include "aircraft_state.h"

// Path interpolation — returns position at goal time (binary search + linear lerp)
gp_vec3 getInterpolatedTargetPosition(PathProvider& pathProvider,
                                       int32_t currentTimeMsec,
                                       gp_scalar offsetSteps);

// Roll angle to target (body-frame YZ plane projection)
gp_scalar executeGetDPhi(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);

// Pitch angle to target (body-frame XZ plane)
gp_scalar executeGetDTheta(PathProvider& pathProvider, AircraftState& aircraftState, gp_scalar arg);

#endif

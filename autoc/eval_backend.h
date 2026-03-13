#ifndef EVAL_BACKEND_H
#define EVAL_BACKEND_H

#include "aircraft_state.h"

// Abstract interface for controller evaluation.
// All controller types (GP tree, bytecode, neural net) implement this.
class ControllerBackend {
public:
    virtual ~ControllerBackend() = default;

    // Run controller, set control commands on AircraftState
    virtual void evaluate(AircraftState& aircraftState, PathProvider& pathProvider) = 0;

    // Return controller type name for logging
    virtual const char* getName() const = 0;
};

#endif

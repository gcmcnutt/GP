#ifndef GP_EVALUATOR_CRCSIM_H
#define GP_EVALUATOR_CRCSIM_H

#include "gp_evaluator_portable.h"
#include <vector>

// CRRCSim-specific extensions
class GPEvaluatorCRRCSim {
public:
    // Evaluate using portable bytecode with CRRCSim path format
    static double evaluateBytecodeForCRRCSim(const GPBytecode* program, int program_size,
                                            std::vector<Path>& path, AircraftState& aircraftState, 
                                            double arg = 0.0);
    
    // Convert boost serialized bytecode to portable format
    // (This would be used when CRRCSim receives bytecode via RPC)
    static void convertFromBoostBytecode(const void* boostBytecode, int boost_size,
                                        std::vector<GPBytecode>& portableBytecode);
};

#endif
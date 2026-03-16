#ifndef EVAL_LOGGER_H
#define EVAL_LOGGER_H

#include "autoc/types.h"

// Shared logging extracted from autoc.cc and autoc-eval.cc.
// Writes per-step data to data.dat and per-generation stats to data.stc.
class EvalLogger {
public:
    void logStepHeader();
    void logStep();
    void logGenerationStats();
    void logBestController();
    void logNNWeightStats();
};

#endif

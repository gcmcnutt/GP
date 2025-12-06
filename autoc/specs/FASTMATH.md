# use much faster integer math for the GP evaluation

Today, the autoc GP evaluator uses a mix of float and double precision math for its calculations. This is done to balance performance and accuracy.  But recent training runs have resulted ig GP evaluations that take over 100msec per evaluation.  Which is causing an issue for the actual flight control responsiveness.

Is it possible to switch to using integer math for the GP evaluation?  This would be similar to how microcontrollers often use fixed point integer math for performance reasons.

## Proposal: FASTMATH
The proposal is to implement a FASTMATH mode for the GP evaluator that uses integer math instead of floating point math. This would involve:
1. Defining a fixed point representation for all real numbers used in the GP evaluation. For example, using 16.16 fixed point format (16 bits for integer part, 16 bits for fractional part).
2. Implementing integer versions of all mathematical operations used in the GP evaluation (addition, subtraction, multiplication, division, trigonometric functions, etc.) that operate on the fixed point representation.
3. Modifying the GP evaluator to use the integer math functions.

## Strategy
Example the code in autoc/gp_evaluator_portable.cc, its header and the consumers in
- ~/GP/autoc/minisim 
- ~/xiao-gp/src/msplink.cpp
- ~/crsim/crrcsim-0.9.13/src/mod_inputdev/inputdev_autoc.cpp

See how the values are injected and how they are consumed.  Examine the various match functions and propose limited precisio and perhaps fast lookup tables.

## Goal
We expect an evaluation loop to execute under 10msec, ideally under 5msec.  This would allow for real-time control of the aircraft using the GP evolved programs.

## Implementation Steps
1. Research fixed point math libraries or implement custom fixed point math functions.
2. Modify the GP evaluator to use fixed point math.
3. Test the modified GP evaluator for correctness and performance.

Important, since the actual system really only has 2 or 3 decimal digits of accuracy, a fixed point (either 16.16 or strict 16bit integer with scaling) should be optimal.


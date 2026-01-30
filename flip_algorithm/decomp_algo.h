#pragma once
#include "input.h"

namespace df {

// runs whole decomposition while(true) loop
// changes in.tri_current, in.step_history, etc.
// returns true iff current triangulation equals lower at the end
bool run_flip_algorithm(df::InputData& in, bool enable_vis);

} // namespace df
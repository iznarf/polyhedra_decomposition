#pragma once

#include "input.h"   // use existing StepRecord structure

namespace df {

void apply_step(const StepRecord& step, InputData& D);

} // namespace df

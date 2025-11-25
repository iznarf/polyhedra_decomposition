#pragma once
#include "input.h"   

namespace df {

void apply_step(const StepRecord& step, InputData& D);

void init_replay(InputData& D);
void replay_ui();

} // namespace df

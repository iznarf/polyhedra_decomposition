#pragma once
#include "input.h"

namespace df {


void apply_edge_flip(df::vertex_id ia,
                df::vertex_id ib,
                df::InputData& D);

void print_flip_history(const df::InputData& D);
} // namespace df
#pragma once
#include "input.h"
#include <array>

namespace df { namespace reg {


bool is_flip_conforming(df::vertex_id ia,
                        df::vertex_id ib,
                        const df::InputData& D);


}} // namespace df::reg

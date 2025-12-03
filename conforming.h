#pragma once
#include "input.h"
#include <array>

namespace df { namespace reg {


bool is_flip_conforming(df::vertex_id ia,
                        df::vertex_id ib,
                        const df::InputData& D, const df::Tri2& tri_current);


}} // namespace df::reg

#pragma once
#include "input.h"
#include <array>

namespace df { namespace reg {
    bool is_insertion_conforming(df::vertex_id id,
                                 const df::InputData& D);
}} // namespace df
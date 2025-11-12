#pragma once
#include "input.h"
#include <array>

namespace df { namespace reg {


bool is_flip_conforming(df::vertex_id ia,
                        df::vertex_id ib,
                        const df::InputData& D,
                        const std::unordered_map<df::vertex_id,int>& local_index_current,  const std::unordered_map<df::vertex_id,int>& local_index_lower);


}} // namespace df::reg

#pragma once
#include "input.h"
#include <vector>
#include <array>

namespace df {
namespace reg {


std::vector<std::array<vertex_id, 2>> find_locally_non_regular_edges(const Tri2& tri);
bool edge_in_target(const df::vertex_id ia, const df::vertex_id ib, const df::InputData& D);

} 
} 



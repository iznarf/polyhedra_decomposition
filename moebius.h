#pragma once

#include <vector>
#include <cstdint>

#include "poset.h"   // pst::Node
#include "poset2.h"  // pst2::Poset2

namespace mob {

// Generic Möbius on a poset represented by upward cover edges:
// cover_out[u] contains v where u < v is a cover relation.
std::int64_t mobius_xy_from_cover(
    const std::vector<std::vector<int>>& cover_out,
    int x, int y);

// Build "upward" cover adjacency for poset1 from pst::Node list.
// Your nodes[u].children are DOWN edges (parent/top -> child/bottom),
// so we reverse them to get bottom -> top (upward).
std::vector<std::vector<int>> build_cover_up_from_poset1_nodes(
    const std::vector<pst::Node>& nodes);

// ImGui UI: compute μ(x,y) for poset1 and poset2 and compare.
void draw_mobius_compare_ui(
    const std::vector<pst::Node>& poset1_nodes,
    const pst2::Poset2& poset2);

} // namespace mob


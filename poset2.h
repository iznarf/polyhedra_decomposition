#pragma once
#include "poset.h"
#include <vector>

namespace pst2 {




struct Poset2 {
    std::vector<std::vector<int>> out2;       // optional: full <=2 DAG
    std::vector<std::vector<int>> cover_out;  // covering edges u -> v where u < v and no w with u < w < v
};

Poset2 build_poset2(const df::InputData& D, const std::vector<pst::Node>& nodes, bool keep_full_out2 = false);
void print_cover_relations(const pst2::Poset2& P2);


} // namespace pst2

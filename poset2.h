#pragma once
#include "poset.h"
#include <vector>

namespace pst2 {

struct Poset2 {
    std::vector<std::vector<int>> out2;       // optional for keeping the full <=2 DAG
    std::vector<std::vector<int>> cover_out;  // covering edges u -> v where u < v and no w with u < w < v
};

Poset2 build_poset2(const df::InputData& D, const std::vector<pst::Node>& nodes, bool keep_full_out2 = false);
int print_cover_relations(const pst2::Poset2& P2);

std::vector<int> interval_xy(const Poset2& P, int x, int y);

std::vector<int> meet_candidates_xy(const Poset2& P, int x, int y); // maximal of [0,x] intersection [0,y]
std::vector<int> join_candidates_xy(const Poset2& P, int x, int y); // minimal of [x,1] intersection [y,1]



} // namespace pst2

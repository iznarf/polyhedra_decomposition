#pragma once
#include "dedekind_cut.h"
#include "poset_utils.h"

#include <vector>

namespace dm_completion {

struct DedekindPoset {
    std::vector<dm_completion::DedekindCut> cuts;

    // cover edges 
    std::vector<std::vector<int>> cover_up;
    std::vector<std::vector<int>> cover_down;

    // topological order of cuts
    std::vector<int> topo_up;
    std::vector<int> topo_down;
    // reachability bitsets: R[a][b] = 1 iff a <= b
    std::vector<pst::Bitset> reachability_up;
    std::vector<pst::Bitset> reachability_down;

    // cached: I' as bitsets over base poset
    std::vector<pst::Bitset> Ip_bit;

    // levels for drawing (longest chain from minima)
    std::vector<int> levels;
    int maxLevel = 0;
};


DedekindPoset build_dedekind_poset(const PosetView& P, std::vector<dm_completion::DedekindCut> cuts);

bool check_complete_lattice(const PosetView& P, bool verbose);


} // namespace dm_completion


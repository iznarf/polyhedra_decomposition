#pragma once
#include "poset2.h"
#include "dedekind_cut.h"
#include "poset_utils.h"

#include <vector>

namespace pst2 {

struct DedekindPoset {
    std::vector<DedekindCut> cuts;

    // cover edges 
    std::vector<std::vector<int>> cover_up;
    std::vector<std::vector<int>> cover_dn;

    // topological order of cuts
    std::vector<int> topo;
    // reachability bitsets: R[a][b] = 1 iff a <= b
    std::vector<pst::Bitset> R;

    // cached: I' as bitsets over base poset
    std::vector<pst::Bitset> Ip_bit;

    // levels for drawing (longest chain from minima)
    std::vector<int> level;
    int maxLevel = 0;
};



DedekindPoset build_dedekind_poset(const Poset2& P2, std::vector<DedekindCut> cuts);

bool check_complete_lattice(const pst2::DedekindPoset& D, bool verbose = true);

} // namespace pst2


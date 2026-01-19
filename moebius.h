#pragma once

#include <vector>
#include <cstdint>

#include "poset.h"
#include "poset2.h"
#include "poset_utils.h"

namespace mob {

using pst::Bitset;

// Cached “view” for a poset given by cover_up edges.
struct CoverUpView {
    const std::vector<std::vector<int>>* cover_up = nullptr;
    std::vector<int> topo;   // topo order of cover_up
    std::vector<Bitset> R;   // R[u][v] = 1 iff v reachable from u via cover_up (i.e. u <= v)

    bool valid() const {
        return cover_up && (int)R.size() == (int)cover_up->size() && !topo.empty();
    }
};

// Build cached topo + reachability (throws if cover_up is cyclic / invalid)
CoverUpView build_view(const std::vector<std::vector<int>>& cover_up);

// Order test using cached reachability
bool leq(const CoverUpView& V, int a, int b);

// Unique min/max helpers (optional; you can also keep them in poset2.h)
int unique_min_cover(const std::vector<std::vector<int>>& cover_down);
int unique_max_cover(const std::vector<std::vector<int>>& cover_up);

// ImGui UI: compare μ(x,y) between P1 and P2.
// (Implementation should internally cache CoverUpView for both posets.)
void draw_mobius_compare_ui(const pst::Poset1& P1, const pst2::Poset2& P2);

} // namespace mob


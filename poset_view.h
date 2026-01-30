#pragma once
#include "poset_utils.h"

#include <vector>

#include "poset.h"
#include "poset2.h"

struct PosetView {
    int n = 0; // number of elements (nodes)

    // a <= b ?
    bool (*leq)(const void* ctx, int a, int b) = nullptr;
    const void* ctx = nullptr;

    const std::vector<std::vector<int>>* cover_up   = nullptr;
    const std::vector<std::vector<int>>* cover_down = nullptr;

    const std::vector<int>* topo_up   = nullptr;
    const std::vector<int>* topo_down = nullptr;

    // reachability bitsets (optional but you have them)
    const std::vector<pst::Bitset>* reachability_down = nullptr; // down[b] contains all a with a <= b
    const std::vector<pst::Bitset>* reachability_up   = nullptr; // up[a] contains all b with a <= b

    // levels for visualization/layout
    const std::vector<int>* levels = nullptr;
};

inline bool leq(const PosetView& P, int a, int b) {
    return P.leq(P.ctx, a, b);
}

// ------------------------------------------------------------
// Adapters: concrete posets -> PosetView
// ------------------------------------------------------------

inline PosetView view_of(const pst::Poset1& P) {
    PosetView v;
    v.n = static_cast<int>(P.levels.size());

    v.ctx = &P;
    v.leq = [](const void* ctx, int a, int b) -> bool {
        auto& P = *static_cast<const pst::Poset1*>(ctx);
        // a <= b  <=>  a is in down-closure of b
        return P.reachability_down[b].test(a);
    };

    v.cover_up   = &P.cover_up;
    v.cover_down = &P.cover_down;

    v.topo_up   = &P.topo_up;
    v.topo_down = &P.topo_down;

    v.reachability_down = &P.reachability_down;
    v.reachability_up   = &P.reachability_up;

    v.levels = &P.levels;
    return v;
}

inline PosetView view_of(const pst2::Poset2& P) {
    PosetView v;
    v.n = static_cast<int>(P.levels.size());

    v.ctx = &P;
    v.leq = [](const void* ctx, int a, int b) -> bool {
        auto& P = *static_cast<const pst2::Poset2*>(ctx);
        return pst2::leq(P, a, b);
    };

    v.cover_up   = &P.cover_up;
    v.cover_down = &P.cover_down;

    v.topo_up   = &P.topo_up;
    v.topo_down = &P.topo_down;

    v.reachability_down = &P.reachability_down;
    v.reachability_up   = &P.reachability_up;

    v.levels = &P.levels;
    return v;
}


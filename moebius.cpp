#include "moebius.h"
#include "poset_view.h"

#include <iostream>

namespace {

// global minimum = node with no lower cover neighbors
static int find_unique_global_min(const PosetView& P) {
    if (!P.cover_down) return -1;

    int mn = -1;
    for (int v = 0; v < P.n; ++v) {
        if ((*P.cover_down)[v].empty()) {
            if (mn != -1) return -1; // not unique
            mn = v;
        }
    }
    return mn;
}

} // namespace

namespace pst {

std::vector<std::int64_t>
mobius_from_global_min(const PosetView& P) {

    std::vector<std::int64_t> mu(P.n, 0);
    if (P.n == 0) return mu;

    if (!P.cover_down || !P.topo_up || !P.reachability_down) {
        std::cerr << "[moebius] PosetView missing required data\n";
        return mu;
    }

    const int mn = find_unique_global_min(P);
    if (mn < 0) {
        std::cerr << "[moebius] no UNIQUE global minimum found\n";
        return mu;
    }

    mu[mn] = 1;

    // Bottom -> top order (already computed by you)
    for (int x : *P.topo_up) {
        if (x == mn) continue;

        const pst::Bitset& downX = (*P.reachability_down)[x];
        std::int64_t sum = 0;

        // iterate strict predecessors a < x
        for (auto a = downX.find_first();
             a != pst::Bitset::npos;
             a = downX.find_next(a))
        {
            if ((int)a == x) continue;
            sum += mu[(int)a];
        }

        mu[x] = -sum;
    }

    return mu;
}

} // namespace pst





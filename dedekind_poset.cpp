#include "dedekind_poset.h"
#include "dedekind_cut.h"
#include "poset_vis_helpers.h"

#include "poset_utils.h"
#include "meet_join.h"   

#include <algorithm>
#include <iostream>


namespace dm_completion {

// strict subset A ⊂ B: every element of A in B, and A != B
static bool strict_subset(const pst::Bitset& A, const pst::Bitset& B) {
    if (A.count() >= B.count()) return false; // early out
    pst::Bitset tmp = A;
    tmp &= ~B;          // tmp = A \ B
    return tmp.none();  // true iff A\B empty
}

// turn Iprime of a cut element into a bitset over base poset nodes
static pst::Bitset cut_Iprime_bitset(int baseN, const dm_completion::DedekindCut& C) {
    pst::Bitset b(baseN);
    for (int x : C.Iprime) {
        if (0 <= x && x < baseN) b.set((size_t)x);
    }
    return b;
}

// ------------------------------------------------------------
// Build Dedekind–MacNeille completion 
// ------------------------------------------------------------
DedekindPoset build_dedekind_poset(const PosetView& P, std::vector<dm_completion::DedekindCut> cuts) {
    DedekindPoset D;
    D.cuts = std::move(cuts);

    const int M     = (int)D.cuts.size(); // DM size
    const int baseN = P.n;                // base poset size (for I' bitsets)

    D.cover_up.assign(M, {});
    D.cover_down.assign(M, {});
    D.Ip_bit.resize(M);

    D.levels.assign(M, 0);
    D.maxLevel = 0;

    if (M == 0) return D;

    // 1) compute bitsets for I'
    for (int i = 0; i < M; ++i) {
        D.Ip_bit[i] = cut_Iprime_bitset(baseN, D.cuts[i]);
    }

    // 2) full strict order graph: out[i] = all j with I'_i ⊂ I'_j
    std::vector<std::vector<int>> out(M);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            if (i == j) continue;
            if (strict_subset(D.Ip_bit[i], D.Ip_bit[j])) out[i].push_back(j);
        }
    }
    pst::sort_unique_adjacency(out);

    // 3) topo + reachability + transitive reduction => cover_up
    D.topo_up = pst::topo_sort_kahn(out);
    D.reachability_up    = pst::compute_reachability(out, D.topo_up);
    D.cover_up = pst::transitive_reduction(out, D.reachability_up);

    // build cover_dn from cover_up
    for (int u = 0; u < M; ++u) {
        for (int v : D.cover_up[u]) {
            if (0 <= v && v < M) D.cover_down[v].push_back(u);
        }
    }

    D.topo_down = pst::topo_sort_kahn(D.cover_down);
    D.reachability_down = pst::compute_reachability(D.cover_down, D.topo_down);


    // 4) levels for visualization (use your existing helper)
    D.levels = pst_vis::compute_levels_longest_from_roots_cover_up(D.cover_up);
    for (int lv : D.levels) D.maxLevel = std::max(D.maxLevel, lv);

    return D;
}

// ------------------------------------------------------------
// lattice check for EVERY poset 
// Finite lattice <=> complete lattice, so we check unique meet+join for every pair
// ------------------------------------------------------------
bool check_complete_lattice(const PosetView& P, bool verbose) {
    const int n = P.n;
    if (n == 0) return true;

    for (int a = 0; a < n; ++a) {
        for (int b = a; b < n; ++b) {

            auto j = meet_join::join(P, a, b);
            if (!j.exists() || !j.unique()) {
                if (verbose) {
                    std::cerr << "[lattice-check] join not unique/existent for ("
                              << a << "," << b << "), candidates=" << j.candidates.size() << "\n";
                }
                return false;
            }

            auto m = meet_join::meet(P, a, b);
            if (!m.exists() || !m.unique()) {
                if (verbose) {
                    std::cerr << "[lattice-check] meet not unique/existent for ("
                              << a << "," << b << "), candidates=" << m.candidates.size() << "\n";
                }
                return false;
            }
        }
    }

    if (verbose) {
        std::cout << "[lattice-check] OK: poset is a lattice (finite => complete lattice).\n";
    }
    return true;
}

} // namespace pst2

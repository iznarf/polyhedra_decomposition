#include "dedekind_poset.h"

#include <algorithm>
#include <queue>

namespace pst2 {

// strict subset A ⊂ B
static bool strict_subset(const pst::Bitset& A, const pst::Bitset& B) {
    if (A.count() >= B.count()) return false;
    pst::Bitset tmp = A;
    tmp &= ~B;          // A \ B
    return tmp.none();  // empty => subset
}

// turn cut's Iprime into a bitset over base poset nodes
static pst::Bitset cut_Iprime_bitset(int baseN, const DedekindCut& C) {
    pst::Bitset b(baseN);
    for (int x : C.Iprime) {
        if (0 <= x && x < baseN) b.set((size_t)x);
    }
    return b;
}

// cover_out edges u -> v meaning u < v (cover)
// level[v] = longest chain length from any minimal element to v
static std::vector<int> levels_from_cover(const std::vector<std::vector<int>>& cover_out) {
    const int M = (int)cover_out.size();

    std::vector<int> indeg(M, 0);
    for (int u = 0; u < M; ++u)
        for (int v : cover_out[u])
            if (0 <= v && v < M) indeg[v]++;

    std::queue<int> q;
    std::vector<int> level(M, 0);

    for (int i = 0; i < M; ++i)
        if (indeg[i] == 0) q.push(i); // minimal elements level 0

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : cover_out[u]) {
            if (v < 0 || v >= M) continue;
            level[v] = std::max(level[v], level[u] + 1);
            if (--indeg[v] == 0) q.push(v);
        }
    }
    return level;
}

DedekindPoset build_dedekind_poset(const Poset2& P2, std::vector<DedekindCut> cuts)
{
    DedekindPoset D;
    D.cuts = std::move(cuts);

    const int M     = (int)D.cuts.size();
    const int baseN = (int)P2.cover_up.size();

    D.cover_up.assign(M, {});
    D.cover_dn.assign(M, {});
    D.Ip_bit.resize(M);

    if (M == 0) return D;

    // 1) bitsets for I'
    for (int i = 0; i < M; ++i) D.Ip_bit[i] = cut_Iprime_bitset(baseN, D.cuts[i]);

    // 2) full strict order: i -> j iff I'_i ⊂ I'_j
    std::vector<std::vector<int>> out(M);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            if (i == j) continue;
            if (strict_subset(D.Ip_bit[i], D.Ip_bit[j])) out[i].push_back(j);
        }
    }
    pst::sort_unique_adjacency(out);

    // 3) cover edges via transitive reduction
    D.topo = pst::topo_sort_kahn(out);
    D.R    = pst::compute_reachability(out, D.topo);
    D.cover_up = pst::transitive_reduction(out, D.R);

    // build cover_dn
    for (int u = 0; u < M; ++u) {
        for (int v : D.cover_up[u]) {
            if (0 <= v && v < M) D.cover_dn[v].push_back(u);
        }
    }

    // 4) levels
    D.level = levels_from_cover(D.cover_up);
    for (int lv : D.level) D.maxLevel = std::max(D.maxLevel, lv);

    return D;
}



// fast order test using reachability bitsets
static inline bool dm_leq(const DedekindPoset& D, int a, int b) {
    if (a == b) return true;
    return D.R[a].test((size_t)b);
}

// returns all common upper bounds of (a,b)
static std::vector<int> common_upper_bounds(const DedekindPoset& D, int a, int b) {
    const int M = (int)D.cuts.size();
    std::vector<int> U;
    U.reserve(M);
    for (int u = 0; u < M; ++u) {
        if (dm_leq(D, a, u) && dm_leq(D, b, u)) U.push_back(u);
    }
    return U;
}

// returns all common lower bounds of (a,b)
static std::vector<int> common_lower_bounds(const DedekindPoset& D, int a, int b) {
    const int M = (int)D.cuts.size();
    std::vector<int> L;
    L.reserve(M);
    for (int l = 0; l < M; ++l) {
        if (dm_leq(D, l, a) && dm_leq(D, l, b)) L.push_back(l);
    }
    return L;
}

// from a set S, keep only the minimal elements (w.r.t. <=)
static std::vector<int> minimal_elements(const DedekindPoset& D, const std::vector<int>& S) {
    std::vector<int> mins;
    for (int x : S) {
        bool is_min = true;
        for (int y : S) {
            if (y == x) continue;
            if (dm_leq(D, y, x)) { // y < x inside S => x not minimal
                is_min = false;
                break;
            }
        }
        if (is_min) mins.push_back(x);
    }
    return mins;
}

// from a set S, keep only the maximal elements (w.r.t. <=)
static std::vector<int> maximal_elements(const DedekindPoset& D, const std::vector<int>& S) {
    std::vector<int> maxs;
    for (int x : S) {
        bool is_max = true;
        for (int y : S) {
            if (y == x) continue;
            if (dm_leq(D, x, y)) { // x < y inside S => x not maximal
                is_max = false;
                break;
            }
        }
        if (is_max) maxs.push_back(x);
    }
    return maxs;
}

// compute join a∨b: unique minimal common upper bound
static bool compute_join(const DedekindPoset& D, int a, int b, int& outJoin) {
    auto U = common_upper_bounds(D, a, b);
    if (U.empty()) return false;

    auto mins = minimal_elements(D, U);
    if (mins.size() != 1) return false;

    outJoin = mins[0];
    return true;
}

// compute meet a∧b: unique maximal common lower bound
static bool compute_meet(const DedekindPoset& D, int a, int b, int& outMeet) {
    auto L = common_lower_bounds(D, a, b);
    if (L.empty()) return false;

    auto maxs = maximal_elements(D, L);
    if (maxs.size() != 1) return false;

    outMeet = maxs[0];
    return true;
}

// Returns true iff the given finite poset is a lattice (all binary meets/joins exist uniquely).
// For a finite poset, "lattice" implies "complete lattice" (every subset has meet/join).
bool check_complete_lattice(const DedekindPoset& D, bool verbose) {
    const int M = (int)D.cuts.size();
    if (M == 0) return true;

    for (int a = 0; a < M; ++a) {
        for (int b = a; b < M; ++b) {
            int j = -1, m = -1;

            if (!compute_join(D, a, b, j)) {
                if (verbose) std::cerr << "[lattice-check] No unique join for (" << a << "," << b << ")\n";
                return false;
            }
            if (!compute_meet(D, a, b, m)) {
                if (verbose) std::cerr << "[lattice-check] No unique meet for (" << a << "," << b << ")\n";
                return false;
            }
        }
    }

    if (verbose) std::cout << "[lattice-check] OK: poset is a lattice (finite => complete lattice).\n";
    return true;
}





} // namespace pst2

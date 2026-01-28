#include "dedekind_poset.h"

#include <algorithm>
#include <queue>

// for 8 and 9 vertices we can precompute the full DM poset and check completeness


namespace pst2 {

// strict subset A ⊂ B: is every element of A in B, and A != B
static bool strict_subset(const pst::Bitset& A, const pst::Bitset& B) {
    if (A.count() >= B.count()) return false; // if |A| >= |B| then A !⊂ B -> early out 
    pst::Bitset tmp = A; 
    tmp &= ~B;          // tmp = A \ B -> bitwise NOT of B. bits are 1 where B does not contain elements of A
    return tmp.none();  // is true iff all bits are 0, so when A \ B is empty => subset
}

// turn Iprime of a cut element into a bitset over base poset nodes
static pst::Bitset cut_Iprime_bitset(int baseN, const DedekindCut& C) {
    pst::Bitset b(baseN);
    for (int x : C.Iprime) {
        if (0 <= x && x < baseN) b.set((size_t)x);
    }
    return b;
}

// cover_out edges u -> v meaning u < v (cover)
// level[v] = longest chain length from any minimal element to v
// minimal elements have no incoming cover edges 
static std::vector<int> levels_from_cover(const std::vector<std::vector<int>>& cover_out) {
    const int M = (int)cover_out.size();

    // first we compute indegrees: number of incoming cover edges
    std::vector<int> indeg(M, 0);
    for (int u = 0; u < M; ++u)
        for (int v : cover_out[u])
            if (0 <= v && v < M) indeg[v]++;

    
    std::queue<int> q;
    // levels
    std::vector<int> level(M, 0);

    // kahn style BFS from minimal elements 
    // first all minimal elements to queue, set level to 0
    for (int i = 0; i < M; ++i)
        if (indeg[i] == 0) q.push(i); // minimal elements level 0

    
    while (!q.empty()) {
        int u = q.front(); q.pop();
        // iterate over outgoing cover edges from u 
        for (int v : cover_out[u]) {
            if (v < 0 || v >= M) continue;
            // set level of v to longest chain to v
            level[v] = std::max(level[v], level[u] + 1);
            // decrease indegree, if 0 add to queue -> all predecessors processed
            if (--indeg[v] == 0) q.push(v);
        }
    }
    return level;
}

DedekindPoset build_dedekind_poset(const Poset2& P2, std::vector<DedekindCut> cuts){
    DedekindPoset D;
    D.cuts = std::move(cuts);

    // number of DM elements
    const int M     = (int)D.cuts.size(); 
    // number of base poset elements to size bitsets correctly
    const int baseN = (int)P2.cover_up.size();

    D.cover_up.assign(M, {});
    D.cover_dn.assign(M, {});
    // cached I' bitsets: for ordering the Dm elements we need to check for inclusion of I' sets
    // stores bitset for each cut's I': at [i] is bitset for cut D.cuts[i].Iprime
    D.Ip_bit.resize(M);

    if (M == 0) return D;

    // 1) compute bitsets for I'
    for (int i = 0; i < M; ++i) D.Ip_bit[i] = cut_Iprime_bitset(baseN, D.cuts[i]);

    // 2) full strict order: i -> j iff I'_i ⊂ I'_j
    // build full order graph: out [i] = all j with i < j
    std::vector<std::vector<int>> out(M);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            if (i == j) continue;
            if (strict_subset(D.Ip_bit[i], D.Ip_bit[j])) out[i].push_back(j);
        }
    }
    pst::sort_unique_adjacency(out);

    // 3) cover edges via transitive reduction
    // we need topo sort of cuts for reachability computation
    // i < j in topo order => I'_i ⊂ I'_j
    D.topo = pst::topo_sort_kahn(out);
    D.R    = pst::compute_reachability(out, D.topo);
    D.cover_up = pst::transitive_reduction(out, D.R);

    // build cover_down edges from cover_up
    for (int u = 0; u < M; ++u) {
        for (int v : D.cover_up[u]) {
            if (0 <= v && v < M) D.cover_dn[v].push_back(u);
        }
    }

    // 4) compute levels for visualization
    D.level = levels_from_cover(D.cover_up);
    for (int lv : D.level) D.maxLevel = std::max(D.maxLevel, lv);

    return D;
}

// -------------------------------------------------------------------------------------------------
// check if the given Dedekind poset is a complete lattice
// since our poset is finite, it suffices to check if it is a lattice (every two elements have unique meets/joins)
// because finite lattice => complete lattice
// -------------------------------------------------------------------------------------------------

// fast order test using reachability bitsets
static inline bool dm_leq(const DedekindPoset& D, int a, int b) {
    if (a == b) return true;
    return D.R[a].test((size_t)b); // D.R[a][b] = 1 iff a <= b
}

// returns all common upper bounds of (a,b)
// ints are indices of Dedekind cuts in D.cuts
static std::vector<int> common_upper_bounds(const DedekindPoset& D, int a, int b) {
    const int M = (int)D.cuts.size();
    std::vector<int> U;
    // U holds indices of upper bounds of a and b
    U.reserve(M);
    // go through all cuts and check if they are >= a and >= b
    for (int u = 0; u < M; ++u) {
        // if yes, add to U
        if (dm_leq(D, a, u) && dm_leq(D, b, u)) U.push_back(u);
    }
    return U;
}

// returns all common lower bounds of (a,b)
static std::vector<int> common_lower_bounds(const DedekindPoset& D, int a, int b) {
    const int M = (int)D.cuts.size();
    std::vector<int> L;
    // L holds indices of lower bounds of a and b
    L.reserve(M);
    // go through all cuts and check if they are <= a and <= b
    for (int l = 0; l < M; ++l) {
        // if yes, add to L
        if (dm_leq(D, l, a) && dm_leq(D, l, b)) L.push_back(l);
    }
    return L;
}

// from a set S, keep only the minimal elements (w.r.t. <=)
// S is a set of indices of Dedekind cuts in D.cuts
static std::vector<int> minimal_elements(const DedekindPoset& D, const std::vector<int>& S) {
    std::vector<int> mins;
    // check each cut in S
    for (int x : S) {
        bool is_min = true;
        // compare to every other cut in S
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
// S is a set of indices of Dedekind cuts in D.cuts
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

// compute join (a,b) : unique minimal common upper bound
// outJoin is set to the index of the join cut if it exists
static bool compute_join(const DedekindPoset& D, int a, int b, int& outJoin) {
    // U holds indices of common upper bounds of a and b
    auto U = common_upper_bounds(D, a, b);
    if (U.empty()) return false;
    // from U, get minimal elements
    auto mins = minimal_elements(D, U);
    // if not unique, no join
    if (mins.size() != 1) return false;

    outJoin = mins[0];
    return true;
}

// compute meet (a,b) : unique maximal common lower bound
// outMeet is set to the index of the meet cut if it exists
static bool compute_meet(const DedekindPoset& D, int a, int b, int& outMeet) {
    // L holds indices of common lower bounds of a and b
    auto L = common_lower_bounds(D, a, b);
    if (L.empty()) return false;
    // from L, get maximal elements
    // if not unique, no meet
    auto maxs = maximal_elements(D, L);
    if (maxs.size() != 1) return false;

    outMeet = maxs[0];
    return true;
}

// returns true iff the given finite poset is a lattice 
// for a finite poset, lattice implies complete lattice 
// therefore we check the existence of unique meet/join for every pair of elements
bool check_complete_lattice(const DedekindPoset& D, bool verbose) {
    const int M = (int)D.cuts.size();
    if (M == 0) return true;


    for (int a = 0; a < M; ++a) {
        for (int b = a; b < M; ++b) {
            // output indices of join and meet 
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

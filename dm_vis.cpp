#include "dm_vis.h"

#include "poset_utils.h" // pst::topo_sort_kahn, pst::compute_reachability, pst::sort_unique_adjacency

#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <iostream>
#include <queue>
#include <vector>

namespace viz_dm {

static bool g_enabled  = true;
static int  g_selected = 0;

static const std::vector<pst2::DedekindCut>* g_cuts_ptr = nullptr;

static polyscope::CurveNetwork* g_cn_edges = nullptr;
static polyscope::PointCloud*   g_pc_nodes = nullptr;

static const char* kDmEdgesName = "DM completion (edges)";
static const char* kDmNodesName = "DM completion (nodes)";

// strict subset A ⊂ B
static bool strict_subset(const pst::Bitset& A, const pst::Bitset& B) {
    if (A.count() >= B.count()) return false;
    pst::Bitset tmp = A;
    tmp &= ~B;          // A \ B
    return tmp.none();  // empty => subset
}

// turn cut's Iprime into a bitset over base poset nodes
static pst::Bitset cut_Iprime_bitset(int baseN, const pst2::DedekindCut& C) {
    pst::Bitset b(baseN);
    for (int x : C.Iprime) {
        if (0 <= x && x < baseN) b.set((size_t)x);
    }
    return b;
}

// remove old polyscope structures if they exist
static void remove_dm_structures() {
    if (polyscope::hasCurveNetwork(kDmEdgesName)) polyscope::removeCurveNetwork(kDmEdgesName, false);
    if (polyscope::hasPointCloud(kDmNodesName))   polyscope::removePointCloud(kDmNodesName, false);
    g_cn_edges = nullptr;
    g_pc_nodes = nullptr;
}

// ------------------------------------------------------------
// DM lattice levels:
// cover_out is Hasse edges u -> v meaning u < v (cover)
// level[v] = longest chain length from any minimal element to v
static std::vector<int> dm_levels_from_cover(const std::vector<std::vector<int>>& cover_out) {
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

















// ------------------------------------------------------------

void build_and_register_dm_graph(const pst2::Poset2& P2,
                                 const std::vector<pst2::DedekindCut>& cuts)
{
    g_cuts_ptr = &cuts;
    g_selected = 0;

    const int M     = (int)cuts.size();        // number of DM elements
    const int baseN = (int)P2.cover_up.size(); // base poset size

    if (M == 0) {
        std::cerr << "[dm_vis] No cuts -> nothing to draw.\n";
        remove_dm_structures();
        return;
    }

    // 1) bitsets for I' (fast inclusion tests)
    std::vector<pst::Bitset> Ip(M);
    for (int i = 0; i < M; ++i) Ip[i] = cut_Iprime_bitset(baseN, cuts[i]);

    // 2) strict order graph: i -> j iff I'_i ⊂ I'_j   (NOT cover, full comparability)
    std::vector<std::vector<int>> out(M);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            if (i == j) continue;
            if (strict_subset(Ip[i], Ip[j])) out[i].push_back(j);
        }
    }
    pst::sort_unique_adjacency(out);

    // 3) cover edges (Hasse diagram) via transitive reduction
    std::vector<int> topo = pst::topo_sort_kahn(out);
    std::vector<pst::Bitset> R = pst::compute_reachability(out, topo);
    std::vector<std::vector<int>> cover_out = pst::transitive_reduction(out, R); 

    // DEBUG: find sizes and identify candidate top/bottom
    int bestFull = -1, bestEmpty = -1;
    size_t bestFullSz = 0, bestEmptySz = (size_t)baseN + 1;

    for (int i = 0; i < M; ++i) {
        size_t s = Ip[i].count();
        if (s > bestFullSz) { bestFullSz = s; bestFull = i; }
        if (s < bestEmptySz) { bestEmptySz = s; bestEmpty = i; }
    }

    std::cout << "[dm_vis] baseN=" << baseN
            << "  max|I'| node=" << bestFull << " size=" << bestFullSz
            << "  min|I'| node=" << bestEmpty << " size=" << bestEmptySz
            << "\n";

    // DEBUG: show indegree/outdegree of these nodes in the cover graph
    auto indeg_cover = std::vector<int>(M,0);
    for (int u=0; u<M; ++u) for (int v: cover_out[u]) indeg_cover[v]++;

    std::cout << "[dm_vis] cover: topCandidate=" << bestFull
            << " indeg=" << indeg_cover[bestFull]
            << " outdeg=" << cover_out[bestFull].size()
            << "\n";
    std::cout << "[dm_vis] cover: botCandidate=" << bestEmpty
            << " indeg=" << indeg_cover[bestEmpty]
            << " outdeg=" << cover_out[bestEmpty].size()
            << "\n";



    // 4) CORRECT LEVELS from cover graph
    std::vector<int> level = dm_levels_from_cover(cover_out);
    int maxLevel = 0;
    for (int lv : level) maxLevel = std::max(maxLevel, lv);


    // 5) layout (GOOD): center each layer + reduce crossings by barycenter sweeps

    // group nodes by level
    std::vector<std::vector<int>> layers(maxLevel + 1);
    for (int i = 0; i < M; ++i) layers[level[i]].push_back(i);

    // build reverse adjacency for "parents" (incoming cover edges)
    std::vector<std::vector<int>> cover_in(M);
    for (int u = 0; u < M; ++u) {
        for (int v : cover_out[u]) {
            if (0 <= v && v < M) cover_in[v].push_back(u);
        }
    }

    // layout params
    const float XSP = 1.4f;
    const float ZSP = 1.4f;
    const glm::vec3 OFFSET(0.f, 0.f, 0.f); // set e.g. (40,0,0) to move aside

    // current x-position for each node (updated after each sweep)
    std::vector<float> xpos(M, 0.f);

    // helper: assign centered x coords for a layer, and update xpos[]
    auto place_layer_centered = [&](int L) {
        auto& vec = layers[L];
        int k = (int)vec.size();
        for (int j = 0; j < k; ++j) {
            int idx = vec[j];
            float x = (j - 0.5f * (k - 1)) * XSP;
            xpos[idx] = x;
        }
    };

    // init: deterministic order then centered placement
    for (int L = 0; L <= maxLevel; ++L) {
        auto& vec = layers[L];
        std::sort(vec.begin(), vec.end());
        place_layer_centered(L);
    }

    // barycenter ordering: reorder each layer by average x of neighbors
    auto barycenter = [&](int node, bool use_parents) -> float {
        const auto& nbrs = use_parents ? cover_in[node] : cover_out[node];
        if (nbrs.empty()) return xpos[node];
        double s = 0.0;
        for (int nb : nbrs) s += xpos[nb];
        return (float)(s / (double)nbrs.size());
    };

    const int SWEEPS = 6;
    for (int it = 0; it < SWEEPS; ++it) {

        // downward sweep: order layer L by parents in layer L-1
        for (int L = 1; L <= maxLevel; ++L) {
            auto& vec = layers[L];
            std::stable_sort(vec.begin(), vec.end(),
                [&](int a, int b) { return barycenter(a, true) < barycenter(b, true); });
            place_layer_centered(L);
        }

        // upward sweep: order layer L by children in layer L+1
        for (int L = maxLevel - 1; L >= 0; --L) {
            auto& vec = layers[L];
            std::stable_sort(vec.begin(), vec.end(),
                [&](int a, int b) { return barycenter(a, false) < barycenter(b, false); });
            place_layer_centered(L);
        }
    }

    // final points in XZ plane (y = 0), z = level
    std::vector<glm::vec3> pts(M, glm::vec3(0.f));
    for (int i = 0; i < M; ++i) {
        pts[i] = OFFSET + glm::vec3(xpos[i], 0.f, (float)level[i] * ZSP);
    }

    // 6) edges list for Polyscope (cover edges)
    std::vector<std::array<int,2>> edges;
    edges.reserve((size_t)M * 4);
    for (int u = 0; u < M; ++u) {
        for (int v : cover_out[u]) edges.push_back({u, v});
    }

    // 7) register / replace: nodes as PointCloud, edges as CurveNetwork
    remove_dm_structures();

    // nodes
    g_pc_nodes = polyscope::registerPointCloud(kDmNodesName, pts);
    g_pc_nodes->setEnabled(g_enabled);
    g_pc_nodes->setPointRadius(0.0060f, true);

    // show cut index as scalar
    {
        std::vector<double> node_id(M);
        for (int i = 0; i < M; ++i) node_id[i] = (double)i;
        g_pc_nodes->addScalarQuantity("cut index", node_id);
    }

    // edges (thin)
    g_cn_edges = polyscope::registerCurveNetwork(kDmEdgesName, pts, edges);
    g_cn_edges->setEnabled(g_enabled);
    g_cn_edges->setRadius(0.00020f, true);

    std::cout << "[dm_vis] DM graph built: nodes=" << M
              << " cover_edges=" << edges.size()
              << " maxLevel=" << maxLevel << "\n";
}
















void set_dm_enabled(bool enabled) {
    g_enabled = enabled;
    if (g_pc_nodes) g_pc_nodes->setEnabled(g_enabled);
    if (g_cn_edges) g_cn_edges->setEnabled(g_enabled);
}

bool dm_enabled() { return g_enabled; }

void set_selected_dm_node(int idx) {
    if (!g_cuts_ptr) { g_selected = 0; return; }
    int M = (int)g_cuts_ptr->size();
    if (M <= 0) { g_selected = 0; return; }
    g_selected = std::max(0, std::min(idx, M - 1));
}

int selected_dm_node() { return g_selected; }

const pst2::DedekindCut* selected_cut() {
    if (!g_cuts_ptr) return nullptr;
    if (g_selected < 0 || g_selected >= (int)g_cuts_ptr->size()) return nullptr;
    return &(*g_cuts_ptr)[g_selected];
}

} // namespace viz_dm

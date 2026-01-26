#include "dm_vis.h"

// dm_vis should only do layout + rendering.
// The DM lattice construction (cover graph + levels) lives in DedekindPoset now.

#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

namespace viz_dm {

static bool g_enabled  = true;
static int  g_selected = 0;

// pointer to DM cuts stored inside DedekindPoset (must outlive visualization)
static const std::vector<pst2::DedekindCut>* g_cuts_ptr = nullptr;

static polyscope::CurveNetwork* g_cn_edges = nullptr;
static polyscope::PointCloud*   g_pc_nodes = nullptr;

static const char* kDmEdgesName = "DM completion (edges)";
static const char* kDmNodesName = "DM completion (nodes)";

// remove old polyscope structures if they exist
static void remove_dm_structures() {
    if (polyscope::hasCurveNetwork(kDmEdgesName)) polyscope::removeCurveNetwork(kDmEdgesName, false);
    if (polyscope::hasPointCloud(kDmNodesName))   polyscope::removePointCloud(kDmNodesName, false);
    g_cn_edges = nullptr;
    g_pc_nodes = nullptr;
}

// ------------------------------------------------------------

void build_and_register_dm_graph(const pst2::DedekindPoset& D)
{
    g_cuts_ptr = &D.cuts;
    g_selected = 0;

    const int M = (int)D.cuts.size(); // number of DM elements

    if (M == 0) {
        std::cerr << "[dm_vis] No cuts -> nothing to draw.\n";
        remove_dm_structures();
        return;
    }

    // DM lattice structure is precomputed in DedekindPoset
    const auto& cover_out = D.cover_up;
    const auto& level     = D.level;
    const int   maxLevel  = D.maxLevel;

    // 5) layout: center each layer + reduce crossings by barycenter sweeps

    // group nodes by level: longest chain from minimal elements
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

// ------------------------------------------------------------

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


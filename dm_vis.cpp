#include "dm_vis.h"
#include "dedekind_poset.h"
#include "dedekind_cut.h"

#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>

#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

#include <imgui.h>
#include <sstream>


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

void build_and_register_dm_graph(const pst2::DedekindPoset& D){
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

    // ------------------------------------------------------------
    // Simple layout :
    // - group nodes by level
    // - sort deterministically within each layer
    // - center each layer horizontally
    // - IMPORTANT: flip z so that "top element" is at the bottom:
    //     z = (maxLevel - level[i]) * ZSP
    // ------------------------------------------------------------

    // group nodes by level
    std::vector<std::vector<int>> layers(maxLevel + 1);
    for (int i = 0; i < M; ++i) {
        int lv = level[i];
        if (0 <= lv && lv <= maxLevel) layers[lv].push_back(i);
    }

    // layout params (MATCH vis_poset.cpp)
    const float XSP = 4.f;  // NODE_SPACING
    const float ZSP = 4.f;  // LEVEL_SPACING
    const glm::vec3 OFFSET(0.f, 0.f, 0.f); // set e.g. (40,0,0) to move aside
    // shift DM vertically by whole levels so roots align
    const int zOffsetLevels = -1; 


    // per-node x coordinate
    std::vector<float> xpos(M, 0.f);

    // place a layer centered at x=0
    auto place_layer_centered = [&](int L) {
        auto& vec = layers[L];
        int k = (int)vec.size();
        for (int j = 0; j < k; ++j) {
            int idx = vec[j];
            float x = (j - 0.5f * (k - 1)) * XSP;
            xpos[idx] = x;
        }
    };

    // deterministic order per layer + center
    for (int L = 0; L <= maxLevel; ++L) {
        auto& vec = layers[L];
        std::sort(vec.begin(), vec.end());
        place_layer_centered(L);
    }

    // final points in XZ plane (y=0)
    // NOTE: FLIPPED z 
    std::vector<glm::vec3> pts(M, glm::vec3(0.f));
    for (int i = 0; i < M; ++i) {
        float z = (float)(maxLevel - level[i] + zOffsetLevels) * ZSP;

        pts[i] = OFFSET + glm::vec3(xpos[i], 0.f, z);
    }

    // edges list for Polyscope (cover edges)
    std::vector<std::array<int,2>> edges;
    edges.reserve((size_t)M * 4);
    for (int u = 0; u < M; ++u) {
        for (int v : cover_out[u]) edges.push_back({u, v});
    }

    // register / replace: nodes as PointCloud, edges as CurveNetwork
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

    // edges (MATCH vis_poset.cpp thickness)
    g_cn_edges = polyscope::registerCurveNetwork(kDmEdgesName, pts, edges);
    g_cn_edges->setEnabled(g_enabled);
    g_cn_edges->setRadius(0.00037f, true);

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

// ------------------------------------------------------------

// --- helpers -------------------------------------------------

static std::string ints_to_string(const std::vector<int>& v) {
    std::ostringstream oss;
    oss << "{";
    for (size_t i = 0; i < v.size(); ++i) {
        oss << v[i];
        if (i + 1 < v.size()) oss << ", ";
    }
    oss << "}";
    return oss.str();
}
static void sync_selected_from_polyscope_selection() {
    if (!g_pc_nodes || !g_cuts_ptr) return;

    if (!polyscope::haveSelection()) return;

    polyscope::PickResult sel = polyscope::getSelection();
    if (!sel.isHit) return;

    // only react to clicks on our DM point cloud
    if (sel.structure != g_pc_nodes) return;

    int idx = (int)sel.localIndex; // for point clouds: point index
    int M   = (int)g_cuts_ptr->size();
    if (0 <= idx && idx < M) g_selected = idx;
}


// --- UI ------------------------------------------------------

void dm_ui() {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("dedekind poset nodes")) return;

    if (!g_cuts_ptr) {
        ImGui::TextUnformatted("DM: no cuts registered.");
        return;
    }
    const int M = (int)g_cuts_ptr->size();
    if (M <= 0) {
        ImGui::TextUnformatted("DM: empty.");
        return;
    }

    // Update selection from Polyscope clicking
    sync_selected_from_polyscope_selection();


    ImGui::Checkbox("enabled", &g_enabled);
    set_dm_enabled(g_enabled);

    // Manual selection too (nice for debugging)
    ImGui::SliderInt("selected cut", &g_selected, 0, M - 1);

    const auto& C = (*g_cuts_ptr)[g_selected];

    ImGui::Separator();
    ImGui::Text("cut index: %d / %d", g_selected, M - 1);

    // Show minF / maxI' (your struct fields)
    ImGui::Text("minF size: %d", (int)C.minimal_elements_F.size());
    ImGui::TextWrapped("minF = %s", ints_to_string(C.minimal_elements_F).c_str());

    ImGui::Text("maxI' size: %d", (int)C.maximal_elements_Iprime.size());
    ImGui::TextWrapped("maxI' = %s", ints_to_string(C.maximal_elements_Iprime).c_str());

    ImGui::Text("I' size: %d", (int)C.Iprime.size());
    ImGui::TextWrapped("I' = %s", ints_to_string(C.Iprime).c_str());


}


} // namespace viz_dm



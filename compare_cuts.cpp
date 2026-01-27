#include "compare_cuts.h"

#include "vis_poset.h"
#include "render_helpers.h"

#include "poset.h"
#include "input.h"
#include "visualization.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <imgui.h>
#include <glm/glm.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// ------------------------------------------------------------
// Access globals living in vis_poset.cpp (same pattern as compare_nodes.cpp)
// ------------------------------------------------------------
extern const df::InputData* g_input_for_poset;
extern pst::Poset1 g_poset1;

namespace {

// We store each overlay mesh (2D + lifted) and a small UI toggle state.
struct OverlayMeshPair {
    int node = -1;
    polyscope::SurfaceMesh* m2 = nullptr;
    polyscope::SurfaceMesh* m3 = nullptr;
    bool show2 = true;
    bool show3 = true;
};

static std::vector<OverlayMeshPair> g_orange_minF;   // min(F) meshes (both cuts)
static std::vector<OverlayMeshPair> g_purple_maxIp;  // max(I') meshes (both cuts)

static bool g_show_orange_2d = true;
static bool g_show_orange_3d = true;
static bool g_show_purple_2d = true;
static bool g_show_purple_3d = true;

// helper
static void clear_overlay_vec(std::vector<OverlayMeshPair>& V) {
    for (auto& p : V) {
        if (p.m2) { polyscope::removeStructure(p.m2, false); p.m2 = nullptr; }
        if (p.m3) { polyscope::removeStructure(p.m3, false); p.m3 = nullptr; }
    }
    V.clear();
}

static void clear_compare_cuts_overlay() {
    clear_overlay_vec(g_orange_minF);
    clear_overlay_vec(g_purple_maxIp);
}

static void apply_enabled_states() {
    for (auto& p : g_orange_minF) {
        if (p.m2) p.m2->setEnabled(g_show_orange_2d && p.show2);
        if (p.m3) p.m3->setEnabled(g_show_orange_3d && p.show3);
    }
    for (auto& p : g_purple_maxIp) {
        if (p.m2) p.m2->setEnabled(g_show_purple_2d && p.show2);
        if (p.m3) p.m3->setEnabled(g_show_purple_3d && p.show3);
    }
}

static bool valid_base_node(int v) {
    const int n = (int)g_poset1.nodes.size();
    return 0 <= v && v < n;
}

// copy of the overlay builder from compare_nodes.cpp (true overlay = no offsets)
// with extra safety checks to avoid polyscope crashes on empty/invalid geometry.
static void build_one_overlay(int node_idx,
                              polyscope::SurfaceMesh*& out2d,
                              polyscope::SurfaceMesh*& out3d,
                              const char* name2d,
                              const char* name3d,
                              float offx,
                              float offy)
{
    out2d = nullptr;
    out3d = nullptr;

    if (!g_input_for_poset) return;
    if (!valid_base_node(node_idx)) return;

    df::Tri2 tri = g_input_for_poset->tri_poset;
    pst::replay_history_poset(tri, g_poset1.nodes[node_idx].history, *g_input_for_poset);

    auto ids      = viz::present_ids(tri);
    if (ids.empty()) return;

    auto to_local = viz::make_local_index(ids);
    auto faces    = viz_helpers::faces_from_triangles(tri, to_local);
    if (faces.empty()) return;

    float S = 20.0f;

    auto V2 = viz_helpers::make_planar_poset_vertices(
        ids, g_input_for_poset->points2d, 0.0f, 20.0f, S);

    auto V3 = viz_helpers::make_lifted_poset_vertices(
        ids, g_input_for_poset->points2d, 0.0f, 20.0f, S, 0.5f * S);

    if (V2.empty() || V3.empty()) return;

    // guard face indices
    const int nV = (int)V2.size();
    for (auto f : faces) {
        if (f[0] < 0 || f[1] < 0 || f[2] < 0) return;
        if (f[0] >= nV || f[1] >= nV || f[2] >= nV) return;
    }

    for (auto& p : V2) { p.x += offx; p.y += offy; }
    for (auto& p : V3) { p.x += offx; p.y += offy; }

    // Extra safety: remove if name already exists (should not happen with unique names,
    // but this avoids any silent crash from collisions)
    if (polyscope::hasSurfaceMesh(name2d)) polyscope::removeSurfaceMesh(name2d, false);
    if (polyscope::hasSurfaceMesh(name3d)) polyscope::removeSurfaceMesh(name3d, false);

    out2d = polyscope::registerSurfaceMesh(name2d, V2, faces);
    out3d = polyscope::registerSurfaceMesh(name3d, V3, faces);

    viz_helpers::add_global_id_quantity(out2d, ids);
    viz_helpers::add_global_id_quantity(out3d, ids);

    out2d->setEdgeWidth(1.0f);
    out3d->setEdgeWidth(1.0f);
}

static void add_group_meshes(std::vector<OverlayMeshPair>& dst,const std::vector<int>& nodes, const char* tag, const glm::vec3& color, int cA, int cB){
    // true overlay: no offsets
    const float offx = 0.0f;
    const float offy = 0.0f;

    for (int i = 0; i < (int)nodes.size(); ++i) {
        const int node = nodes[i];
        if (!valid_base_node(node)) {
            std::cout << "[compare_cuts] skip invalid base node idx " << node << "\n";
            continue;
        }

        OverlayMeshPair pair;
        pair.node = node;

        // IMPORTANT: unique names to prevent polyscope crashes
        char n2d[256], n3d[256];
        std::snprintf(n2d, sizeof(n2d), "cutcmp A%d B%d %s node%d 2D", cA, cB, tag, node);
        std::snprintf(n3d, sizeof(n3d), "cutcmp A%d B%d %s node%d 3D", cA, cB, tag, node);

        build_one_overlay(node, pair.m2, pair.m3, n2d, n3d, offx, offy);

        if (!pair.m2 || !pair.m3) {
            // failed to build -> skip pushing; avoids dangling UI toggles
            continue;
        }

        pair.m2->setSurfaceColor(color);
        pair.m3->setSurfaceColor(color);

        pair.m2->setTransparency(0.35f);
        pair.m3->setTransparency(0.35f);

        dst.push_back(pair);
    }
}

static std::vector<int> concat_unique(std::vector<int> a, const std::vector<int>& b) {
    a.insert(a.end(), b.begin(), b.end());
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    return a;
}

} // namespace


namespace viz_poset {

void compare_cuts_ui(const pst2::DedekindPoset& D) {

    if (!ImGui::CollapsingHeader("compare cuts (overlay)", 0)) return;

    if (!g_input_for_poset) {
        ImGui::TextUnformatted("compare cuts needs g_input_for_poset (null).");
        return;
    }

    const int M = (int)D.cuts.size();
    if (M <= 0) {
        ImGui::TextUnformatted("(no Dedekind cuts computed yet)");
        return;
    }

    static char bufA[32] = "0";
    static char bufB[32] = "0";

    ImGui::InputText("cut A", bufA, IM_ARRAYSIZE(bufA));
    ImGui::InputText("cut B", bufB, IM_ARRAYSIZE(bufB));

    const int cA = std::atoi(bufA);
    const int cB = std::atoi(bufB);

    auto valid_cut = [&](int c) { return 0 <= c && c < M; };
    if (!valid_cut(cA) || !valid_cut(cB)) {
        ImGui::Text("Valid cut indices: 0..%d", std::max(0, M - 1));
    }

    ImGui::Separator();

    // group toggles
    bool changed = false;

    ImGui::TextUnformatted("Groups:");
    changed |= ImGui::Checkbox("min(F) 2D (orange)", &g_show_orange_2d); ImGui::SameLine();
    changed |= ImGui::Checkbox("min(F) lifted (orange)", &g_show_orange_3d);

    changed |= ImGui::Checkbox("max(I') 2D (purple)", &g_show_purple_2d); ImGui::SameLine();
    changed |= ImGui::Checkbox("max(I') lifted (purple)", &g_show_purple_3d);

    if (changed) apply_enabled_states();

    ImGui::Separator();

    static bool hide_poset_when_comparing = true;
    ImGui::Checkbox("hide poset when comparing", &hide_poset_when_comparing);

    if (ImGui::Button("build cut overlay")) {

        clear_compare_cuts_overlay();
        if (!valid_cut(cA) || !valid_cut(cB)) return;

        std::vector<int> minF_nodes;
        std::vector<int> maxIp_nodes;

        {
            const auto& A = D.cuts[cA];
            const auto& B = D.cuts[cB];

            minF_nodes  = concat_unique(A.minimal_elements_F, B.minimal_elements_F);
            maxIp_nodes = concat_unique(A.maximal_elements_Iprime, B.maximal_elements_Iprime);
        }

        // orange: min(F)
        add_group_meshes(
            g_orange_minF,
            minF_nodes,
            "minF",
            glm::vec3(1.00f, 0.55f, 0.10f),
            cA, cB
        );

        // purple: max(I')
        add_group_meshes(
            g_purple_maxIp,
            maxIp_nodes,
            "maxIp",
            glm::vec3(0.65f, 0.20f, 1.00f),
            cA, cB
        );

        apply_enabled_states();

        if (hide_poset_when_comparing) {
            viz_poset::set_poset_enabled(false);
        }
    }

    ImGui::SameLine();
    if (ImGui::Button("clear cut overlay")) {
        clear_compare_cuts_overlay();
        viz_poset::apply_poset_visibility_from_masks();
    }

    ImGui::Separator();

    // ------------------------------------------------------------
    // Per-mesh toggles
    // ------------------------------------------------------------
    ImGui::TextUnformatted("Per-mesh toggles:");

    if (!g_orange_minF.empty()) {
        ImGui::TextUnformatted("min(F) (orange):");
        for (int i = 0; i < (int)g_orange_minF.size(); ++i) {
            auto& p = g_orange_minF[i];

            char lab3[128];
            std::snprintf(lab3, sizeof(lab3), "min(F) node %d lifted##minF3_%d", p.node, i);
            bool ch3 = ImGui::Checkbox(lab3, &p.show3);

            ImGui::SameLine();
            char lab2[128];
            std::snprintf(lab2, sizeof(lab2), "2D##minF2_%d", i);
            bool ch2 = ImGui::Checkbox(lab2, &p.show2);

            if (ch2 || ch3) apply_enabled_states();
        }
    }

    if (!g_purple_maxIp.empty()) {
        ImGui::TextUnformatted("max(I') (purple):");
        for (int i = 0; i < (int)g_purple_maxIp.size(); ++i) {
            auto& p = g_purple_maxIp[i];

            char lab3[128];
            std::snprintf(lab3, sizeof(lab3), "max(I') node %d lifted##maxIp3_%d", p.node, i);
            bool ch3 = ImGui::Checkbox(lab3, &p.show3);

            ImGui::SameLine();
            char lab2[128];
            std::snprintf(lab2, sizeof(lab2), "2D##maxIp2_%d", i);
            bool ch2 = ImGui::Checkbox(lab2, &p.show2);

            if (ch2 || ch3) apply_enabled_states();
        }
    }
}

} // namespace viz_poset


#include "compare_nodes.h"

#include "vis_poset.h"
#include "render_helpers.h"

#include "poset.h"
#include "input.h"
#include "visualization.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <imgui.h>
#include <glm/glm.hpp>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>

// ------------------------------------------------------------
// Access globals living in vis_poset.cpp
// ------------------------------------------------------------
extern const df::InputData* g_input_for_poset;
extern pst::Poset1 g_poset1;

// compare UI state
extern char g_cmp_a_buf[32];
extern char g_cmp_b_buf[32];

extern bool g_cmp_show_A_2d;
extern bool g_cmp_show_A_3d;
extern bool g_cmp_show_B_2d;
extern bool g_cmp_show_B_3d;

extern bool g_compare_mode;

// compare meshes (single A/B)
extern polyscope::SurfaceMesh* g_cmpA_2d;
extern polyscope::SurfaceMesh* g_cmpB_2d;
extern polyscope::SurfaceMesh* g_cmpA_3d;
extern polyscope::SurfaceMesh* g_cmpB_3d;

// compare meshes for sets
static std::vector<polyscope::SurfaceMesh*> g_cmpA_set_2d;
static std::vector<polyscope::SurfaceMesh*> g_cmpA_set_3d;
static std::vector<polyscope::SurfaceMesh*> g_cmpB_set_2d;
static std::vector<polyscope::SurfaceMesh*> g_cmpB_set_3d;

// ------------------------------------------------------------
// Internal helpers
// ------------------------------------------------------------
namespace {

void clear_mesh_vec(std::vector<polyscope::SurfaceMesh*>& v) {
    for (auto* m : v) if (m) m->remove();
    v.clear();
}

void clear_compare_meshes() {
    // single A/B
    if (g_cmpA_2d) { g_cmpA_2d->remove(); g_cmpA_2d = nullptr; }
    if (g_cmpB_2d) { g_cmpB_2d->remove(); g_cmpB_2d = nullptr; }
    if (g_cmpA_3d) { g_cmpA_3d->remove(); g_cmpA_3d = nullptr; }
    if (g_cmpB_3d) { g_cmpB_3d->remove(); g_cmpB_3d = nullptr; }

    // sets
    clear_mesh_vec(g_cmpA_set_2d);
    clear_mesh_vec(g_cmpA_set_3d);
    clear_mesh_vec(g_cmpB_set_2d);
    clear_mesh_vec(g_cmpB_set_3d);
}

void set_compare_enabled() {
    // single A/B
    if (g_cmpA_2d) g_cmpA_2d->setEnabled(g_cmp_show_A_2d);
    if (g_cmpA_3d) g_cmpA_3d->setEnabled(g_cmp_show_A_3d);
    if (g_cmpB_2d) g_cmpB_2d->setEnabled(g_cmp_show_B_2d);
    if (g_cmpB_3d) g_cmpB_3d->setEnabled(g_cmp_show_B_3d);

    // sets
    for (auto* m : g_cmpA_set_2d) if (m) m->setEnabled(g_cmp_show_A_2d);
    for (auto* m : g_cmpA_set_3d) if (m) m->setEnabled(g_cmp_show_A_3d);
    for (auto* m : g_cmpB_set_2d) if (m) m->setEnabled(g_cmp_show_B_2d);
    for (auto* m : g_cmpB_set_3d) if (m) m->setEnabled(g_cmp_show_B_3d);
}

static void build_one_overlay(int node_idx,
                              polyscope::SurfaceMesh*& out2d,
                              polyscope::SurfaceMesh*& out3d,
                              const char* name2d,
                              const char* name3d,
                              float offx,
                              float offy)
{
    df::Tri2 tri = g_input_for_poset->tri_poset;
    pst::replay_history_poset(tri, g_poset1.nodes[node_idx].history, *g_input_for_poset);

    auto ids      = viz::present_ids(tri);
    auto to_local = viz::make_local_index(ids);
    auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

    float S = 20.0f;

    auto V2 = viz_helpers::make_planar_poset_vertices(
        ids, g_input_for_poset->points2d, 0.0f, 20.0f, S);

    auto V3 = viz_helpers::make_lifted_poset_vertices(
        ids, g_input_for_poset->points2d, 0.0f, 20.0f, S, 0.5f * S);

    for (auto& p : V2) { p.x += offx; p.y += offy; }
    for (auto& p : V3) { p.x += offx; p.y += offy; }

    out2d = polyscope::registerSurfaceMesh(name2d, V2, faces);
    out3d = polyscope::registerSurfaceMesh(name3d, V3, faces);

    viz_helpers::add_global_id_quantity(out2d, ids);
    viz_helpers::add_global_id_quantity(out3d, ids);

    out2d->setEdgeWidth(1.0f);
    out3d->setEdgeWidth(1.0f);
}




} // namespace



// ============================================================
// PUBLIC API
// ============================================================
namespace viz_poset {

void clear_compare_overlay() {
    clear_compare_meshes();
    g_compare_mode = false;
    viz_poset::apply_poset_visibility_from_masks();
}

// ------------------------------------------------------------
// Original compare-nodes UI 
// ------------------------------------------------------------
void compare_nodes_ui() {

    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("compare nodes (overlay)",
                                ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }

    static bool hide_poset_when_comparing = true;
    ImGui::Checkbox("hide poset when comparing", &hide_poset_when_comparing);

    ImGui::InputText("x", g_cmp_a_buf, IM_ARRAYSIZE(g_cmp_a_buf));
    ImGui::InputText("y", g_cmp_b_buf, IM_ARRAYSIZE(g_cmp_b_buf));

    bool changed = false;
    changed |= ImGui::Checkbox("x 2D mesh", &g_cmp_show_A_2d);
    ImGui::SameLine();
    changed |= ImGui::Checkbox("x 3D mesh", &g_cmp_show_A_3d);

    changed |= ImGui::Checkbox("y 2D mesh", &g_cmp_show_B_2d);
    ImGui::SameLine();
    changed |= ImGui::Checkbox("y 3D mesh", &g_cmp_show_B_3d);

    if (changed) set_compare_enabled();

    const int a = std::atoi(g_cmp_a_buf);
    const int b = std::atoi(g_cmp_b_buf);
    const int n = (int)g_poset1.nodes.size();

    auto valid_idx = [&](int v) { return 0 <= v && v < n; };
    const bool a_ok = valid_idx(a);
    const bool b_ok = valid_idx(b);

    if (!a_ok || !b_ok) {
        ImGui::Text("Valid node indices: 0..%d", std::max(0, n - 1));
    }

    if (ImGui::Button("build compare overlay")) {

        clear_compare_meshes();

        if (!g_input_for_poset || !a_ok || !b_ok) return;

        auto build_one = [&](int node_idx,
                             polyscope::SurfaceMesh*& out2d,
                             polyscope::SurfaceMesh*& out3d,
                             const std::string& tag) {

            df::Tri2 tri = g_input_for_poset->tri_poset;
            pst::replay_history_poset(
                tri, g_poset1.nodes[node_idx].history, *g_input_for_poset);

            auto ids      = viz::present_ids(tri);
            auto to_local = viz::make_local_index(ids);
            auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

            float S = 20.0f;

            auto V2 = viz_helpers::make_planar_poset_vertices(
                ids, g_input_for_poset->points2d, 0.0f, 20.0f, S);

            auto V3 = viz_helpers::make_lifted_poset_vertices(
                ids, g_input_for_poset->points2d, 0.0f, 20.0f, S, 0.5f * S);

            out2d = polyscope::registerSurfaceMesh("compare " + tag + " 2D", V2, faces);
            out3d = polyscope::registerSurfaceMesh("compare " + tag + " lifted", V3, faces);

            viz_helpers::add_global_id_quantity(out2d, ids);
            viz_helpers::add_global_id_quantity(out3d, ids);

            out2d->setEdgeWidth(1.0f);
            out3d->setEdgeWidth(1.0f);
        };

        build_one(a, g_cmpA_2d, g_cmpA_3d, "A");
        build_one(b, g_cmpB_2d, g_cmpB_3d, "B");

        set_compare_enabled();
        g_compare_mode = true;

        if (hide_poset_when_comparing) {
            viz_poset::set_poset_enabled(false);
        }
    }

    ImGui::SameLine();
    if (ImGui::Button("clear compare overlay")) {
        clear_compare_meshes();
        g_compare_mode = false;
        viz_poset::apply_poset_visibility_from_masks();
    }
}

// ------------------------------------------------------------
// compare overlay from SETS (used by DM completion UI)
// ------------------------------------------------------------
void build_compare_overlay_from_sets(const std::vector<int>& A,
                                    const std::vector<int>& B,
                                    const char* tagA,
                                    const char* tagB,
                                    bool hide_poset_when_comparing)
{
    clear_compare_meshes();

    if (!g_input_for_poset) {
        std::cerr << "[vis_poset] compare overlay sets: g_input_for_poset is null\n";
        return;
    }

    const int n = (int)g_poset1.nodes.size();
    auto valid_idx = [&](int v) { return 0 <= v && v < n; };

    // true overlay: no offsets
    const float offx = 0.0f;
    const float offy = 0.0f;

    // build A set (minF) - BLUE
    for (int i = 0; i < (int)A.size(); ++i) {
        int node = A[i];
        if (!valid_idx(node)) continue;

        polyscope::SurfaceMesh* m2 = nullptr;
        polyscope::SurfaceMesh* m3 = nullptr;

        char n2d[128], n3d[128];
        std::snprintf(n2d, sizeof(n2d), "%s #%d (node %d) 2D", tagA, i, node);
        std::snprintf(n3d, sizeof(n3d), "%s #%d (node %d) lifted", tagA, i, node);

        build_one_overlay(node, m2, m3, n2d, n3d, offx, offy);

        if (m2) m2->setSurfaceColor(glm::vec3(0.20f, 0.45f, 1.00f));
        if (m3) m3->setSurfaceColor(glm::vec3(0.20f, 0.45f, 1.00f));

        g_cmpA_set_2d.push_back(m2);
        g_cmpA_set_3d.push_back(m3);
    }

    // build B set (maxI') - RED
    for (int i = 0; i < (int)B.size(); ++i) {
        int node = B[i];
        if (!valid_idx(node)) continue;

        polyscope::SurfaceMesh* m2 = nullptr;
        polyscope::SurfaceMesh* m3 = nullptr;

        char n2d[128], n3d[128];
        std::snprintf(n2d, sizeof(n2d), "%s #%d (node %d) 2D", tagB, i, node);
        std::snprintf(n3d, sizeof(n3d), "%s #%d (node %d) lifted", tagB, i, node);

        build_one_overlay(node, m2, m3, n2d, n3d, offx, offy);

        if (m2) m2->setSurfaceColor(glm::vec3(1.00f, 0.25f, 0.25f));
        if (m3) m3->setSurfaceColor(glm::vec3(1.00f, 0.25f, 0.25f));

        g_cmpB_set_2d.push_back(m2);
        g_cmpB_set_3d.push_back(m3);
    }

    set_compare_enabled();

    g_compare_mode = true;
    if (hide_poset_when_comparing) {
        viz_poset::set_poset_enabled(false);
    }
}


} // namespace viz_poset




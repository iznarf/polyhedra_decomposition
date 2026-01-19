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
#include <cstdio>   // std::snprintf
#include <utility>  // std::swap

// Access globals living in vis_poset.cpp
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

// compare meshes
extern polyscope::SurfaceMesh* g_cmpA_2d;
extern polyscope::SurfaceMesh* g_cmpB_2d;
extern polyscope::SurfaceMesh* g_cmpA_3d;
extern polyscope::SurfaceMesh* g_cmpB_3d;

namespace {

// compare helpers 
void clear_compare_meshes() {
    if (g_cmpA_2d) { g_cmpA_2d->remove(); g_cmpA_2d = nullptr; }
    if (g_cmpB_2d) { g_cmpB_2d->remove(); g_cmpB_2d = nullptr; }
    if (g_cmpA_3d) { g_cmpA_3d->remove(); g_cmpA_3d = nullptr; }
    if (g_cmpB_3d) { g_cmpB_3d->remove(); g_cmpB_3d = nullptr; }
}

void set_compare_enabled() {
    if (g_cmpA_2d) g_cmpA_2d->setEnabled(g_cmp_show_A_2d);
    if (g_cmpA_3d) g_cmpA_3d->setEnabled(g_cmp_show_A_3d);
    if (g_cmpB_2d) g_cmpB_2d->setEnabled(g_cmp_show_B_2d);
    if (g_cmpB_3d) g_cmpB_3d->setEnabled(g_cmp_show_B_3d);
}

} // namespace



namespace viz_poset {

void compare_nodes_ui() {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("compare nodes (overlay)", ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }

    // Optional UX toggle: keep or hide the big poset while comparing
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

        if (!g_input_for_poset) {
            std::cerr << "[vis_poset] compare overlay: g_input_for_poset is null\n";
            return;
        }
        if (!a_ok || !b_ok) {
            std::cerr << "[vis_poset] compare overlay: invalid node indices A=" << a << " B=" << b << "\n";
            return;
        }

        auto build_one = [&](int node_idx,
                             polyscope::SurfaceMesh*& out2d,
                             polyscope::SurfaceMesh*& out3d,
                             const std::string& tag) {

            // reconstruct triangulation at node_idx
            df::Tri2 tri = g_input_for_poset->tri_poset;
            pst::replay_history_poset(tri, g_poset1.nodes[node_idx].history, *g_input_for_poset);

            auto ids      = viz::present_ids(tri);
            auto to_local = viz::make_local_index(ids);
            auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

            float S = 20.0f; // keep your fixed choice

            // keep your current offset choice
            auto V2 = viz_helpers::make_planar_poset_vertices(ids, g_input_for_poset->points2d,
                                                              0.0f, 20.0f, S);

            auto V3 = viz_helpers::make_lifted_poset_vertices(ids, g_input_for_poset->points2d,
                                                              0.0f, 20.0f, S, 0.5f * S);

            // Fixed names so rebuilding doesn’t accumulate weirdly
            std::string name2d = "compare " + tag + " 2D";
            std::string name3d = "compare " + tag + " lifted";

            out2d = polyscope::registerSurfaceMesh(name2d, V2, faces);
            out3d = polyscope::registerSurfaceMesh(name3d, V3, faces);

            viz_helpers::add_global_id_quantity(out2d, ids);
            viz_helpers::add_global_id_quantity(out3d, ids);

            out2d->setEdgeWidth(1.0f);
            out2d->setEdgeColor(glm::vec3(0,0,0));

            out3d->setEdgeWidth(1.0f);
            out3d->setEdgeColor(glm::vec3(0,0,0));
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

        // restore poset view (respects focus/interval if active)
        viz_poset::apply_poset_visibility_from_masks();
    }
}

} // namespace viz_poset



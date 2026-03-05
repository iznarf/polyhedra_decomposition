#include "poset_vis_ui.h"
#include "input.h"
#include "poset.h"
#include "poset2.h"
#include "poset_vis.h"

#include <imgui.h>
#include <glm/glm.hpp>

#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>

// global input data (from main)
extern df::InputData g_in;

extern pst::Poset1  g_P1; // global poset1 (from main)
extern pst2::Poset2 g_P2; // global poset2 (from main)
extern pst::Poset_just_flips g_P1_flips;

extern bool g_has_P1;       // indicates if g_P1 is already built
extern bool g_has_P2;       // indicates if g_P2 is already built
extern bool g_has_P1_flips; // indicates if g_P1_flips is already built

namespace pst_vis_ui {

void poset_ui() {

    // toggles for P1/P2 triangulation mesh visibility
    // P1 2D, 3D, edges, grid
    static bool showP1_2D = false;
    static bool showP1_3D = false;
    static bool showP1_edges = false;
    static bool showP1_grid = false;

    // P2 2D, 3D, edges, grid
    static bool showP2_2D = false;
    static bool showP2_3D = false;
    static bool showP2_edges = false;
    static bool showP2_grid = false;

    // P1 flips 2D, 3D, edges, grid
    static bool showP1_flips_2D = false;
    static bool showP1_flips_3D = false;
    static bool showP1_flips_edges = false;
    static bool showP1_flips_grid = false;

    ImGui::SeparatorText("Poset visibility");
    ImGui::Checkbox("P1 2D meshes", &showP1_2D);
    ImGui::SameLine();
    ImGui::Checkbox("P1 3D meshes", &showP1_3D);
    ImGui::SameLine();
    ImGui::Checkbox("P1 grid", &showP1_grid);
    ImGui::SameLine();
    ImGui::Checkbox("P1 edges", &showP1_edges);

    ImGui::Checkbox("P2 2D meshes", &showP2_2D);
    ImGui::SameLine();
    ImGui::Checkbox("P2 3D meshes", &showP2_3D);
    ImGui::SameLine();
    ImGui::Checkbox("P2 grid", &showP2_grid);
    ImGui::SameLine();
    ImGui::Checkbox("P2 edges", &showP2_edges);

    ImGui::Checkbox("P1 flips 2D meshes", &showP1_flips_2D);
    ImGui::SameLine();
    ImGui::Checkbox("P1 flips 3D meshes", &showP1_flips_3D);
    ImGui::SameLine();
    ImGui::Checkbox("P1 flips grid", &showP1_flips_grid);
    ImGui::SameLine();
    ImGui::Checkbox("P1 flips edges", &showP1_flips_edges);

    // apply every frame
    pst_vis::set_trivis_enabled(pst_vis::g_triVis_P1, showP1_2D, showP1_3D);
    pst_vis::set_trivis_enabled(pst_vis::g_triVis_P2, showP2_2D, showP2_3D);
    pst_vis::set_trivis_enabled(pst_vis::g_triVis_P1_flips, showP1_flips_2D, showP1_flips_3D);

    // ---- P1 grid/edges ----
    if (polyscope::hasPointCloud("P1 grid nodes")) {
        polyscope::getPointCloud("P1 grid nodes")->setEnabled(showP1_grid);
    }
    if (polyscope::hasCurveNetwork("P1 edges")) {
        polyscope::getCurveNetwork("P1 edges")->setEnabled(showP1_edges);
    }

    // ---- P1 flips grid/edges ----
    if (polyscope::hasPointCloud("P1 flips grid nodes")) {
        polyscope::getPointCloud("P1 flips grid nodes")->setEnabled(showP1_flips_grid);
    }
    if (polyscope::hasCurveNetwork("P1 flips edges")) {
        polyscope::getCurveNetwork("P1 flips edges")->setEnabled(showP1_flips_edges);
    }

    // ---- P2 grid/edges ----
    if (polyscope::hasPointCloud("P2 grid nodes")) {
        polyscope::getPointCloud("P2 grid nodes")->setEnabled(showP2_grid);
    }
    if (polyscope::hasCurveNetwork("P2 edges")) {
        polyscope::getCurveNetwork("P2 edges")->setEnabled(showP2_edges);
    }

    ImGui::Separator();

    if (ImGui::Button("Visualize poset 1")) {
        if (g_has_P1) {
            pst_vis::register_poset1(
                g_in, g_P1,
                glm::vec3(0.f, 0.f, 0.f),
                0.2f, 0.2f,
                glm::vec3(1.0f, 0.0f, 0.0f)
            );
            showP1_2D = true;
            showP1_3D = true;
            showP1_grid = true;
            showP1_edges = true;
        }
    }

    if (ImGui::Button("Visualize poset 2")) {
        if (g_has_P1 && g_has_P2) {
            pst_vis::register_poset2(
                g_in, g_P1, g_P2,
                glm::vec3(0.f, 0.f, 0.f),
                0.2f, 0.2f,
                glm::vec3(13.0f / 255.0f, 100.0f / 255.0f, 13.0f / 255.0f)
            );
            showP2_2D = true;
            showP2_3D = true;
            showP2_grid = true;
            showP2_edges = true;
        }
    }

    if (ImGui::Button("Visualize poset 1 just flips")) {
        if (g_has_P1_flips) {
            pst_vis::register_poset1_flips(
                g_in, g_P1_flips,
                glm::vec3(0.f, 0.f, 0.f),
                0.3f, 0.3f,
                glm::vec3(0.9f, 0.2f, 0.9f)
            );

            showP1_flips_2D = true;
            showP1_flips_3D = true;
            showP1_flips_grid = true;
            showP1_flips_edges = true;
        }
    }
}

} // namespace pst_vis_ui

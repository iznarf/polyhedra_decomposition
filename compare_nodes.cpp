#include "compare_nodes.h"

#include "render_helpers.h"
#include "poset.h"
#include "poset2.h"
#include "input.h"
#include "visualization.h"

#include "envelope_vis.h"   // uses: compute_and_register_from_P1_nodes + build_mesh_pair_from_P1_node

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>

#include <imgui.h>
#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// ------------------------------------------------------------
// Globals
// ------------------------------------------------------------
extern df::InputData g_in;

extern pst::Poset1  g_P1;
extern pst2::Poset2 g_P2;
extern bool g_has_P1;
extern bool g_has_P2;

namespace viz_poset {

namespace {

// ------------------------------------------------------------
// UI state
// ------------------------------------------------------------
static char g_a_buf[32] = "0";
static char g_b_buf[32] = "0";

static bool g_show_A_2d = true;
static bool g_show_A_3d = true;
static bool g_show_B_2d = true;
static bool g_show_B_3d = true;

static bool g_show_E_2d = true;
static bool g_show_E_3d = true;

// overlay meshes
static polyscope::SurfaceMesh* g_A2 = nullptr;
static polyscope::SurfaceMesh* g_A3 = nullptr;
static polyscope::SurfaceMesh* g_B2 = nullptr;
static polyscope::SurfaceMesh* g_B3 = nullptr;

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------
static void remove_mesh(polyscope::SurfaceMesh*& m) {
  if (m) { m->remove(); m = nullptr; }
}

static void clear_compare_overlay_internal() {
  remove_mesh(g_A2);
  remove_mesh(g_A3);
  remove_mesh(g_B2);
  remove_mesh(g_B3);
}

static void apply_enabled_toggles_overlays_only() {
  if (g_A2) g_A2->setEnabled(g_show_A_2d);
  if (g_A3) g_A3->setEnabled(g_show_A_3d);
  if (g_B2) g_B2->setEnabled(g_show_B_2d);
  if (g_B3) g_B3->setEnabled(g_show_B_3d);
}

static void apply_enabled_toggles_all() {
  apply_enabled_toggles_overlays_only();
  envelope_vis::set_enabled(g_show_E_2d, g_show_E_3d);
}

static std::string prefix_for(int whichPoset) {
  return (whichPoset == 1) ? "P1 " : "P2 ";
}

} // namespace

// ============================================================
// UI
// ============================================================
void compare_nodes_ui() {

  ImGui::SetNextItemOpen(false, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("Compare", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  // choose poset
  static int which = 2;
  if (!g_has_P2) which = 1;

  const bool canP1 = g_has_P1;
  const bool canP2 = g_has_P2;

  ImGui::TextUnformatted("Target poset:");
  ImGui::BeginDisabled(!canP1);
  ImGui::RadioButton("P1", &which, 1);
  ImGui::EndDisabled();

  ImGui::SameLine();

  ImGui::BeginDisabled(!canP2);
  ImGui::RadioButton("P2", &which, 2);
  ImGui::EndDisabled();

  if ((which == 1 && !canP1) || (which == 2 && !canP2)) {
    ImGui::TextUnformatted("Selected poset not available.");
    return;
  }

  // Inputs
  ImGui::InputText("x", g_a_buf, IM_ARRAYSIZE(g_a_buf));
  ImGui::InputText("y", g_b_buf, IM_ARRAYSIZE(g_b_buf));

  bool changed = false;
  changed |= ImGui::Checkbox("x 2D", &g_show_A_2d); ImGui::SameLine();
  changed |= ImGui::Checkbox("x 3D", &g_show_A_3d);
  changed |= ImGui::Checkbox("y 2D", &g_show_B_2d); ImGui::SameLine();
  changed |= ImGui::Checkbox("y 3D", &g_show_B_3d);
  changed |= ImGui::Checkbox("env 2D", &g_show_E_2d); ImGui::SameLine();
  changed |= ImGui::Checkbox("env 3D", &g_show_E_3d);

  if (changed) apply_enabled_toggles_all();

  const int a = std::atoi(g_a_buf);
  const int b = std::atoi(g_b_buf);

  // P1 and P2 share node IDs, so we always validate against g_P1 here
  const int n = (int)g_P1.nodes.size();
  auto valid_idx = [&](int v) { return 0 <= v && v < n; };

  const bool a_ok = valid_idx(a);
  const bool b_ok = valid_idx(b);

  if (!a_ok || !b_ok) {
    ImGui::Text("Valid node indices: 0..%d", std::max(0, n - 1));
  }

  // Buttons row 1
  if (ImGui::Button("Build overlay")) {

    clear_compare_overlay_internal();

    if (!a_ok || !b_ok) {
      std::cout << "[compare_nodes] cannot build overlay (invalid ids)\n";
    } else {
      const std::string pref = prefix_for(which);

      // ONE shared function for triangulation meshes
      {
        auto pair = envelope_vis::build_mesh_pair_from_P1_node(
            pref + "compare x 2D",
            pref + "compare x 3D",
            g_in, g_P1, a,
            /*lift_min*/0.0f, /*lift_max*/2.0f, /*scale*/1.0f,
            /*edge_width*/1.0f,
            /*add_global_id_quantity*/true);
        g_A2 = pair.m2;
        g_A3 = pair.m3;
      }

      {
        auto pair = envelope_vis::build_mesh_pair_from_P1_node(
            pref + "compare y 2D",
            pref + "compare y 3D",
            g_in, g_P1, b,
            /*lift_min*/0.0f, /*lift_max*/2.0f, /*scale*/1.0f,
            /*edge_width*/1.0f,
            /*add_global_id_quantity*/true);
        g_B2 = pair.m2;
        g_B3 = pair.m3;
      }

      apply_enabled_toggles_all();
    }
  }

  ImGui::SameLine();

  // NOTE: button name kept as-is even though this computes min envelope of {x,y}
  if (ImGui::Button("Compute min envelope of min(F)")) {

    if (!a_ok || !b_ok) {
      std::cout << "[compare_nodes] cannot compute envelope (invalid ids)\n";
    } else {

      envelope_vis::Options opt;
      opt.show_2d   = g_show_E_2d;
      opt.show_3d   = g_show_E_3d;
      opt.scale     = 1.0f;
      opt.lift_min  = 0.0f;
      opt.lift_max  = 2.0f;
      opt.edge_width = 1.0f;
      opt.reserve_tris = 2000;

      const std::string pref = prefix_for(which);

      const bool ok = envelope_vis::compute_and_register_from_P1_nodes(
          pref + "min envelope 2D",
          pref + "min envelope 3D",
          g_in, g_P1,
          std::vector<int>{a, b},
          opt);

      if (!ok) {
        std::cout << "[compare_nodes] envelope computation failed\n";
      } else {
        apply_enabled_toggles_all();
      }
    }
  }

  // buttons
  if (ImGui::Button("Clear overlay")) {
    clear_compare_overlay_internal();
  }

  ImGui::SameLine();

  if (ImGui::Button("Clear envelope")) {
    envelope_vis::clear();

    if (polyscope::hasCurveNetwork("envelope diagram")) {
      polyscope::removeStructure("envelope diagram");
    }
  }

  ImGui::SameLine();

  if (ImGui::Button("Clear all")) {
    clear_compare_overlay_internal();
    envelope_vis::clear();

    if (polyscope::hasCurveNetwork("envelope diagram")) {
      polyscope::removeStructure("envelope diagram");
    }
  }
}

} // namespace viz_poset







#include "dedekind_envelope_ui.h"

#include "envelope_vis.h"

#include "render_helpers.h"
#include "poset.h"
#include "poset2.h"
#include "input.h"
#include "visualization.h"

#include "dedekind_cut.h"
#include "dedekind_poset.h"
#include "dedekind_completion.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>

#include <imgui.h>
#include <glm/glm.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

extern df::InputData g_in;

extern pst::Poset1  g_P1;
extern pst2::Poset2 g_P2;
extern bool g_has_P1;
extern bool g_has_P2;

extern dm_completion::DedekindPoset g_completion_P1;
extern dm_completion::DedekindPoset g_completion_P2;
extern bool g_hasCompletion_P1;
extern bool g_hasCompletion_P2;

namespace dm_completion {

namespace {

// ------------------------------------------------------------
// UI state
// ------------------------------------------------------------
static char g_cut_buf[32] = "0";

// 1 = P1 completion, 2 = P2 completion
static int g_completionWhich = 1;

static bool g_show_env_2d = true;
static bool g_show_env_3d = true;

static bool g_show_cut_nodes_2d = true;
static bool g_show_cut_nodes_3d = true;

// cut node meshes (just for toggling + clearing)
static std::vector<polyscope::SurfaceMesh*> g_cut_nodes_2d;
static std::vector<polyscope::SurfaceMesh*> g_cut_nodes_3d;

// ------------------------------------------------------------
// helpers
// ------------------------------------------------------------
static void clear_cut_node_meshes() {
  for (auto* m : g_cut_nodes_2d) if (m) m->remove();
  for (auto* m : g_cut_nodes_3d) if (m) m->remove();
  g_cut_nodes_2d.clear();
  g_cut_nodes_3d.clear();
}

// get selected completion pointer (or nullptr)
static const DedekindPoset* selected_completion() {
  if (g_completionWhich == 1) return g_hasCompletion_P1 ? &g_completion_P1 : nullptr;
  return g_hasCompletion_P2 ? &g_completion_P2 : nullptr;
}

} // namespace

// ============================================================
// UI
// ============================================================
void dedekind_envelope_ui() {

  ImGui::SetNextItemOpen(false, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("Envelope", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  // --- choose which completion we browse ---
  ImGui::TextUnformatted("Completion source:");
  ImGui::BeginDisabled(!g_hasCompletion_P1);
  ImGui::RadioButton("P1 completion", &g_completionWhich, 1);
  ImGui::EndDisabled();

  ImGui::SameLine();

  ImGui::BeginDisabled(!g_hasCompletion_P2);
  ImGui::RadioButton("P2 completion", &g_completionWhich, 2);
  ImGui::EndDisabled();

  const DedekindPoset* comp = selected_completion();
  if (!comp) {
    ImGui::TextUnformatted("No completion available for this choice yet.");
    return;
  }

  const int C = (int)comp->cuts.size();

  ImGui::InputText("cut id", g_cut_buf, IM_ARRAYSIZE(g_cut_buf));
  const int cut_id = std::atoi(g_cut_buf);

  // ---- toggles ----
  bool env2dChanged = ImGui::Checkbox("env 2D", &g_show_env_2d); ImGui::SameLine();
  bool env3dChanged = ImGui::Checkbox("env 3D", &g_show_env_3d);
  if (env2dChanged || env3dChanged) {
    envelope_vis::set_enabled(g_show_env_2d, g_show_env_3d);
  }

  bool cut2dChanged = ImGui::Checkbox("cut nodes 2D", &g_show_cut_nodes_2d); ImGui::SameLine();
  bool cut3dChanged = ImGui::Checkbox("cut nodes 3D", &g_show_cut_nodes_3d);

  if (cut2dChanged)
    for (auto* m : g_cut_nodes_2d) if (m) m->setEnabled(g_show_cut_nodes_2d);
  if (cut3dChanged)
    for (auto* m : g_cut_nodes_3d) if (m) m->setEnabled(g_show_cut_nodes_3d);

  // ---- info ----
  if (!(0 <= cut_id && cut_id < C)) {
    ImGui::Text("Valid cut indices: 0..%d", std::max(0, C - 1));
  } else {
    ImGui::Text("min(F) size: %d", (int)comp->cuts[cut_id].minimal_elements_F.size());
  }

  // ---- compute ----
  if (ImGui::Button("Compute min envelope of min(F)")) {

    if (0 <= cut_id && cut_id < C) {

      const auto& nodes = comp->cuts[cut_id].minimal_elements_F;

      // --------------------------------------------------------
      // Build cut-node meshes using ONE shared function
      // --------------------------------------------------------
      clear_cut_node_meshes();

      const std::string pref =
          std::string(g_completionWhich == 1 ? "P1" : "P2") +
          " Cut " + std::to_string(cut_id);

      for (int node_id : nodes) {
        // same scaling as compare_nodes
        auto pair = envelope_vis::build_mesh_pair_from_P1_node(
            pref + " node " + std::to_string(node_id) + " 2D",
            pref + " node " + std::to_string(node_id) + " 3D",
            g_in, g_P1, node_id,
            /*lift_min*/0.0f, /*lift_max*/2.0f, /*scale*/1.0f,
            /*edge_width*/1.0f,
            /*add_global_id_quantity*/false);

        if (pair.m2) {
          pair.m2->setEnabled(g_show_cut_nodes_2d);
          g_cut_nodes_2d.push_back(pair.m2);
        }
        if (pair.m3) {
          pair.m3->setEnabled(g_show_cut_nodes_3d);
          g_cut_nodes_3d.push_back(pair.m3);
        }
      }

      // --------------------------------------------------------
      // Compute envelope (shared pipeline)
      // --------------------------------------------------------
      envelope_vis::Options opt;
      opt.show_2d = g_show_env_2d;
      opt.show_3d = g_show_env_3d;

      // match compare_nodes scaling
      opt.scale = 1.0f;
      opt.lift_min = 0.0f;
      opt.lift_max = 2.0f;

      opt.edge_width = 1.0f;
      opt.reserve_tris = 4000;
      opt.winner_quantity_name = "winner node";

      // node ids from P1 and P2 are the same, just cut IDs differ
      const bool ok = envelope_vis::compute_and_register_from_P1_nodes(
          "envelope 2D",
          "envelope 3D",
          g_in, g_P1,
          nodes,
          opt);

      if (!ok) {
        std::cout << "[dedekind_envelope_ui] envelope computation failed\n";
      } else {
        envelope_vis::set_enabled(g_show_env_2d, g_show_env_3d);
      }
    }
  }

  ImGui::SameLine();

  if (ImGui::Button("Clear")) {
    envelope_vis::clear();
    clear_cut_node_meshes();
    if (polyscope::hasCurveNetwork("envelope diagram")) {
      polyscope::removeStructure("envelope diagram");
    }
  }
}

} // namespace dm_completion


#include "compare_nodes.h"

#include "render_helpers.h"   // viz_helpers::...
#include "poset.h"
#include "poset2.h"
#include "input.h"
#include "visualization.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <imgui.h>
#include <glm/glm.hpp>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// ------------------------------------------------------------
// Globals (same pattern as your other UIs)
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

static void apply_enabled_toggles() {
  if (g_A2) g_A2->setEnabled(g_show_A_2d);
  if (g_A3) g_A3->setEnabled(g_show_A_3d);
  if (g_B2) g_B2->setEnabled(g_show_B_2d);
  if (g_B3) g_B3->setEnabled(g_show_B_3d);
}

static std::string prefix_for(int whichPoset) {
  return (whichPoset == 1) ? "P1 " : "P2 ";
}

static void build_one_mesh_pair_from_P1_history(int node_idx,
                                                polyscope::SurfaceMesh*& out2d,
                                                polyscope::SurfaceMesh*& out3d,
                                                const std::string& name2d,
                                                const std::string& name3d)
{
  if (!(0 <= node_idx && node_idx < (int)g_P1.nodes.size())) {
    std::cout << "[compare_nodes] ERROR: invalid node id " << node_idx << "\n";
    return;
  }

  // replay triangulation from P1 history
  df::Tri2 tri = g_in.tri_poset;
  pst::replay_history_poset(tri, g_P1.nodes[node_idx].history, g_in);

  // build mesh data
  auto ids      = viz::present_ids(tri);
  auto to_local = viz::make_local_index(ids);
  auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

  float S = 20.0f;
  auto V2 = viz_helpers::make_planar_poset_vertices(
      ids, g_in.points2d, 0.0f, 20.0f, S);

  auto V3 = viz_helpers::make_lifted_poset_vertices(
      ids, g_in.points2d, 0.0f, 20.0f, S, 0.5f * S);

  // register overlays (true overlay: no offsets)
  out2d = polyscope::registerSurfaceMesh(name2d, V2, faces);
  out3d = polyscope::registerSurfaceMesh(name3d, V3, faces);

  // vertex id quantity (what you already do)
  viz_helpers::add_global_id_quantity(out2d, ids);
  viz_helpers::add_global_id_quantity(out3d, ids);

  out2d->setEdgeWidth(1.0f);
  out3d->setEdgeWidth(1.0f);
}

} // namespace


// ============================================================
// UI
// ============================================================
void compare_nodes_ui() {

  ImGui::SetNextItemOpen(false, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("Compare", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  // Choose poset (same style as your other UIs)
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

  if (changed) apply_enabled_toggles();

  const int a = std::atoi(g_a_buf);
  const int b = std::atoi(g_b_buf);

  // Since you said P1 and P2 share node ids, validate against P1 node universe.
  const int n = (int)g_P1.nodes.size();
  auto valid_idx = [&](int v) { return 0 <= v && v < n; };

  const bool a_ok = valid_idx(a);
  const bool b_ok = valid_idx(b);

  if (!a_ok || !b_ok) {
    ImGui::Text("Valid node indices: 0..%d", std::max(0, n - 1));
  }

  // Buttons
  if (ImGui::Button("Build overlay")) {

    clear_compare_overlay_internal();

    if (!a_ok || !b_ok) {
      std::cout << "[compare_nodes] cannot build overlay (invalid ids)\n";
    } else {
      const std::string pref = prefix_for(which);

      build_one_mesh_pair_from_P1_history(
          a, g_A2, g_A3,
          pref + "compare A 2D",
          pref + "compare A lifted");

      build_one_mesh_pair_from_P1_history(
          b, g_B2, g_B3,
          pref + "compare B 2D",
          pref + "compare B lifted");

      apply_enabled_toggles();
    }
  }

  ImGui::SameLine();

  if (ImGui::Button("Clear overlay")) {
    clear_compare_overlay_internal();
  }
}

} // namespace viz_poset





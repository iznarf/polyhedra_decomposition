#include "dedekind_completion.h"

#include "poset_view.h"
#include "dedekind_cut.h"
#include "dedekind_poset.h"

#include "poset_vis.h"

#include <imgui.h>
#include <glm/glm.hpp>

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>

#include <string>
#include <vector>
#include <iostream>

extern pst::Poset1  g_P1;
extern pst2::Poset2 g_P2;

// ------------------------------------------------------------
// store completion results for BOTH posets
// ------------------------------------------------------------
dm_completion::DedekindPoset g_completion_P1;
dm_completion::DedekindPoset g_completion_P2;
bool g_hasCompletion_P1 = false;
bool g_hasCompletion_P2 = false;

// status strings (separate so switching doesn't delete the other one)
static std::string g_status_P1;
static std::string g_status_P2;

// which completion are we currently VISUALIZING in this UI?
// 0 = P1, 1 = P2
static int  g_selectedPoset = 0;
static bool g_showCompletionNodes = false;
static bool g_showCompletionEdges = false;

// match other poset spacing defaults
static float g_compXSpacing = 0.2f;
static float g_compZSpacing = 0.2f;

// fixed center for completion visualization
static constexpr glm::vec3 center = glm::vec3(0.f, 0.f, 0.f);

// ------------------------------
// views
// ------------------------------
static PosetView view_of_completion(const dm_completion::DedekindPoset& D) {
  PosetView V;
  V.n = (int)D.cuts.size();
  V.cover_up          = &D.cover_up;
  V.cover_down        = &D.cover_down;
  V.topo_up           = &D.topo_up;
  V.topo_down         = &D.topo_down;
  V.reachability_up   = &D.reachability_up;
  V.reachability_down = &D.reachability_down;
  V.levels            = &D.levels;
  return V;
}

static PosetView current_input_view() {
  return (g_selectedPoset == 0) ? view_of(g_P1) : view_of(g_P2);
}

static dm_completion::DedekindPoset& current_completion() {
  return (g_selectedPoset == 0) ? g_completion_P1 : g_completion_P2;
}

static const dm_completion::DedekindPoset& current_completion_const() {
  return (g_selectedPoset == 0) ? g_completion_P1 : g_completion_P2;
}

static bool& current_hasCompletion_flag() {
  return (g_selectedPoset == 0) ? g_hasCompletion_P1 : g_hasCompletion_P2;
}

static std::string& current_status_string() {
  return (g_selectedPoset == 0) ? g_status_P1 : g_status_P2;
}

static std::string completion_nodes_name() {
  return (g_selectedPoset == 0) ? "P1 completion nodes" : "P2 completion nodes";
}
static std::string completion_edges_name() {
  return (g_selectedPoset == 0) ? "P1 completion edges" : "P2 completion edges";
}

// ------------------------------
// cleanup
// ------------------------------
static void remove_completion_vis() {
  const char* nodeNames[] = {"P1 completion nodes", "P2 completion nodes"};
  const char* edgeNames[] = {"P1 completion edges", "P2 completion edges"};
  for (auto* n : nodeNames) if (polyscope::hasPointCloud(n))   polyscope::removeStructure(n);
  for (auto* e : edgeNames) if (polyscope::hasCurveNetwork(e)) polyscope::removeStructure(e);
}

// ------------------------------
// build vis for CURRENT selected poset
// ------------------------------
static void build_completion_vis_for_selected() {

  const bool has = (g_selectedPoset == 0) ? g_hasCompletion_P1 : g_hasCompletion_P2;
  if (!has) return;

  const auto& D = current_completion_const();
  const int M = (int)D.cuts.size();
  if (M <= 0) return;

  // OLD color matches selected input poset:
  //   P1 -> red, P2 -> green
  const glm::vec3 oldColor = (g_selectedPoset == 0)
      ? glm::vec3(1.f, 0.f, 0.f)   // P1 red
      : glm::vec3(13.0f / 255.0f, 100.0f / 255.0f, 13.0f / 255.0f);  // P2 green

  // NEW completion nodes always blue
  const glm::vec3 newColor = glm::vec3(0.f, 0.f, 1.f);

  // build per-node colors using rule: |max(I')| == 1  <=> old/original
  std::vector<glm::vec3> colors;
  colors.reserve((size_t)M);

  for (int i = 0; i < M; ++i) {
    const auto& maxIp = D.cuts[i].maximal_elements_Iprime;
    const bool isOld = ((int)maxIp.size() == 1);
    colors.push_back(isOld ? oldColor : newColor);
  }

  const std::string node_names = completion_nodes_name();
  const std::string edge_names = completion_edges_name();

  // nodes (grid)
  std::vector<glm::vec3> pos =
      pst_vis::register_poset_as_grid_colored(
          node_names,
          M,
          D.levels,
          center,
          g_compXSpacing,
          g_compZSpacing,
          colors,
          0.01f);

  // edges in black
  pst_vis::register_cover_up_as_edges(edge_names, pos, D.cover_up, glm::vec3(0.f, 0.f, 0.f));

  // toggles
  if (polyscope::hasPointCloud(node_names))
    polyscope::getPointCloud(node_names)->setEnabled(g_showCompletionNodes);
  if (polyscope::hasCurveNetwork(edge_names))
    polyscope::getCurveNetwork(edge_names)->setEnabled(g_showCompletionEdges);
}

// ============================================================
// UI
// ============================================================
void dedekind_completion_ui() {

  if (!ImGui::CollapsingHeader("Completion", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  ImGui::SeparatorText("Dedekind completion");

  ImGui::TextUnformatted("Input poset:");
  int beforeSel = g_selectedPoset;
  ImGui::RadioButton("P1", &g_selectedPoset, 0);
  ImGui::SameLine();
  ImGui::RadioButton("P2", &g_selectedPoset, 1);
  ImGui::Spacing();

  // switching which completion you're looking at:
  // do NOT delete the other completion anymore — only hide visuals
  if (beforeSel != g_selectedPoset) {
    g_showCompletionNodes = false;
    g_showCompletionEdges = false;
    remove_completion_vis(); // remove both, we'll rebuild selected when needed
  }

  if (ImGui::Button("Compute completion")) {

    std::string& status = current_status_string();
    status.clear();

    PosetView P = current_input_view();

    std::vector<dm_completion::DedekindCut> cuts = dm_completion::compute_all_cuts(P);

    dm_completion::DedekindPoset& D = current_completion();
    D = dm_completion::build_dedekind_poset(P, std::move(cuts));

    bool& hasFlag = current_hasCompletion_flag();
    hasFlag = true;

    PosetView PC = view_of_completion(D);
    bool ok = dm_completion::check_complete_lattice(PC, true);

    status += "Completion nodes: " + std::to_string((int)D.cuts.size()) + "\n";
    status += ok ? "OK: completion is a lattice (finite => complete lattice)\n"
                 : "FAIL: completion is NOT a lattice\n";
    status += "Max level: " + std::to_string(D.maxLevel) + "\n";
  }

  // show the status of the currently selected completion
  std::string& status = current_status_string();
  if (!status.empty()) {
    ImGui::Separator();
    ImGui::TextUnformatted(status.c_str());
  }

  const bool hasSelected = (g_selectedPoset == 0) ? g_hasCompletion_P1 : g_hasCompletion_P2;
  if (!hasSelected) return;

  ImGui::SeparatorText("Completion visualization");

  /*
  ImGui::SliderFloat("x spacing", &g_compXSpacing, 0.02f, 1.0f, "%.3f");
  ImGui::SliderFloat("z spacing", &g_compZSpacing, 0.02f, 1.0f, "%.3f");
  */

  if (ImGui::Button("Show completion")) {
    remove_completion_vis(); // remove both old ones to avoid duplicates
    g_showCompletionNodes = true;
    g_showCompletionEdges = true;
    build_completion_vis_for_selected();
  }

  ImGui::Checkbox("Show completion nodes", &g_showCompletionNodes);
  ImGui::SameLine();
  ImGui::Checkbox("Show completion edges", &g_showCompletionEdges);

  // live-enable/disable without rebuild
  const std::string node_names = completion_nodes_name();
  const std::string edge_names = completion_edges_name();
  if (polyscope::hasPointCloud(node_names))
    polyscope::getPointCloud(node_names)->setEnabled(g_showCompletionNodes);
  if (polyscope::hasCurveNetwork(edge_names))
    polyscope::getCurveNetwork(edge_names)->setEnabled(g_showCompletionEdges);
}


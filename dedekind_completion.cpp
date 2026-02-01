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

static dm_completion::DedekindPoset g_completion;
static bool g_hasCompletion = false;
static std::string g_status;

static int  g_selectedPoset = 0; // 0=P1, 1=P2
static bool g_showCompletionNodes = false;
static bool g_showCompletionEdges = false;

// match your other poset spacing defaults
static float g_compXSpacing = 0.2f;
static float g_compZSpacing = 0.2f;

// fixed center for completion visualization
// shifted to the right of original posets
static constexpr glm::vec3 kCenter = glm::vec3(2.5f, 0.f, 0.f);

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

static std::string completion_nodes_name() {
  return (g_selectedPoset == 0) ? "P1 completion nodes" : "P2 completion nodes";
}
static std::string completion_edges_name() {
  return (g_selectedPoset == 0) ? "P1 completion edges" : "P2 completion edges";
}

static void remove_completion_vis() {
    const char* nodeNames[] = {"P1 completion nodes", "P2 completion nodes"};
    const char* edgeNames[] = {"P1 completion edges", "P2 completion edges"};
    for (auto* n : nodeNames) if (polyscope::hasPointCloud(n))   polyscope::removeStructure(n);
    for (auto* e : edgeNames) if (polyscope::hasCurveNetwork(e)) polyscope::removeStructure(e);
}

static void build_completion_vis() {
  if (!g_hasCompletion) return;

  const int M = (int)g_completion.cuts.size();
  if (M <= 0) return;

  // OLD color matches selected input poset:
  //   P1 -> red, P2 -> green
  const glm::vec3 oldColor = (g_selectedPoset == 0)
      ? glm::vec3(1.f, 0.f, 0.f)   // P1 red
      : glm::vec3(0.f, 1.f, 0.f);  // P2 green

  // NEW completion nodes always blue
  const glm::vec3 newColor = glm::vec3(0.f, 0.f, 1.f);

  // build per-node colors using rule: |max(I')| == 1  <=> old/original
  std::vector<glm::vec3> colors;
  colors.reserve((size_t)M);

  for (int i = 0; i < M; ++i) {
    const auto& maxIp = g_completion.cuts[i].maximal_elements_Iprime;
    const bool isOld = ((int)maxIp.size() == 1);
    colors.push_back(isOld ? oldColor : newColor);
  }

  // nodes (grid)
  const std::string nn = completion_nodes_name();
  const std::string en = completion_edges_name();

  std::vector<glm::vec3> pos =
      pst_vis::register_poset_as_grid_colored(
          nn,
          M,
          g_completion.levels,
          kCenter,
          g_compXSpacing,
          g_compZSpacing,
          colors,
          0.01f);

  // edges in black 
  pst_vis::register_cover_up_as_edges(en, pos, g_completion.cover_up, glm::vec3(0.f, 0.f, 0.f));

  // toggles
  if (polyscope::hasPointCloud(nn))   polyscope::getPointCloud(nn)->setEnabled(g_showCompletionNodes);
  if (polyscope::hasCurveNetwork(en)) polyscope::getCurveNetwork(en)->setEnabled(g_showCompletionEdges);
}

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

  // switching poset: drop old visuals + result
  if (beforeSel != g_selectedPoset) {
    g_hasCompletion = false;
    g_status.clear();
    g_showCompletionNodes = false;
    g_showCompletionEdges = false;
    remove_completion_vis();
  }

  if (ImGui::Button("Compute completion")) {

    g_status.clear();

    PosetView P = current_input_view();
    std::vector<dm_completion::DedekindCut> cuts = dm_completion::compute_all_cuts(P);
    g_completion = dm_completion::build_dedekind_poset(P, std::move(cuts));
    g_hasCompletion = true;

    PosetView PC = view_of_completion(g_completion);
    bool ok = dm_completion::check_complete_lattice(PC, true);

    g_status += "Completion nodes: " + std::to_string((int)g_completion.cuts.size()) + "\n";
    g_status += ok ? "OK: completion is a lattice (finite => complete lattice)\n"
                   : "FAIL: completion is NOT a lattice\n";
    g_status += "Max level: " + std::to_string(g_completion.maxLevel) + "\n";
  }

  if (!g_status.empty()) {
    ImGui::Separator();
    ImGui::TextUnformatted(g_status.c_str());
  }

  if (!g_hasCompletion) return;

  ImGui::SeparatorText("Completion visualization");

  ImGui::SliderFloat("x spacing", &g_compXSpacing, 0.02f, 1.0f, "%.3f");
  ImGui::SliderFloat("z spacing", &g_compZSpacing, 0.02f, 1.0f, "%.3f");

  if (ImGui::Button("Show completion")) {
    remove_completion_vis();
    // if you clicked build, enable by default
    g_showCompletionNodes = true;
    g_showCompletionEdges = true;
    build_completion_vis();
  }

    ImGui::Checkbox("Show completion nodes", &g_showCompletionNodes);
    ImGui::SameLine();
    ImGui::Checkbox("Show completion edges", &g_showCompletionEdges);

  // live-enable/disable without rebuild
  const std::string nn = completion_nodes_name();
  const std::string en = completion_edges_name();
  if (polyscope::hasPointCloud(nn))
    polyscope::getPointCloud(nn)->setEnabled(g_showCompletionNodes);
  if (polyscope::hasCurveNetwork(en))
    polyscope::getCurveNetwork(en)->setEnabled(g_showCompletionEdges);
}

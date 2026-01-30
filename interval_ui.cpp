#include "interval_ui.h"

#include "poset.h"
#include "poset2.h"
#include "poset_utils.h"
#include "poset_vis.h"

#include "poset_view.h"
#include "poset_interval.h"

#include <imgui.h>
#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>

#include <glm/glm.hpp>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <utility>

extern pst::Poset1  g_P1;
extern pst2::Poset2 g_P2;
extern bool g_has_P1;
extern bool g_has_P2;

namespace viz_poset {

namespace {

// ------------------------------------------------------------
// Overlay helpers
// ------------------------------------------------------------

static void clear_interval_overlay(int whichPoset /*1 or 2*/) {
  const std::string base = (whichPoset == 1) ? "P1 interval" : "P2 interval";

  if (polyscope::hasPointCloud(base + " nodes"))
    polyscope::removeStructure(base + " nodes");
  if (polyscope::hasCurveNetwork(base + " edges"))
    polyscope::removeStructure(base + " edges");
}

// Overlay interval nodes+edges using existing grid positions (global indexing)
static void register_interval_overlay(int whichPoset /*1 or 2*/,
                                      const PosetView& V,
                                      const pst::Bitset& mask,
                                      glm::vec3 colorNodes,
                                      glm::vec3 colorEdges) {
  // Convert mask -> list of global node ids
  std::vector<int> nodes = pst_interval::nodes_from_bitset(mask);

  // Induced cover edges (global ids)
  std::vector<std::pair<int,int>> edgesGlobal = pst_interval::induced_cover_up_edges(V, mask);

  // Need global positions that match node ids
  const auto& allPos = (whichPoset == 1) ? pst_vis::g_gridPos_P1 : pst_vis::g_gridPos_P2;

  if ((int)allPos.size() != V.n) {
    std::cout << "[interval_ui] ERROR: grid positions not computed for this poset.\n"
              << "             Click 'Visualize poset 1/2' first.\n";
    return;
  }

  // Build sub positions in the same order as `nodes`
  std::vector<glm::vec3> subPos;
  subPos.reserve(nodes.size());
  for (int v : nodes) subPos.push_back(allPos[v]);

  // global -> local map
  std::vector<int> idx(V.n, -1);
  for (int i = 0; i < (int)nodes.size(); ++i) idx[nodes[i]] = i;

  // Local edge list for curve network
  std::vector<std::array<int,2>> edgesLocal;
  edgesLocal.reserve(edgesGlobal.size());
  for (auto [u, v] : edgesGlobal) {
    int iu = (0 <= u && u < V.n) ? idx[u] : -1;
    int iv = (0 <= v && v < V.n) ? idx[v] : -1;
    if (iu >= 0 && iv >= 0) edgesLocal.push_back({iu, iv});
  }

  const std::string base = (whichPoset == 1) ? "P1 interval" : "P2 interval";

  // Replace old overlay if present
  clear_interval_overlay(whichPoset);

  // Nodes overlay
  auto* pc = polyscope::registerPointCloud(base + " nodes", subPos);
  pc->setEnabled(true);
  pc->setPointRadius(0.016f, false);   // bigger than base
  pc->setPointColor(colorNodes);

  // Edges overlay
  if (!edgesLocal.empty()) {
    auto* cn = polyscope::registerCurveNetwork(base + " edges", subPos, edgesLocal);
    cn->setEnabled(true);
    cn->setRadius(0.0040f, false);     // thicker than base
    cn->setColor(colorEdges);
  }
}

} // namespace


// ------------------------------------------------------------
// Interval UI
// ------------------------------------------------------------
void interval_ui() {

  ImGui::SetNextItemOpen(true, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("Interval", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  // Choose which poset we target
  static int which = 2; // default P2
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

  // Build view + levels pointer
  PosetView V;
  const std::vector<int>* levels = nullptr;

  if (which == 1) {
    V = view_of(g_P1);
    levels = &g_P1.levels;
  } else {
    V = view_of(g_P2);
    levels = &g_P2.levels;
  }

  const int N = V.n;
  if (!levels || (int)levels->size() != N) {
    ImGui::TextUnformatted("Poset levels missing or invalid.");
    return;
  }

  // base grid must have been visualized already (positions stored)
  bool baseReady = false;
  if (which == 1) baseReady = ((int)pst_vis::g_gridPos_P1.size() == N);
  else            baseReady = ((int)pst_vis::g_gridPos_P2.size() == N);

  if (!baseReady) {
    ImGui::TextUnformatted("Grid not ready: click 'Visualize poset 1/2' first.");
  }


  // Inputs
  static int  x = 0, y = 0;
  static char xbuf[32] = "0";
  static char ybuf[32] = "0";

  ImGui::InputText("x", xbuf, IM_ARRAYSIZE(xbuf));
  ImGui::InputText("y", ybuf, IM_ARRAYSIZE(ybuf));
  x = std::atoi(xbuf);
  y = std::atoi(ybuf);

  auto valid = [&](int v) { return 0 <= v && v < N; };
  bool x_ok = valid(x);
  bool y_ok = valid(y);

  if (!x_ok || !y_ok) {
    ImGui::Text("Valid indices: 0 .. %d", N - 1);
  }

  if (ImGui::Button("swap x/y")) {
    std::swap(x, y);
    std::snprintf(xbuf, IM_ARRAYSIZE(xbuf), "%d", x);
    std::snprintf(ybuf, IM_ARRAYSIZE(ybuf), "%d", y);
    x_ok = valid(x);
    y_ok = valid(y);
  }

  ImGui::SameLine();

  if (ImGui::Button("clear overlay")) {
    clear_interval_overlay(which);
  }

  ImGui::Separator();

    // Overlay colors (pick something that is NOT the base color)
  // base: P1 is red, P2 is green in your poset_vis_ui
  glm::vec3 overlayNodes = (which == 1) ? glm::vec3(0.2f, 0.8f, 1.0f) : glm::vec3(1.0f, 0.2f, 0.8f);
  glm::vec3 overlayEdges = overlayNodes;

  // Disable interval buttons if the base grid isn't ready
  ImGui::BeginDisabled(!baseReady);

  if (ImGui::Button("visualize interval [x,y]")) {

    if (!x_ok || !y_ok) {
      std::cout << "[interval_ui] invalid x/y\n";
    } else {
      pst::Bitset mask = pst_interval::interval_xy(V, x, y);
      auto nodes = pst_interval::nodes_from_bitset(mask);

      if (nodes.empty()) {
        std::cout << "[interval_ui] interval [" << x << "," << y << "] is empty (likely x !<= y)\n";
        clear_interval_overlay(which);
      } else {
        register_interval_overlay(which, V, mask, overlayNodes, overlayEdges);
      }
    }
  }

  ImGui::SameLine();

  if (ImGui::Button("visualize down(x)")) {
    if (!x_ok) {
      std::cout << "[interval_ui] invalid x\n";
    } else {
      pst::Bitset mask = pst_interval::lower_set(V, x);
      auto nodes = pst_interval::nodes_from_bitset(mask);

      if (nodes.empty()) {
        std::cout << "[interval_ui] down(" << x << ") is empty\n";
        clear_interval_overlay(which);
      } else {
        register_interval_overlay(which, V, mask, overlayNodes, overlayEdges);
      }
    }
  }

  ImGui::SameLine();

  if (ImGui::Button("visualize up(x)")) {
    if (!x_ok) {
      std::cout << "[interval_ui] invalid x\n";
    } else {
      pst::Bitset mask = pst_interval::upper_set(V, x);
      auto nodes = pst_interval::nodes_from_bitset(mask);

      if (nodes.empty()) {
        std::cout << "[interval_ui] up(" << x << ") is empty\n";
        clear_interval_overlay(which);
      } else {
        register_interval_overlay(which, V, mask, overlayNodes, overlayEdges);
      }
    }
  }

  ImGui::EndDisabled();


  // Status
  ImGui::Separator();
  ImGui::Text("poset: %s    nodes: %d", (which == 1 ? "P1" : "P2"), N);

  const std::string base = (which == 1) ? "P1 interval" : "P2 interval";
  ImGui::Text("overlay: %s / %s",
              polyscope::hasPointCloud(base + " nodes") ? "nodes ON" : "nodes OFF",
              polyscope::hasCurveNetwork(base + " edges") ? "edges ON" : "edges OFF");
}

} // namespace viz_poset


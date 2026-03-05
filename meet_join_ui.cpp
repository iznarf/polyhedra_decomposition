#include "meet_join_ui.h"

#include "poset.h"
#include "poset2.h"
#include "poset_view.h"
#include "poset_interval.h"
#include "poset_vis.h"
#include "meet_join.h"

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
extern pst::Poset_just_flips g_P1_flips;
extern bool g_has_P1;
extern bool g_has_P2;
extern bool g_has_P1_flips;

namespace viz_poset {

namespace {

// ---------- Overlay names ----------
static std::string overlayBase(int whichPoset) {
  if (whichPoset == 1) return "P1 meet/join";
  if (whichPoset == 2) return "P2 meet/join";
  return "P1 flips meet/join"; // whichPoset == 3
}

// ---------- Helpers ----------
static void clear_meet_join_overlay(int whichPoset) {
  const std::string base = overlayBase(whichPoset);

  if (polyscope::hasPointCloud(base + " inputs"))
    polyscope::removeStructure(base + " inputs");
  if (polyscope::hasPointCloud(base + " set nodes"))
    polyscope::removeStructure(base + " set nodes");
  if (polyscope::hasCurveNetwork(base + " set edges"))
    polyscope::removeStructure(base + " set edges");

  // NEW: candidate overlays
  if (polyscope::hasPointCloud(base + " meet candidates"))
    polyscope::removeStructure(base + " meet candidates");
  if (polyscope::hasPointCloud(base + " join candidates"))
    polyscope::removeStructure(base + " join candidates");
}

// Build an induced subgraph overlay (nodes+cover edges inside mask)
static void register_set_overlay(int whichPoset,
                                 const PosetView& V,
                                 const pst::Bitset& mask,
                                 glm::vec3 setColor)
{
  // Convert mask -> list of global node ids
  std::vector<int> nodes = pst_interval::nodes_from_bitset(mask);

  // Cover edges induced inside mask (global ids)
  std::vector<std::pair<int,int>> edgesGlobal = pst_interval::induced_cover_up_edges(V, mask);

  const auto& allPos =
    (whichPoset == 1) ? pst_vis::g_gridPos_P1 :
    (whichPoset == 2) ? pst_vis::g_gridPos_P2 :
                        pst_vis::g_gridPos_P1_flips;
  if ((int)allPos.size() != V.n) {
    std::cout << "[meet_join_ui] ERROR: grid positions not ready. Visualize poset first.\n";
    return;
  }

  // global -> local map
  std::vector<int> idx(V.n, -1);

  std::vector<glm::vec3> subPos;
  subPos.reserve(nodes.size());
  for (int i = 0; i < (int)nodes.size(); ++i) {
    int gid = nodes[i];
    idx[gid] = i;
    subPos.push_back(allPos[gid]);
  }

  std::vector<std::array<int,2>> edgesLocal;
  edgesLocal.reserve(edgesGlobal.size());
  for (auto [u, v] : edgesGlobal) {
    int iu = (0 <= u && u < V.n) ? idx[u] : -1;
    int iv = (0 <= v && v < V.n) ? idx[v] : -1;
    if (iu >= 0 && iv >= 0) edgesLocal.push_back({iu, iv});
  }

  const std::string base = overlayBase(whichPoset);

  // Replace if exists
  if (polyscope::hasPointCloud(base + " set nodes"))
    polyscope::removeStructure(base + " set nodes");
  if (polyscope::hasCurveNetwork(base + " set edges"))
    polyscope::removeStructure(base + " set edges");

  // nodes of the set
  auto* pc = polyscope::registerPointCloud(base + " set nodes", subPos);
  pc->setEnabled(true);
  pc->setPointRadius(0.016f, false);
  pc->setPointColor(setColor);

  // edges of the induced cover graph
  if (!edgesLocal.empty()) {
    auto* cn = polyscope::registerCurveNetwork(base + " set edges", subPos, edgesLocal);
    cn->setEnabled(true);
    cn->setRadius(0.0040f, false);
    cn->setColor(setColor);
  }
}

// Highlight the two input nodes as their own tiny point cloud (yellow)
static void register_inputs_overlay(int whichPoset,
                                    const PosetView& V,
                                    int x, int y,
                                    glm::vec3 colorInputs)
{
  const auto& allPos =
    (whichPoset == 1) ? pst_vis::g_gridPos_P1 :
    (whichPoset == 2) ? pst_vis::g_gridPos_P2 :
                        pst_vis::g_gridPos_P1_flips;
  if ((int)allPos.size() != V.n) return;

  std::vector<glm::vec3> pts;
  pts.reserve(2);
  pts.push_back(allPos[x]);
  if (y != x) pts.push_back(allPos[y]);

  const std::string base = overlayBase(whichPoset);

  if (polyscope::hasPointCloud(base + " inputs"))
    polyscope::removeStructure(base + " inputs");

  auto* pc = polyscope::registerPointCloud(base + " inputs", pts);
  pc->setEnabled(true);
  pc->setPointRadius(0.022f, false); // slightly bigger than set nodes
  pc->setPointColor(colorInputs);
}

// NEW: highlight candidate nodes (pink)
static void register_candidates_overlay(int whichPoset,
                                        const PosetView& V,
                                        const std::vector<int>& cand,
                                        const std::string& suffix,
                                        glm::vec3 color)
{
  const auto& allPos =
    (whichPoset == 1) ? pst_vis::g_gridPos_P1 :
    (whichPoset == 2) ? pst_vis::g_gridPos_P2 :
                        pst_vis::g_gridPos_P1_flips;
  if ((int)allPos.size() != V.n) return;

  std::vector<glm::vec3> pts;
  pts.reserve(cand.size());
  for (int id : cand) {
    if (0 <= id && id < V.n) pts.push_back(allPos[id]);
  }

  const std::string base = overlayBase(whichPoset);
  const std::string name = base + " " + suffix;

  if (polyscope::hasPointCloud(name))
    polyscope::removeStructure(name);

  auto* pc = polyscope::registerPointCloud(name, pts);
  pc->setEnabled(true);
  pc->setPointRadius(0.028f, false); // biggest
  pc->setPointColor(color);
}

// Common upper bounds of {x,y} = Up(x) ∩ Up(y)
static pst::Bitset common_upper_bounds_mask(const PosetView& V, int x, int y) {
  pst::Bitset Ux = pst_interval::upper_set(V, x);
  pst::Bitset Uy = pst_interval::upper_set(V, y);
  Ux &= Uy;
  return Ux;
}

// Common lower bounds of {x,y} = Down(x) ∩ Down(y)
static pst::Bitset common_lower_bounds_mask(const PosetView& V, int x, int y) {
  pst::Bitset Dx = pst_interval::lower_set(V, x);
  pst::Bitset Dy = pst_interval::lower_set(V, y);
  Dx &= Dy;
  return Dx;
}

// Helper to show ids in one line
static void imgui_ids_inline(const std::vector<int>& ids) {
  for (int i = 0; i < (int)ids.size(); ++i) {
    if (i > 0) {
      ImGui::SameLine();
      ImGui::TextUnformatted(",");
      ImGui::SameLine();
    }
    ImGui::Text("%d", ids[i]);
  }
}

} // namespace


void meet_join_ui() {

  ImGui::SetNextItemOpen(true, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("Meet & Join", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  // Choose poset
  static int which = 2;
  if (!g_has_P2) which = 1;

  const bool canP1 = g_has_P1;
  const bool canP2 = g_has_P2;
  const bool canP1flips = g_has_P1_flips;

  ImGui::TextUnformatted("Target poset:");

  ImGui::BeginDisabled(!canP1);
  ImGui::RadioButton("P1", &which, 1);
  ImGui::EndDisabled();

  ImGui::SameLine();

  ImGui::BeginDisabled(!canP2);
  ImGui::RadioButton("P2", &which, 2);
  ImGui::EndDisabled();

  ImGui::SameLine();

  ImGui::BeginDisabled(!canP1flips);
  ImGui::RadioButton("P1 just flips", &which, 3);
  ImGui::EndDisabled();

  if ((which == 1 && !canP1) || (which == 2 && !canP2) || (which == 3 && !canP1flips)) {
    ImGui::TextUnformatted("Selected poset not available.");
    return;
  }

  // view
  PosetView V;
  if (which == 1)      V = view_of(g_P1);
  else if (which == 2) V = view_of(g_P2);
  else                 V = view_of(g_P1_flips);

  const int N = V.n;

  // base ready?
  bool baseReady = false;
  if (which == 1)      baseReady = ((int)pst_vis::g_gridPos_P1.size() == N);
  else if (which == 2) baseReady = ((int)pst_vis::g_gridPos_P2.size() == N);
  else                 baseReady = ((int)pst_vis::g_gridPos_P1_flips.size() == N);

  if (!baseReady) {
    ImGui::TextUnformatted("Grid not ready: click 'Visualize poset 1/2' first.");
  }

  // Inputs x,y
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
    clear_meet_join_overlay(which);
  }

  ImGui::Separator();

  // Colors
  glm::vec3 colorInputs     = glm::vec3(1.0f, 1.0f, 0.0f); // yellow
  glm::vec3 colorSet        = glm::vec3(1.0f, 1.0f, 1.0f); // white
  glm::vec3 colorCandidates = glm::vec3(1.0f, 0.2f, 0.8f); // pink

  ImGui::BeginDisabled(!baseReady);

  if (ImGui::Button("common upper bounds & join candidates")) {
    if (!x_ok || !y_ok) {
      std::cout << "[meet_join_ui] invalid x/y\n";
    } else {
      clear_meet_join_overlay(which);

      pst::Bitset U = common_upper_bounds_mask(V, x, y);

      // Always show inputs
      register_inputs_overlay(which, V, x, y, colorInputs);

      // Show induced subgraph on the bound set
      auto nodes = pst_interval::nodes_from_bitset(U);
      if (nodes.empty()) {
        std::cout << "[meet_join_ui] common upper bounds empty for (" << x << "," << y << ")\n";
      } else {
        register_set_overlay(which, V, U, colorSet);
      }

      // NEW: show join candidates in pink (unique or not)
      auto j = meet_join::join(V, x, y);
      if (j.exists()) {
        register_candidates_overlay(which, V, j.candidates, "join candidates", colorCandidates);
      }
    }
  }

  ImGui::SameLine();

  if (ImGui::Button("common lower bounds & meet candidates")) {
    if (!x_ok || !y_ok) {
      std::cout << "[meet_join_ui] invalid x/y\n";
    } else {
      clear_meet_join_overlay(which);

      pst::Bitset D = common_lower_bounds_mask(V, x, y);

      register_inputs_overlay(which, V, x, y, colorInputs);

      auto nodes = pst_interval::nodes_from_bitset(D);
      if (nodes.empty()) {
        std::cout << "[meet_join_ui] common lower bounds empty for (" << x << "," << y << ")\n";
      } else {
        register_set_overlay(which, V, D, colorSet);
      }

      // NEW: show meet candidates in pink (unique or not)
      auto m = meet_join::meet(V, x, y);
      if (m.exists()) {
        register_candidates_overlay(which, V, m.candidates, "meet candidates", colorCandidates);
      }
    }
  }

  ImGui::EndDisabled();

  // Updated meet/join text: show ids even when not unique
  ImGui::Separator();
  if (x_ok && y_ok) {
    auto m = meet_join::meet(V, x, y);
    auto j = meet_join::join(V, x, y);

    // meet
    if (!m.exists()) {
      ImGui::TextUnformatted("meet: none");
    } else if (m.unique()) {
      ImGui::Text("meet: %d", m.candidates[0]);
    } else {
      ImGui::TextUnformatted("meet: (not unique)");
      ImGui::SameLine();
      imgui_ids_inline(m.candidates);
    }

    // join
    if (!j.exists()) {
      ImGui::TextUnformatted("join: none");
    } else if (j.unique()) {
      ImGui::Text("join: %d", j.candidates[0]);
    } else {
      ImGui::TextUnformatted("join: (not unique)");
      ImGui::SameLine();
      imgui_ids_inline(j.candidates);
    }
  }
}

} // namespace viz_poset


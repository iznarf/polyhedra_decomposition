#include "moebius_ui.h"

#include "poset.h"
#include "poset2.h"
#include "poset_view.h"
#include "poset_vis.h"
#include "moebius.h"

#include <imgui.h>
#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>

#include <glm/glm.hpp>

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

extern pst::Poset1  g_P1;
extern pst2::Poset2 g_P2;
extern bool g_has_P1;
extern bool g_has_P2;

namespace viz_poset {

namespace {

// ---------- Overlay name ----------
static std::string overlayName(int whichPoset) {
  return (whichPoset == 1) ? "P1 moebius" : "P2 moebius";
}

// ---------- Helpers ----------
static void clear_moebius_overlay(int whichPoset) {
  const std::string name = overlayName(whichPoset);
  if (polyscope::hasPointCloud(name))
    polyscope::removeStructure(name);
}

static bool base_grid_ready(int whichPoset, const PosetView& V) {
  if (whichPoset == 1) return ((int)pst_vis::g_gridPos_P1.size() == V.n);
  else                 return ((int)pst_vis::g_gridPos_P2.size() == V.n);
}

static glm::vec3 color_for_mu(std::int64_t mu) {
  // +1 green, 0 gray, -1 red, else blue
  if (mu == 1)  return glm::vec3(0.2f, 0.9f, 0.2f);
  if (mu == 0)  return glm::vec3(0.6f, 0.6f, 0.6f);
  if (mu == -1) return glm::vec3(0.9f, 0.2f, 0.2f);
  return glm::vec3(0.2f, 0.4f, 0.9f);
}

static void register_moebius_overlay(int whichPoset,
                                     const PosetView& V,
                                     const std::vector<std::int64_t>& mu)
{
  const auto& allPos = (whichPoset == 1) ? pst_vis::g_gridPos_P1
                                         : pst_vis::g_gridPos_P2;

  if ((int)allPos.size() != V.n) {
    std::cout << "[moebius_ui] ERROR: grid positions not ready. Visualize poset first.\n";
    return;
  }
  if ((int)mu.size() != V.n) {
    std::cout << "[moebius_ui] ERROR: mu size mismatch\n";
    return;
  }

  std::vector<glm::vec3> colors(V.n);
  for (int i = 0; i < V.n; ++i) colors[i] = color_for_mu(mu[i]);

  const std::string name = overlayName(whichPoset);

  if (polyscope::hasPointCloud(name))
    polyscope::removeStructure(name);

  // Overlay uses same positions as base grid
  auto* pc = polyscope::registerPointCloud(name, allPos);
  pc->setEnabled(true);
  pc->setPointRadius(0.020f, false);

  // Put colors in a quantity so Polyscope renders per-point colors
  pc->addColorQuantity("mu(min,x)", colors)->setEnabled(true);
}

} // namespace

void moebius_ui() {

  ImGui::SetNextItemOpen(true, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("Moebius", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  // Choose poset (same pattern as meet/join)
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

  // view
  PosetView V;
  if (which == 1) V = view_of(g_P1);
  else            V = view_of(g_P2);

  // base ready?
  bool ready = base_grid_ready(which, V);
  if (!ready) {
    ImGui::TextUnformatted("Grid not ready: click 'Visualize poset 1/2' first.");
  }

  // Legend text (as requested)
  ImGui::Separator();
  ImGui::TextUnformatted("Legend for moeb(0,x):");
  ImGui::TextColored(ImVec4(0.2f, 0.9f, 0.2f, 1.0f), " +1  -> green");
  ImGui::TextColored(ImVec4(0.6f, 0.6f, 0.6f, 1.0f), "  0  -> gray");
  ImGui::TextColored(ImVec4(0.9f, 0.2f, 0.2f, 1.0f), " -1  -> red");
  ImGui::TextColored(ImVec4(0.2f, 0.4f, 0.9f, 1.0f), " other -> blue");

  ImGui::Spacing();

  if (ImGui::Button("clear overlay")) {
    clear_moebius_overlay(which);
  }

  ImGui::SameLine();

  ImGui::BeginDisabled(!ready);
  if (ImGui::Button("compute moeb(0,x)")) {
    clear_moebius_overlay(which);
    auto mu = pst::mobius_from_global_min(V);
    register_moebius_overlay(which, V, mu);
  }
  ImGui::EndDisabled();
}

} // namespace viz_poset


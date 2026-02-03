#include "dedekind_completion.h"

#include "poset_view.h"
#include "dedekind_cut.h"
#include "dedekind_poset.h"

#include "poset_vis.h"
#include "completion_vis.h"   // <--- NEW

#include <imgui.h>
#include <glm/glm.hpp>

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>

#include <string>
#include <vector>
#include <iostream>

extern df::InputData g_in;   // <--- NEW: needed for completion_vis::build

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

// which completion are we currently COMPUTING in this UI?
// 0 = P1, 1 = P2
static int  g_selectedPoset = 0;

// match other poset spacing defaults
static float g_compXSpacing = 0.2f;
static float g_compZSpacing = 0.2f;

// fixed center for completion visualization
static constexpr glm::vec3 center = glm::vec3(0.f, 0.f, 0.f);

// ------------------------------------------------------------
// visualization toggles (separate for P1 completion and P2 completion)
// ------------------------------------------------------------
static bool g_vis_grid[2]  = {false, false};
static bool g_vis_edges[2] = {false, false};
static bool g_vis_2d[2]    = {false, false};

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

static bool& current_hasCompletion_flag() {
  return (g_selectedPoset == 0) ? g_hasCompletion_P1 : g_hasCompletion_P2;
}

static std::string& current_status_string() {
  return (g_selectedPoset == 0) ? g_status_P1 : g_status_P2;
}

// ============================================================
// UI
// ============================================================
void dedekind_completion_ui() {

  if (!ImGui::CollapsingHeader("Completion", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  ImGui::SeparatorText("Dedekind completion");

  ImGui::TextUnformatted("Input poset for computing completion:");
  ImGui::RadioButton("P1", &g_selectedPoset, 0);
  ImGui::SameLine();
  ImGui::RadioButton("P2", &g_selectedPoset, 1);
  ImGui::Spacing();

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

  // show both status blocks (nice when you compute both)
  if (!g_status_P1.empty() || !g_status_P2.empty()) {
    ImGui::SeparatorText("Status");

    if (!g_status_P1.empty()) {
      ImGui::TextUnformatted("P1 completion:");
      ImGui::TextUnformatted(g_status_P1.c_str());
    }
    if (!g_status_P2.empty()) {
      ImGui::TextUnformatted("P2 completion:");
      ImGui::TextUnformatted(g_status_P2.c_str());
    }
  }

  ImGui::SeparatorText("Completion visualization");

  // ------------------------------------------------------------------
  // P1 completion controls
  // ------------------------------------------------------------------
  ImGui::TextUnformatted("P1 completion:");
  ImGui::BeginDisabled(!g_hasCompletion_P1);
  if (ImGui::Button("Visualize completion of P1")) {
    // build everything for P1 completion
    completion_vis::build(
      0,
      g_in,
      g_P1,
      g_P2,
      g_completion_P1,
      center,
      g_compXSpacing,
      g_compZSpacing
    );

    // default enable grid+edges; keep 2D as user set
    g_vis_grid[0]  = true;
    g_vis_edges[0] = true;
    // g_vis_2d[0] unchanged
    completion_vis::set_enabled(0, g_vis_grid[0], g_vis_edges[0], g_vis_2d[0]);
  }
  ImGui::EndDisabled();

  bool changedP1 = false;
  ImGui::BeginDisabled(!g_hasCompletion_P1);
  changedP1 |= ImGui::Checkbox("P1 grid", &g_vis_grid[0]); ImGui::SameLine();
  changedP1 |= ImGui::Checkbox("P1 edges", &g_vis_edges[0]); ImGui::SameLine();
  changedP1 |= ImGui::Checkbox("P1 2D meshes", &g_vis_2d[0]);
  if (changedP1) completion_vis::set_enabled(0, g_vis_grid[0], g_vis_edges[0], g_vis_2d[0]);

  if (ImGui::Button("Clear P1 completion")) {
    completion_vis::clear(0);
  }
  ImGui::EndDisabled();

  ImGui::Separator();

  // ------------------------------------------------------------------
  // P2 completion controls
  // ------------------------------------------------------------------
  ImGui::TextUnformatted("P2 completion:");
  ImGui::BeginDisabled(!g_hasCompletion_P2);
  if (ImGui::Button("Visualize completion of P2")) {
    completion_vis::build(
      1,
      g_in,
      g_P1,
      g_P2,
      g_completion_P2,
      center,
      g_compXSpacing,
      g_compZSpacing
    );

    g_vis_grid[1]  = true;
    g_vis_edges[1] = true;
    completion_vis::set_enabled(1, g_vis_grid[1], g_vis_edges[1], g_vis_2d[1]);
  }
  ImGui::EndDisabled();

  bool changedP2 = false;
  ImGui::BeginDisabled(!g_hasCompletion_P2);
  changedP2 |= ImGui::Checkbox("P2 grid", &g_vis_grid[1]); ImGui::SameLine();
  changedP2 |= ImGui::Checkbox("P2 edges", &g_vis_edges[1]); ImGui::SameLine();
  changedP2 |= ImGui::Checkbox("P2 2D meshes", &g_vis_2d[1]);
  if (changedP2) completion_vis::set_enabled(1, g_vis_grid[1], g_vis_edges[1], g_vis_2d[1]);

  if (ImGui::Button("Clear P2 completion")) {
    completion_vis::clear(1);
  }
  ImGui::EndDisabled();

  // Optional: clear both
  ImGui::Separator();
  if (ImGui::Button("Clear ALL completion visuals")) {
    completion_vis::clear_all();
  }
}

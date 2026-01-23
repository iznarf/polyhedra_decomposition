#include "vis_poset.h"

#include "visualization.h"
#include "geometry_utils.h"
#include "poset_utils.h"
#include "render_helpers.h"
#include "node_coloring.h"
#include "interval_ui.h"
#include "meet_join_ui.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>

#include <glm/glm.hpp>
#include <glm/mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <imgui.h>

#include <vector>
#include <array>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <limits>

using glm::vec3;
using df::vertex_id;

// ================== GLOBALS (external linkage) ==================

std::vector<polyscope::SurfaceMesh*> g_poset_meshes_2d;
std::vector<polyscope::SurfaceMesh*> g_poset_meshes_3d;

bool g_show_poset_2d = true;
bool g_show_poset_3d = true;

std::vector<glm::vec3> g_node_centers;

polyscope::CurveNetwork* g_poset1_edges = nullptr;
bool g_show_poset1_edges = true;

pst::Poset1 g_poset1;

polyscope::CurveNetwork* g_poset2_edges = nullptr;
bool g_show_poset2_edges = true;

const df::InputData* g_input_for_poset = nullptr;

pst2::Poset2 g_poset2;

// interval/focus UI globals (as in your snippet)
int  g_interval_x = 0;
int  g_interval_y = 0;
bool g_interval_active = false;
std::vector<char> g_interval_mask;
char g_interval_x_buf[32] = "0";
char g_interval_y_buf[32] = "0";

char g_pair_a_buf[32] = "0";
char g_pair_b_buf[32] = "0";

bool g_focus_active = false;
std::vector<char> g_focus_mask;

char g_cmp_a_buf[32] = "0";
char g_cmp_b_buf[32] = "0";

bool g_cmp_show_A_2d = false;
bool g_cmp_show_A_3d = true;
bool g_cmp_show_B_2d = false;
bool g_cmp_show_B_3d = true;

bool g_compare_mode = false;

polyscope::SurfaceMesh* g_cmpA_2d = nullptr;
polyscope::SurfaceMesh* g_cmpB_2d = nullptr;
polyscope::SurfaceMesh* g_cmpA_3d = nullptr;
polyscope::SurfaceMesh* g_cmpB_3d = nullptr;

// ===== layout globals declared extern in header =====
namespace viz_poset {
PosetLayout g_layout_p1;
PosetLayout g_layout_p2;
ActiveLayout g_active_layout = ActiveLayout::P1;
} // namespace viz_poset

// ===================================================

namespace viz_poset {

// ------------------- internal helpers -------------------

static void remove_poset_meshes() {
  for (auto* m : g_poset_meshes_2d) if (m) m->remove();
  for (auto* m : g_poset_meshes_3d) if (m) m->remove();
  g_poset_meshes_2d.clear();
  g_poset_meshes_3d.clear();
}

static void remove_edge_networks() {
  if (g_poset1_edges) { g_poset1_edges->remove(); g_poset1_edges = nullptr; }
  if (g_poset2_edges) { g_poset2_edges->remove(); g_poset2_edges = nullptr; }
}

static PosetLayout build_layout_from_cover_down(
    const std::vector<std::vector<int>>& cover_down,
    int n,
    int root,
    float level_spacing,
    float node_spacing)
{
  PosetLayout L;
  L.level = pst::compute_levels_longest_from_root_cover_down(cover_down, root);

  L.max_level = 0;
  for (int lv : L.level) L.max_level = std::max(L.max_level, lv);

  L.byLevel.assign(L.max_level + 1, {});
  for (int i = 0; i < n; ++i) {
    int lv = (i < (int)L.level.size() ? L.level[i] : 0);
    if (lv < 0) lv = 0;
    L.byLevel[lv].push_back(i);
  }

  L.centers.assign(n, glm::vec3(0));

  for (int lv = 0; lv <= L.max_level; ++lv) {
    auto& vec = L.byLevel[lv];
    int k = (int)vec.size();
    float row_z = lv * level_spacing;

    for (int j = 0; j < k; ++j) {
      int node_idx = vec[j];
      float x = (j - 0.5f * (k - 1)) * node_spacing;
      L.centers[node_idx] = glm::vec3(x, 0.f, row_z);
    }
  }

  return L;
}

static const PosetLayout& active_layout() {
  return (g_active_layout == ActiveLayout::P1) ? g_layout_p1 : g_layout_p2;
}

// Move all meshes and update g_node_centers
static void apply_layout_to_meshes(const PosetLayout& L) {
  const int n = (int)g_poset1.nodes.size();
  if ((int)L.centers.size() != n) return;

  g_node_centers = L.centers;

  for (int i = 0; i < n; ++i) {
    glm::vec3 c = L.centers[i];
    glm::mat4 T = glm::translate(glm::mat4(1.f), c);

    if (i < (int)g_poset_meshes_2d.size() && g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setTransform(T);
    if (i < (int)g_poset_meshes_3d.size() && g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setTransform(T);
  }
}

static void rebuild_poset1_edges_filtered() {
  const int n = (int)g_poset1.nodes.size();
  if (n == 0) return;

  if (g_poset1_edges) { g_poset1_edges->remove(); g_poset1_edges = nullptr; }
  if ((int)g_node_centers.size() != n) return;

  std::vector<glm::uvec2> edges;

  size_t m = 0;
  for (int u = 0; u < n; ++u) if (node_visible(u)) m += g_poset1.cover_down[u].size();
  edges.reserve(m);

  for (int u = 0; u < n; ++u) {
    if (!node_visible(u)) continue;
    for (int v : g_poset1.cover_down[u]) {
      if (v < 0 || v >= n) continue;
      if (!node_visible(v)) continue;
      edges.emplace_back((unsigned)u, (unsigned)v);
    }
  }

  if (edges.empty()) return;

  g_poset1_edges = polyscope::registerCurveNetwork("poset1 edges", g_node_centers, edges);
  g_poset1_edges->setEnabled(g_show_poset1_edges);
  g_poset1_edges->setRadius(0.00037f, true);
  g_poset1_edges->setColor(glm::vec3(1.0f, 0.0f, 0.0f));
}

static void rebuild_poset2_edges_filtered() {
  const int n = (int)g_poset2.cover_down.size();
  if (n == 0) return;

  if (g_poset2_edges) { g_poset2_edges->remove(); g_poset2_edges = nullptr; }
  if ((int)g_node_centers.size() != n) return;

  std::vector<glm::uvec2> edges;

  size_t m = 0;
  for (int u = 0; u < n; ++u) if (node_visible(u)) m += g_poset2.cover_down[u].size();
  edges.reserve(m);

  for (int u = 0; u < n; ++u) {
    if (!node_visible(u)) continue;
    for (int v : g_poset2.cover_down[u]) {
      if (v < 0 || v >= n) continue;
      if (!node_visible(v)) continue;
      edges.emplace_back((unsigned)u, (unsigned)v);
    }
  }

  if (edges.empty()) return;

  g_poset2_edges = polyscope::registerCurveNetwork("poset2 edges", g_node_centers, edges);
  g_poset2_edges->setEnabled(g_show_poset2_edges);
  g_poset2_edges->setRadius(0.00037f, true);
  g_poset2_edges->setColor(glm::vec3(0.0f, 0.8f, 0.0f));
}

static void rebuild_both_edge_networks() {
  rebuild_poset1_edges_filtered();
  rebuild_poset2_edges_filtered();
}

// ------------------- API functions (declared in header) -------------------

bool node_visible(int idx) {
  const int n = (int)g_poset1.nodes.size();
  if (idx < 0 || idx >= n) return false;

  if (g_focus_active) {
    return idx < (int)g_focus_mask.size() && g_focus_mask[idx];
  }
  if (g_interval_active) {
    return idx < (int)g_interval_mask.size() && g_interval_mask[idx];
  }
  return true;
}

void set_poset_enabled(bool enabled) {
  for (auto* m : g_poset_meshes_2d) if (m) m->setEnabled(enabled && g_show_poset_2d);
  for (auto* m : g_poset_meshes_3d) if (m) m->setEnabled(enabled && g_show_poset_3d);

  if (g_poset1_edges) g_poset1_edges->setEnabled(enabled && g_show_poset1_edges);
  if (g_poset2_edges) g_poset2_edges->setEnabled(enabled && g_show_poset2_edges);
}

void apply_poset_visibility_from_masks() {
  const int n = (int)g_poset1.nodes.size();
  for (int i = 0; i < n; ++i) {
    bool on = node_visible(i);
    if (i < (int)g_poset_meshes_2d.size() && g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setEnabled(on && g_show_poset_2d);
    if (i < (int)g_poset_meshes_3d.size() && g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setEnabled(on && g_show_poset_3d);
  }
  rebuild_both_edge_networks();
}

// Keep names for compatibility with other code:
void rebuild_downflip_network_filtered() {
  // historically this was “poset1 edges”; keep it, but now we rebuild both because centers changed / masks changed
  rebuild_both_edge_networks();
}

void register_poset1_node_meshes(const df::InputData& D, const pst::Poset1& P1) {
    const int n = (int)P1.nodes.size();
    if (n == 0) return;

    g_poset_meshes_2d.assign(n, nullptr);
    g_poset_meshes_3d.assign(n, nullptr);

    // We register meshes around origin; layout is applied via transforms
    const float TRI_SCALE_PLAN  = 0.7f;
    const float TRI_SCALE_LIFT  = 0.7f;
    const float LIFT_HEIGHT_SCL = 0.5f;

    for (int node_idx = 0; node_idx < n; ++node_idx) {
      const auto& node = P1.nodes[node_idx];

      df::Tri2 tri = D.tri_poset;
      pst::replay_history_poset(tri, node.history, D);

      auto ids      = viz::present_ids(tri);
      auto to_local = viz::make_local_index(ids);
      auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

      auto V2 = viz_helpers::make_planar_poset_vertices(ids, D.points2d, 0.f, 0.f, TRI_SCALE_PLAN);
      auto V3 = viz_helpers::make_lifted_poset_vertices(ids, D.points2d, 0.f, 0.f, TRI_SCALE_LIFT, LIFT_HEIGHT_SCL);

      std::string name2d = "poset node " + std::to_string(node_idx) + " 2D";
      std::string name3d = "poset node " + std::to_string(node_idx) + " lifted";

      auto* m2 = polyscope::registerSurfaceMesh(name2d, V2, faces);
      viz_helpers::add_global_id_quantity(m2, ids);
      m2->setEnabled(g_show_poset_2d);
      m2->setSurfaceColor(glm::vec3(0.6f, 0.8f, 1.0f));
      m2->setEdgeWidth(1.0f);
      m2->setEdgeColor(glm::vec3(0, 0, 0));

      auto* m3 = polyscope::registerSurfaceMesh(name3d, V3, faces);
      viz_helpers::add_global_id_quantity(m3, ids);
      m3->setEnabled(g_show_poset_3d);
      m3->setSurfaceColor(glm::vec3(0.2f, 0.4f, 0.8f));
      m3->setTransparency(0.6f);
      m3->setEdgeWidth(1.0f);
      m3->setEdgeColor(glm::vec3(0, 0, 0));

      g_poset_meshes_2d[node_idx] = m2;
      g_poset_meshes_3d[node_idx] = m3;
    }
}

void register_poset1_edges(const pst::Poset1& /*P1*/) {
  rebuild_poset1_edges_filtered();
}

// Keep your old API, but now this stores P2 and rebuilds green edges.
void register_poset2_cover_edges(const pst2::Poset2& P2) {
    g_poset2 = P2;

    // optional: ensure cover_up exists (if used elsewhere)
    const int n = (int)g_poset2.cover_down.size();
    if ((int)g_poset2.cover_up.size() != n) {
      g_poset2.cover_up.assign(n, {});
      for (int u = 0; u < n; ++u) {
        for (int v : g_poset2.cover_down[u]) {
          if (v < 0 || v >= n) continue;
          g_poset2.cover_up[v].push_back(u);
        }
      }
      pst::sort_unique_adjacency(g_poset2.cover_up);
    }

    // recompute layout_p2 if we already have poset1 nodes
    if (!g_poset1.nodes.empty() && (int)g_poset1.nodes.size() == n) {
      const float LEVEL_SPACING = 4.f;
      const float NODE_SPACING  = 4.f;
      g_layout_p2 = build_layout_from_cover_down(g_poset2.cover_down, n, 0, LEVEL_SPACING, NODE_SPACING);

      // re-apply current layout (might be P2) and rebuild edges
      apply_layout_to_meshes(active_layout());
      rebuild_both_edge_networks();
    } else {
      rebuild_poset2_edges_filtered();
    }
}

void register_poset(const df::InputData& D, const pst::Poset1& P1) {
    remove_poset_meshes();
    remove_edge_networks();

    g_poset1 = P1;
    g_input_for_poset = &D;

    if (g_poset1.nodes.empty()) {
      std::cout << "[vis_poset] no nodes to visualize.\n";
      return;
    }

    // meshes
    register_poset1_node_meshes(D, g_poset1);

    // layouts
    const int n = (int)g_poset1.nodes.size();
    const float LEVEL_SPACING = 4.f;
    const float NODE_SPACING  = 4.f;

    g_layout_p1 = build_layout_from_cover_down(g_poset1.cover_down, n, 0, LEVEL_SPACING, NODE_SPACING);

    // if poset2 already present and compatible, build its layout too
    if ((int)g_poset2.cover_down.size() == n) {
      g_layout_p2 = build_layout_from_cover_down(g_poset2.cover_down, n, 0, LEVEL_SPACING, NODE_SPACING);
    } else {
      g_layout_p2 = PosetLayout{};
      g_layout_p2.centers = g_layout_p1.centers; // fallback so switching doesn’t crash
    }

    // apply active layout to meshes + centers
    apply_layout_to_meshes(active_layout());

    // edges (both)
    rebuild_both_edge_networks();

    std::cout << "[vis_poset] registered " << g_poset_meshes_2d.size()
              << " poset nodes (2D+3D meshes).\n";
}



static void debug_compare_levels(const PosetLayout& A, const PosetLayout& B) {
    const int n = (int)A.level.size();
    if ((int)B.level.size() != n) {
        std::cout << "[layout] level size mismatch: A=" << n << " B=" << B.level.size() << "\n";
        return;
    }

    int same = 0;
    int diff = 0;
    int max_abs_diff = 0;
    int first_diff = -1;

    for (int i = 0; i < n; ++i) {
        int da = A.level[i];
        int db = B.level[i];
        if (da == db) {
            same++;
        } else {
            diff++;
            if (first_diff < 0) first_diff = i;
            max_abs_diff = std::max(max_abs_diff, std::abs(da - db));
        }
    }

    std::cout << "[layout] levels: same=" << same << "/" << n
              << " diff=" << diff
              << " max_abs_diff=" << max_abs_diff;

    if (first_diff >= 0) {
        std::cout << " first_diff_node=" << first_diff
                  << " (P1=" << A.level[first_diff]
                  << " P2=" << B.level[first_diff] << ")";
    }
    std::cout << "\n";
}


void poset_ui() {
  ImGui::SetNextItemOpen(false, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("poset visualization", ImGuiTreeNodeFlags_DefaultOpen)) return;

  if (g_poset_meshes_2d.empty() && g_poset_meshes_3d.empty()) {
    ImGui::Text("poset: no meshes registered");
    return;
  }

  if (ImGui::Checkbox("show all poset 2D meshes", &g_show_poset_2d)) {
    apply_poset_visibility_from_masks();
  }

  if (ImGui::Checkbox("show all poset 3D meshes", &g_show_poset_3d)) {
    apply_poset_visibility_from_masks();
  }

  if (ImGui::Checkbox("show P1 edges (red)", &g_show_poset1_edges)) {
    if (g_poset1_edges) g_poset1_edges->setEnabled(g_show_poset1_edges);
  }

  if (ImGui::Checkbox("show P2 edges (green)", &g_show_poset2_edges)) {
    if (g_poset2_edges) g_poset2_edges->setEnabled(g_show_poset2_edges);
  }

  ImGui::Separator();
  ImGui::Text("Layout (node positions):");
  if (ImGui::RadioButton("Layout by P1", g_active_layout == ActiveLayout::P1)) {
    g_active_layout = ActiveLayout::P1;
    apply_layout_to_meshes(active_layout());
    rebuild_both_edge_networks();
  }
  ImGui::SameLine();
  if (ImGui::RadioButton("Layout by P2", g_active_layout == ActiveLayout::P2)) {
    g_active_layout = ActiveLayout::P2;
    apply_layout_to_meshes(active_layout());
    rebuild_both_edge_networks();
  }
  if (ImGui::Button("debug: compare P1 vs P2 levels")) {
    debug_compare_levels(g_layout_p1, g_layout_p2);
    }

}

const pst::Poset1& get_poset1() { return g_poset1; }
const pst2::Poset2& get_poset2() { return g_poset2; }

} // namespace viz_poset

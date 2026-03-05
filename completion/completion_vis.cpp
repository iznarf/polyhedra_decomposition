#include "completion_vis.h"

#include "poset_vis.h"

#include "input.h"
#include "poset.h"
#include "poset2.h"
#include "visualization.h"
#include "render_helpers.h"

#include "envelope.h"
#include "dedekind_poset.h"

#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>
#include <polyscope/surface_mesh.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace completion_vis {

namespace {

// -------------------- names --------------------
static std::string grid_name(int which) {
  if (which == 0) return "P1 completion nodes";
  if (which == 1) return "P2 completion nodes";
  return "P1 just flips completion nodes";
}

static std::string edge_name(int which) {
  if (which == 0) return "P1 completion edges";
  if (which == 1) return "P2 completion edges";
  return "P1 just flips completion edges";
}

static std::string tri2d_name(int which, int cutId) {
  const char* p = (which == 0) ? "P1 completion 2D "
              : (which == 1) ? "P2 completion 2D "
                             : "P1 just flips completion 2D ";
  return std::string(p) + std::to_string(cutId);
}

static std::string diag_name(int which, int cutId) {
  const char* p = (which == 0) ? "P1 completion diag "
              : (which == 1) ? "P2 completion diag "
                             : "P1 just flips completion diag ";
  return std::string(p) + std::to_string(cutId);
}

// -------------------- 2D overlay structure name lists --------------------
static std::vector<std::string> g_2d_structs_P1;
static std::vector<std::string> g_2d_structs_P2;
static std::vector<std::string> g_2d_structs_P1_flips;

// -------------------- grid positions per completion poset --------------------
static std::vector<glm::vec3> g_pos_P1;
static std::vector<glm::vec3> g_pos_P2;
static std::vector<glm::vec3> g_pos_P1_flips;

// -------------------- enabled flags per completion poset --------------------
static bool g_showGrid[3]  = {true, true, true};
static bool g_showEdges[3] = {true, true, true};
static bool g_show2D[3]    = {true, true, true};

static std::vector<std::string>& list_2d(int which) {
  if (which == 0) return g_2d_structs_P1;
  if (which == 1) return g_2d_structs_P2;
  return g_2d_structs_P1_flips;
}

static std::vector<glm::vec3>& pos_list(int which) {
  if (which == 0) return g_pos_P1;
  if (which == 1) return g_pos_P2;
  return g_pos_P1_flips;
}

// -------------------- helpers --------------------
static void remove_structure_if_exists(const std::string& nm) {
  if (polyscope::hasSurfaceMesh(nm))  polyscope::removeStructure(nm);
  if (polyscope::hasCurveNetwork(nm)) polyscope::removeStructure(nm);
}

static void clear_2d(int which) {
  auto& names = list_2d(which);
  for (const auto& nm : names) remove_structure_if_exists(nm);
  names.clear();
}

static void set_2d_enabled(int which, bool enabled) {
  auto& names = list_2d(which);
  for (const auto& nm : names) {
    if (polyscope::hasSurfaceMesh(nm))  polyscope::getSurfaceMesh(nm)->setEnabled(enabled);
    if (polyscope::hasCurveNetwork(nm)) polyscope::getCurveNetwork(nm)->setEnabled(enabled);
  }
}

// -------------------- P1 planar mesh helper --------------------
static bool build_planar_triangulation_from_P1_node_transformed(
    const df::InputData& D,
    const pst::Poset1& P1,
    int node_id,
    const glm::vec3& target_center,
    float cell_scale,
    std::vector<glm::vec3>& V2_out,
    std::vector<std::array<int,3>>& F_out)
{
  if (!(0 <= node_id && node_id < (int)P1.nodes.size())) return false;

  df::Tri2 tri = D.tri_poset;
  pst::replay_history_poset(tri, P1.nodes[node_id].history, D);

  auto ids      = viz::present_ids(tri);
  auto to_local = viz::make_local_index(ids);
  F_out         = viz_helpers::faces_from_triangles(tri, to_local);

  float S = 1.0f;
  V2_out = viz_helpers::make_planar_poset_vertices(ids, D.points2d, 0.0f, 2.0f, S);

  if (V2_out.empty() || F_out.empty()) return false;

  glm::vec2 c(0.f);
  for (auto& p : V2_out) c += glm::vec2(p.x, p.z);
  c /= (float)V2_out.size();

  for (auto& p : V2_out) {
    glm::vec2 q = glm::vec2(p.x, p.z) - c;
    q *= cell_scale;
    p.x = q.x + target_center.x;
    p.z = q.y + target_center.z;
    p.y = target_center.y;
  }
  return true;
}

// -------------------- P1 flips planar mesh helper --------------------
static bool build_planar_triangulation_from_P1_flips_node_transformed(
    const df::InputData& D,
    const pst::Poset_just_flips& P,
    int node_id,
    const glm::vec3& target_center,
    float cell_scale,
    std::vector<glm::vec3>& V2_out,
    std::vector<std::array<int,3>>& F_out)
{
  if (!(0 <= node_id && node_id < (int)P.nodes.size())) return false;

  df::Tri2 tri = D.tri_poset_just_flips; 
  pst::replay_history_poset(tri, P.nodes[node_id].history, D);

  auto ids      = viz::present_ids(tri);
  auto to_local = viz::make_local_index(ids);
  F_out         = viz_helpers::faces_from_triangles(tri, to_local);

  float S = 1.0f;
  V2_out = viz_helpers::make_planar_poset_vertices(ids, D.points2d, 0.0f, 2.0f, S);

  if (V2_out.empty() || F_out.empty()) return false;

  glm::vec2 c(0.f);
  for (auto& p : V2_out) c += glm::vec2(p.x, p.z);
  c /= (float)V2_out.size();

  for (auto& p : V2_out) {
    glm::vec2 q = glm::vec2(p.x, p.z) - c;
    q *= cell_scale;
    p.x = q.x + target_center.x;
    p.z = q.y + target_center.z;
    p.y = target_center.y;
  }
  return true;
}

// -------------------- P1 lifted triangles helper --------------------
static void append_lifted_triangles_from_P1_nodes(
    const df::InputData& D,
    const pst::Poset1& P1,
    const std::vector<int>& node_ids,
    std::vector<env_max::InputTriangle>& tris_out)
{
  int mesh_id_counter = 0;

  for (int node_id : node_ids) {
    if (!(0 <= node_id && node_id < (int)P1.nodes.size())) continue;

    df::Tri2 tri = D.tri_poset;
    pst::replay_history_poset(tri, P1.nodes[node_id].history, D);

    auto ids      = viz::present_ids(tri);
    auto to_local = viz::make_local_index(ids);
    auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

    float S = 1.0f;
    auto V3 = viz_helpers::make_lifted_poset_vertices(ids, D.points2d, 0.0f, 2.0f, S, 0.5f * S);

    const int mesh_id = ++mesh_id_counter;
    int face_counter = 0;

    for (const auto& f : faces) {
      const int i0 = f[0], i1 = f[1], i2 = f[2];
      if (i0 < 0 || i1 < 0 || i2 < 0) continue;
      if (i0 >= (int)V3.size() || i1 >= (int)V3.size() || i2 >= (int)V3.size()) continue;

      env_max::InputTriangle t;
      // world (x,y,z) -> envelope (x,z,y)
      t.p0 = { (double)V3[i0].x, (double)V3[i0].z, (double)V3[i0].y };
      t.p1 = { (double)V3[i1].x, (double)V3[i1].z, (double)V3[i1].y };
      t.p2 = { (double)V3[i2].x, (double)V3[i2].z, (double)V3[i2].y };

      t.tag.mesh_id = mesh_id;
      t.tag.face_id = face_counter++;
      tris_out.push_back(t);
    }
  }
}

// -------------------- P1 flips lifted triangles helper --------------------
static void append_lifted_triangles_from_P1_flips_nodes(
    const df::InputData& D,
    const pst::Poset_just_flips& P,
    const std::vector<int>& node_ids,
    std::vector<env_max::InputTriangle>& tris_out)
{
  int mesh_id_counter = 0;

  for (int node_id : node_ids) {
    if (!(0 <= node_id && node_id < (int)P.nodes.size())) continue;

    df::Tri2 tri = D.tri_poset_just_flips;
    pst::replay_history_poset(tri, P.nodes[node_id].history, D);

    auto ids      = viz::present_ids(tri);
    auto to_local = viz::make_local_index(ids);
    auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

    float S = 1.0f;
    auto V3 = viz_helpers::make_lifted_poset_vertices(ids, D.points2d, 0.0f, 2.0f, S, 0.5f * S);

    const int mesh_id = ++mesh_id_counter;
    int face_counter = 0;

    for (const auto& f : faces) {
      const int i0 = f[0], i1 = f[1], i2 = f[2];
      if (i0 < 0 || i1 < 0 || i2 < 0) continue;
      if (i0 >= (int)V3.size() || i1 >= (int)V3.size() || i2 >= (int)V3.size()) continue;

      env_max::InputTriangle t;
      t.p0 = { (double)V3[i0].x, (double)V3[i0].z, (double)V3[i0].y };
      t.p1 = { (double)V3[i1].x, (double)V3[i1].z, (double)V3[i1].y };
      t.p2 = { (double)V3[i2].x, (double)V3[i2].z, (double)V3[i2].y };

      t.tag.mesh_id = mesh_id;
      t.tag.face_id = face_counter++;
      tris_out.push_back(t);
    }
  }
}

// -------------------- grid + edges --------------------
static void build_grid_and_edges(int which,
                                 const dm_completion::DedekindPoset& C,
                                 const glm::vec3& center,
                                 float xSpacing,
                                 float zSpacing)
{
  const int M = (int)C.cuts.size();
  if (M <= 0) return;

  std::vector<glm::vec3> colors;
  colors.reserve((size_t)M);

  // old node colors per which
  const glm::vec3 oldColor =
      (which == 0) ? glm::vec3(1.f, 0.f, 0.f) // P1 red
    : (which == 1) ? glm::vec3(13.0f/255.0f, 100.0f/255.0f, 13.0f/255.0f) // P2 green
                   : glm::vec3(0.9f, 0.2f, 0.9f); // P1 flips purple (adjust if you want)

  const glm::vec3 newColor = glm::vec3(0.f, 0.f, 1.f); // blue

  for (int i = 0; i < M; ++i) {
    const bool isOld = ((int)C.cuts[i].minimal_elements_F.size() == 1);
    colors.push_back(isOld ? oldColor : newColor);
  }

  const std::string gname = grid_name(which);
  const std::string ename = edge_name(which);

  std::vector<glm::vec3> pos = pst_vis::register_poset_as_grid_colored(
      gname, M, C.levels, center, xSpacing, zSpacing, colors, 0.01f);

  pos_list(which) = pos;

  pst_vis::register_cover_up_as_edges(ename, pos, C.cover_up, glm::vec3(0.f, 0.f, 0.f));

  if (polyscope::hasPointCloud(gname))   polyscope::getPointCloud(gname)->setEnabled(g_showGrid[which]);
  if (polyscope::hasCurveNetwork(ename)) polyscope::getCurveNetwork(ename)->setEnabled(g_showEdges[which]);
}

// -------------------- 2D overlay --------------------
static void build_2d_overlay(int which,
                            const df::InputData& D,
                            const pst::Poset1& P1,
                            const pst::Poset_just_flips& P1_flips,
                            const dm_completion::DedekindPoset& C)
{
  clear_2d(which);

  auto& pos = pos_list(which);
  if ((int)pos.size() != (int)C.cuts.size()) return;

  auto& names = list_2d(which);

  const float cell_scale = 0.08f;
  const float y_plane    = 0.0f;
  const float radius     = 0.0013f;
  const glm::vec3 diagColor = glm::vec3(0.f, 0.f, 0.f); // black

  for (int cutId = 0; cutId < (int)C.cuts.size(); ++cutId) {

    // choose which cut feature you want to visualize
    //const auto& minF = C.cuts[cutId].minimal_elements_F;
    const auto& minF = C.cuts[cutId].maximal_elements_Iprime;

    const glm::vec3 center = pos[cutId];

    if ((int)minF.size() == 1) {

      const int node_id = minF[0];

      std::vector<glm::vec3> V2;
      std::vector<std::array<int,3>> F;

      bool ok = false;
      if (which == 2) {
        ok = build_planar_triangulation_from_P1_flips_node_transformed(
            D, P1_flips, node_id, center, cell_scale, V2, F);
      } else {
        ok = build_planar_triangulation_from_P1_node_transformed(
            D, P1, node_id, center, cell_scale, V2, F);
      }

      if (ok) {
        const std::string nm = tri2d_name(which, cutId);
        remove_structure_if_exists(nm);

        auto* m = polyscope::registerSurfaceMesh(nm, V2, F);
        m->setEdgeWidth(1.0f);
        m->setEnabled(g_show2D[which]);
        m->setSurfaceColor(glm::vec3(153.0f/255.0f, 204.0f/255.0f, 255.0f/255.0f));
        m->setEdgeColor(glm::vec3(0.f, 0.f, 0.f));

        names.push_back(nm);
      }

    } else if ((int)minF.size() > 1) {

      const std::string nm = diag_name(which, cutId);
      remove_structure_if_exists(nm);

      std::vector<env_max::InputTriangle> tris;
      tris.reserve(2000);

      if (which == 2) append_lifted_triangles_from_P1_flips_nodes(D, P1_flips, minF, tris);
      else            append_lifted_triangles_from_P1_nodes(D, P1, minF, tris);

      const bool ok = env_max::compute_and_register_lower_envelope_diagram_edges_2d(
          nm, tris, diagColor, radius, y_plane, center.x, center.z, cell_scale);

      if (ok) {
        if (polyscope::hasCurveNetwork(nm))
          polyscope::getCurveNetwork(nm)->setEnabled(g_show2D[which]);
        names.push_back(nm);
      } else {
        std::cout << "[completion_vis] diagram build failed for cut " << cutId << "\n";
      }
    }
  }
}

} // namespace (anonymous)

// ------------------------------------------------------------
// public API
// ------------------------------------------------------------
void build(int which,
           const df::InputData& D,
           const pst::Poset1& P1,
           const pst2::Poset2& P2,
           const pst::Poset_just_flips& P1_flips,
           const dm_completion::DedekindPoset& C,
           const glm::vec3& center,
           float xSpacing,
           float zSpacing)
{
  (void)P2;

  build_grid_and_edges(which, C, center, xSpacing, zSpacing);
  build_2d_overlay(which, D, P1, P1_flips, C);

  set_enabled(which, g_showGrid[which], g_showEdges[which], g_show2D[which]);
}

void set_enabled(int which, bool showGrid, bool showEdges, bool show2D) {
  g_showGrid[which]  = showGrid;
  g_showEdges[which] = showEdges;
  g_show2D[which]    = show2D;

  const std::string gname = grid_name(which);
  const std::string ename = edge_name(which);

  if (polyscope::hasPointCloud(gname))   polyscope::getPointCloud(gname)->setEnabled(showGrid);
  if (polyscope::hasCurveNetwork(ename)) polyscope::getCurveNetwork(ename)->setEnabled(showEdges);

  set_2d_enabled(which, show2D);
}

void clear(int which) {
  clear_2d(which);

  const std::string gname = grid_name(which);
  const std::string ename = edge_name(which);

  if (polyscope::hasPointCloud(gname))   polyscope::removeStructure(gname);
  if (polyscope::hasCurveNetwork(ename)) polyscope::removeStructure(ename);

  pos_list(which).clear();
}

void clear_all() {
  clear(0);
  clear(1);
  clear(2);
}

} // namespace completion_vis

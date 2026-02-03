
#include "envelope_vis.h"

#include "render_helpers.h"
#include "poset.h"
#include "input.h"
#include "visualization.h"

#include "envelope.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <glm/glm.hpp>

#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>

namespace envelope_vis {

namespace {

// ------------------------------------------------------------
// Internal state (so both UIs share one envelope instance)
// ------------------------------------------------------------
static polyscope::SurfaceMesh* g_env2 = nullptr; // planar (XZ)
static polyscope::SurfaceMesh* g_env3 = nullptr; // lifted (XYZ)

// toggles (optional; both UIs can drive them)
static bool g_show_env_2d = true;
static bool g_show_env_3d = true;

static void remove_mesh(polyscope::SurfaceMesh*& m) {
  if (m) {
    m->remove();
    m = nullptr;
  }
}

static void apply_env_toggles() {
  if (g_env2) g_env2->setEnabled(g_show_env_2d);
  if (g_env3) g_env3->setEnabled(g_show_env_3d);
}

// ------------------------------------------------------------
// Envelope input from a triangulation node in P1
// (replays history, builds lifted vertices, converts to InputTriangle)
// ------------------------------------------------------------
static void append_lifted_triangles_from_P1_history(
    const df::InputData& in,
    const pst::Poset1& P,
    int node_idx,
    int mesh_id,
    const Options& opt,
    std::vector<env_max::InputTriangle>& out)
{
  if (!(0 <= node_idx && node_idx < (int)P.nodes.size())) return;

  // replay triangulation from poset history
  df::Tri2 tri = in.tri_poset;
  pst::replay_history_poset(tri, P.nodes[node_idx].history, in);

  auto ids      = viz::present_ids(tri);
  auto to_local = viz::make_local_index(ids);
  auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

  // lift vertices to 3D (world coords)
  const float S = opt.scale;
  auto V3 = viz_helpers::make_lifted_poset_vertices(
      ids, in.points2d,
      opt.lift_min, opt.lift_max,
      S, 0.5f * S);

  // convert to envelope coords: world (x,y,z) -> envelope (x,z,y)
  int face_counter = 0;
  for (const auto& f : faces) {
    const int i0 = f[0], i1 = f[1], i2 = f[2];
    if (i0 < 0 || i1 < 0 || i2 < 0) continue;
    if (i0 >= (int)V3.size() || i1 >= (int)V3.size() || i2 >= (int)V3.size()) continue;

    const glm::vec3 P0 = V3[i0];
    const glm::vec3 P1 = V3[i1];
    const glm::vec3 P2 = V3[i2];

    env_max::InputTriangle t;
    t.p0 = { (double)P0.x, (double)P0.z, (double)P0.y };
    t.p1 = { (double)P1.x, (double)P1.z, (double)P1.y };
    t.p2 = { (double)P2.x, (double)P2.z, (double)P2.y };

    t.tag.mesh_id = mesh_id;
    t.tag.face_id = face_counter++;

    out.push_back(t);
  }
}

// ------------------------------------------------------------
// Register envelope meshes (2D + 3D)
// out.V is in envelope coords (x,z,y) so we convert back:
// lifted world (x,y,z) = (x, yheight, z)
// planar world (x,y,z) = (x, 0, z)
// ------------------------------------------------------------
static void register_envelope_meshes(
    const std::string& name2d,
    const std::string& name3d,
    const env_max::Mesh& out,
    const Options& opt)
{
  // clear old
  remove_mesh(g_env2);
  remove_mesh(g_env3);

  std::vector<glm::vec3> V_lifted;
  std::vector<glm::vec3> V_planar;
  V_lifted.reserve(out.V.size());
  V_planar.reserve(out.V.size());

  for (const auto& p : out.V) {
    const float x = (float)p[0];
    const float z = (float)p[1];
    const float y = (float)p[2];

    V_lifted.emplace_back(x, y, z);
    V_planar.emplace_back(x, 0.0f, z);
  }

  g_env2 = polyscope::registerSurfaceMesh(name2d, V_planar, out.F);
  g_env3 = polyscope::registerSurfaceMesh(name3d, V_lifted, out.F);

  g_env2->setEdgeWidth(opt.edge_width);
  g_env3->setEdgeWidth(opt.edge_width);

  // show which input mesh/node won each envelope triangle
  if (out.tri_tag.size() == out.F.size()) {
    std::vector<double> winner(out.F.size(), 0.0);
    for (size_t i = 0; i < out.tri_tag.size(); ++i) {
      winner[i] = (double)out.tri_tag[i][0]; // mesh_id
    }
    g_env2->addFaceScalarQuantity(opt.winner_quantity_name, winner);
    g_env3->addFaceScalarQuantity(opt.winner_quantity_name, winner);
  }

  // enable/disable
  g_show_env_2d = opt.show_2d;
  g_show_env_3d = opt.show_3d;
  apply_env_toggles();
}

} // namespace

// ============================================================
// Public API
// ============================================================

void clear() {
  remove_mesh(g_env2);
  remove_mesh(g_env3);
}

void set_enabled(bool show2d, bool show3d) {
  g_show_env_2d = show2d;
  g_show_env_3d = show3d;
  apply_env_toggles();
}

bool compute_and_register_from_P1_nodes(
    const std::string& name2d,
    const std::string& name3d,
    const df::InputData& in,
    const pst::Poset1& P,
    const std::vector<int>& node_ids,
    const Options& opt)
{
  // build envelope input
  std::vector<env_max::InputTriangle> tris;
  tris.reserve(opt.reserve_tris);

  int mesh_id = 0;
  for (int node_idx : node_ids) {
    append_lifted_triangles_from_P1_history(in, P, node_idx, mesh_id++, opt, tris);
  }

  if (tris.empty()) return false;

  // compute envelope
  env_max::Mesh out;
  const bool ok = env_max::compute_lower_envelope(tris, out);
  if (!ok) return false;

  // visualize
  register_envelope_meshes(name2d, name3d, out, opt);
  return true;
}



MeshPair build_mesh_pair_from_P1_node(
    const std::string& name2d,
    const std::string& name3d,
    const df::InputData& in,
    const pst::Poset1& P1,
    int node_id,
    float lift_min,
    float lift_max,
    float scale,
    float edge_width,
    bool add_global_id_quantity)
{
  MeshPair out;

  if (!(0 <= node_id && node_id < (int)P1.nodes.size())) {
    std::cout << "[build_mesh_pair_from_P1_node] invalid node id " << node_id << "\n";
    return out;
  }

  df::Tri2 tri = in.tri_poset;
  pst::replay_history_poset(tri, P1.nodes[node_id].history, in);

  auto ids      = viz::present_ids(tri);
  auto to_local = viz::make_local_index(ids);
  auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

  auto V2 = viz_helpers::make_planar_poset_vertices(ids, in.points2d, lift_min, lift_max, scale);
  auto V3 = viz_helpers::make_lifted_poset_vertices(ids, in.points2d, lift_min, lift_max, scale, 0.5f * scale);

  out.m2 = polyscope::registerSurfaceMesh(name2d, V2, faces);
  out.m3 = polyscope::registerSurfaceMesh(name3d, V3, faces);

  out.m2->setEdgeWidth(edge_width);
  out.m3->setEdgeWidth(edge_width);

  if (add_global_id_quantity) {
    viz_helpers::add_global_id_quantity(out.m2, ids);
    viz_helpers::add_global_id_quantity(out.m3, ids);
  }

  return out;
}






} // namespace envelope_vis


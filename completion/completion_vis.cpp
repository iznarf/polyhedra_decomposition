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

// Names
static std::string grid_name(int which)  { return which == 0 ? "P1 completion nodes" : "P2 completion nodes"; }
static std::string edge_name(int which)  { return which == 0 ? "P1 completion edges" : "P2 completion edges"; }
static std::string tri2d_name(int which, int cutId)  { return (which == 0 ? "P1 completion 2D " : "P2 completion 2D ") + std::to_string(cutId); }
static std::string diag_name(int which, int cutId)   { return (which == 0 ? "P1 completion diag " : "P2 completion diag ") + std::to_string(cutId); }

// Storage for toggle groups (2D meshes = surface meshes + curve networks)
static std::vector<std::string> g_2d_structs_P1;
static std::vector<std::string> g_2d_structs_P2;

// Store grid positions (needed to place per-cut overlays)
static std::vector<glm::vec3> g_pos_P1;
static std::vector<glm::vec3> g_pos_P2;

// Last known enabled states for toggles 
static bool g_showGrid[2]  = {true, true};
static bool g_showEdges[2] = {true, true};
static bool g_show2D[2]    = {true, true};

static std::vector<std::string>& list_2d(int which) {
  return (which == 0) ? g_2d_structs_P1 : g_2d_structs_P2;
}
static std::vector<glm::vec3>& pos_list(int which) {
  return (which == 0) ? g_pos_P1 : g_pos_P2;
}

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

// build a planar triangulation mesh for a P1 node id, then center+scale+translate into one completion grid cell
static bool build_planar_triangulation_from_P1_node_transformed(
    const df::InputData& D,
    const pst::Poset1& P1,
    int node_id,
    const glm::vec3& target_center,
    float cell_scale,
    std::vector<glm::vec3>& V2_out,
    std::vector<std::array<int,3>>& F_out){
    
    // sanity check: valid node id
    if (!(0 <= node_id && node_id < (int)P1.nodes.size())) return false;

    // replay history to get triangulation at this node
    df::Tri2 tri = D.tri_poset;
    pst::replay_history_poset(tri, P1.nodes[node_id].history, D);

    auto ids      = viz::present_ids(tri);
    auto to_local = viz::make_local_index(ids);
    F_out         = viz_helpers::faces_from_triangles(tri, to_local);

    // scaling of planar vertices 
    float S = 1.0f;
    V2_out = viz_helpers::make_planar_poset_vertices(ids, D.points2d, 0.0f, 2.0f, S);

    if (V2_out.empty() || F_out.empty()) return false;

    // center in x-z plane -> polyscope convetion (x,y,z) with y=0
    glm::vec2 c(0.f);
    for (auto& p : V2_out) c += glm::vec2(p.x, p.z);
    c /= (float)V2_out.size();

    for (auto& p : V2_out) {
        glm::vec2 q = glm::vec2(p.x, p.z) - c;
        q *= cell_scale;
        p.x = q.x + target_center.x;
        p.z = q.y + target_center.z;
        p.y = target_center.y; // usually 0
    }

    return true;
}

// collect lifted triangles (in envelope coords) from P1 nodes for building an envelope 
static void append_lifted_triangles_from_P1_nodes(
    const df::InputData& D,
    const pst::Poset1& P1,
    const std::vector<int>& node_ids,
    std::vector<env_max::InputTriangle>& tris_out){

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

// builds grid and edges for one completion poset 
static void build_grid_and_edges(int which,
                                 const dm_completion::DedekindPoset& C,
                                 const glm::vec3& center,
                                 float xSpacing,
                                 float zSpacing){
    const int M = (int)C.cuts.size();
    if (M <= 0) return;

    // determine old/new per number of elements in min(F) or max(I')
    // old node <=> |min(F)| == 1
    std::vector<glm::vec3> colors;
    colors.reserve((size_t)M);

    // color old nodes in red (P1) or green (P2), new nodes in blue
    const glm::vec3 oldColor = (which == 0)
        ? glm::vec3(1.f, 0.f, 0.f) // P1 red
        : glm::vec3(13.0f/255.0f, 100.0f/255.0f, 13.0f/255.0f); // P2 green

    const glm::vec3 newColor = glm::vec3(0.f, 0.f, 1.f); // blue

    for (int i = 0; i < M; ++i) {
        const bool isOld = ((int)C.cuts[i].minimal_elements_F.size() == 1);
        colors.push_back(isOld ? oldColor : newColor);
    }

    const std::string gname = grid_name(which);
    const std::string ename = edge_name(which);

    // compute grid nodes positions (returns grid positions)
    // this could be easier but we keeps this version for now 
    std::vector<glm::vec3> pos = pst_vis::register_poset_as_grid_colored(
        gname,
        M,
        C.levels,
        center,
        xSpacing,
        zSpacing,
        colors,
        0.01f);

    pos_list(which) = pos;

    // edges
    pst_vis::register_cover_up_as_edges(ename, pos, C.cover_up, glm::vec3(0.f, 0.f, 0.f));

    // apply enabled flags
    if (polyscope::hasPointCloud(gname)) polyscope::getPointCloud(gname)->setEnabled(g_showGrid[which]);
    if (polyscope::hasCurveNetwork(ename)) polyscope::getCurveNetwork(ename)->setEnabled(g_showEdges[which]);
}


// build 2D overlays for all cuts in the completion poset
static void build_2d_overlay(int which,
                            const df::InputData& D,
                            const pst::Poset1& P1,
                            const dm_completion::DedekindPoset& C){
    // clear old meshes first 
    clear_2d(which);

    // santiy check: which poset should be visualized 
    auto& pos = pos_list(which);
    if ((int)pos.size() != (int)C.cuts.size()) return;

    auto& names = list_2d(which);

    const float cell_scale = 0.08f;   // scaling for each mesh 
    const float y_plane    = 0.0f;  
    const float radius     = 0.0013f; // radius for envelope diagram edges
    const glm::vec3 diagColor = glm::vec3(0.f, 0.f, 1.f); // diagram color is black 

    for (int cutId = 0; cutId < (int)C.cuts.size(); ++cutId) {

        // get min(F) for this cut 
        //const auto& minF = C.cuts[cutId].minimal_elements_F;
        const auto& minF = C.cuts[cutId].maximal_elements_Iprime;
        const glm::vec3 center = pos[cutId];

        // start identifying what to draw: triangulation or diagram

        // old node: |min(F)| == 1 -> triangulation of that node
        if ((int)minF.size() == 1) {
        
            const int node_id = minF[0];

            std::vector<glm::vec3> V2;
            std::vector<std::array<int,3>> F;

            // build triangulation mesh for this node
            if (build_planar_triangulation_from_P1_node_transformed(D, P1, node_id, center, cell_scale, V2, F)) {
                const std::string nm = tri2d_name(which, cutId);
                remove_structure_if_exists(nm);

                auto* m = polyscope::registerSurfaceMesh(nm, V2, F);
                m->setEdgeWidth(1.0f);
                m->setEnabled(g_show2D[which]);
                // set color of 2D meshes to white 
                m->setSurfaceColor(glm::vec3(153.0f / 255.0f, 204.0f / 255.0f, 255.0f / 255.0f)); 
                m->setEdgeColor(glm::vec3(0.0f, 0.0f, 0.0f));

                names.push_back(nm);
            }
        // new node: |min(F)| > 1 -> envelope diagram of all min(F) nodes
        } else if ((int)minF.size() > 1) {
            // new: build diagram of the set minF
            const std::string nm = diag_name(which, cutId);
            remove_structure_if_exists(nm);

            // collect lifted triangles from all min(F) nodes
            std::vector<env_max::InputTriangle> tris;
            tris.reserve(2000);
            append_lifted_triangles_from_P1_nodes(D, P1, minF, tris);

            

            // build and register envelope diagram edges
            const bool ok = env_max::compute_and_register_lower_envelope_diagram_edges_2d(
                nm,
                tris,
                diagColor,
                radius,
                y_plane,
                center.x,
                center.z,
                cell_scale);

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


} // namespace

// ------------------------------------------------------------
// public API
// ------------------------------------------------------------
void build(int which,
           const df::InputData& D,
           const pst::Poset1& P1,
           const pst2::Poset2& P2,
           const dm_completion::DedekindPoset& C,
           const glm::vec3& center,
           float xSpacing,
           float zSpacing)
{
  (void)P2; // not needed for now, but passed for symmetry/future

  // keep existing enabled settings
  build_grid_and_edges(which, C, center, xSpacing, zSpacing);
  build_2d_overlay(which, D, P1, C);

  // enforce enabled states for all groups
  set_enabled(which, g_showGrid[which], g_showEdges[which], g_show2D[which]);
}

void set_enabled(int which, bool showGrid, bool showEdges, bool show2D) {
  g_showGrid[which]  = showGrid;
  g_showEdges[which] = showEdges;
  g_show2D[which]    = show2D;

  const std::string gname = grid_name(which);
  const std::string ename = edge_name(which);

  if (polyscope::hasPointCloud(gname))     polyscope::getPointCloud(gname)->setEnabled(showGrid);
  if (polyscope::hasCurveNetwork(ename))   polyscope::getCurveNetwork(ename)->setEnabled(showEdges);

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
}

} // namespace completion_vis

#include "dedekind_envelope_ui.h"

#include "envelope.h"

#include "render_helpers.h"
#include "poset.h"
#include "poset2.h"
#include "input.h"
#include "visualization.h"

#include "dedekind_cut.h"
#include "dedekind_poset.h"
#include "dedekind_completion.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <imgui.h>
#include <glm/glm.hpp>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

extern df::InputData g_in;

extern pst::Poset1  g_P1;
extern pst2::Poset2 g_P2;
extern bool g_has_P1;
extern bool g_has_P2;

extern dm_completion::DedekindPoset g_completion;

namespace dm_completion {

namespace {

// ------------------------------------------------------------
// UI state
// ------------------------------------------------------------
static char g_cut_buf[32] = "0";

static bool g_show_env_2d = true;
static bool g_show_env_3d = true;

static bool g_show_cut_nodes_2d = true;
static bool g_show_cut_nodes_3d = true;

// envelope meshes
static polyscope::SurfaceMesh* g_env2 = nullptr;
static polyscope::SurfaceMesh* g_env3 = nullptr;

// cut node meshes
static std::vector<polyscope::SurfaceMesh*> g_cut_nodes_2d;
static std::vector<polyscope::SurfaceMesh*> g_cut_nodes_3d;

// ------------------------------------------------------------
// helpers
// ------------------------------------------------------------
static void remove_mesh(polyscope::SurfaceMesh*& m) {
    if (m) { m->remove(); m = nullptr; }
}

static void clear_env() {
    remove_mesh(g_env2);
    remove_mesh(g_env3);
}

static void clear_cut_node_meshes() {
    for (auto* m : g_cut_nodes_2d) if (m) m->remove();
    for (auto* m : g_cut_nodes_3d) if (m) m->remove();
    g_cut_nodes_2d.clear();
    g_cut_nodes_3d.clear();
}

static void apply_env_toggles() {
    if (g_env2) g_env2->setEnabled(g_show_env_2d);
    if (g_env3) g_env3->setEnabled(g_show_env_3d);
}

// ------------------------------------------------------------
// Envelope input from triangulation
// ------------------------------------------------------------
static void append_lifted_triangles_from_P1_history(int node_idx,int mesh_id,std::vector<env_max::InputTriangle>& out){
    if (!(0 <= node_idx && node_idx < (int)g_P1.nodes.size())) return;

    // fix: we already built the triangulation before, we do not have to build it here again
    // give function tri as input! 

    df::Tri2 tri = g_in.tri_poset; // start from poset base triangulation for replay. this is upper triangulation 
    pst::replay_history_poset(tri, g_P1.nodes[node_idx].history, g_in);

    // build mesh data
    auto ids      = viz::present_ids(tri); // vertex ids used by triangultion
    auto to_local = viz::make_local_index(ids); // mapping from vertex id to local index
    auto faces    = viz_helpers::faces_from_triangles(tri, to_local); // triangle faces as index triples


    // lift vertices to 3D 
    float S = 20.0f;    // scaling factor, no math meaning here
    auto V3 = viz_helpers::make_lifted_poset_vertices(ids, g_in.points2d, 0.0f, 20.0f, S, 0.5f * S);

    // now convert to CGAL input surface patch and store where it came from
    int face_counter = 0;
    // triangle faces 
    for (const auto& f : faces) {
        const int i0 = f[0], i1 = f[1], i2 = f[2];
        if (i0 < 0 || i1 < 0 || i2 < 0) continue; // invalid index check 
        if (i0 >= (int)V3.size() || i1 >= (int)V3.size() || i2 >= (int)V3.size()) continue; // invalid index check

        // get triangle vertex positions in 3D
        const glm::vec3 P0 = V3[i0];
        const glm::vec3 P1 = V3[i1];
        const glm::vec3 P2 = V3[i2];

        // build envelope input triangle
        env_max::InputTriangle t;

        // important: convert to envelope coordinates: world (x,y,z) -> envelope (x,z,y)
        // this is CGAL / polyscope coordinate system difference
        t.p0 = { (double)P0.x, (double)P0.z, (double)P0.y };
        t.p1 = { (double)P1.x, (double)P1.z, (double)P1.y };
        t.p2 = { (double)P2.x, (double)P2.z, (double)P2.y };

        // tag the triangle with its origin triangulation node
        t.tag.mesh_id = mesh_id; // which triangulation this belongs to 
        t.tag.face_id = face_counter++; // to which face in the triangulation it belongs 

        // global envelope input triangle list
        out.push_back(t);
    }
}

// ------------------------------------------------------------
// Register envelope meshes
// ------------------------------------------------------------
static void register_envelope_meshes(const std::string& name2d,const std::string& name3d,const env_max::Mesh& out){
    clear_env();

    std::vector<glm::vec3> V_lifted;
    std::vector<glm::vec3> V_planar;

    for (const auto& p : out.V) {
        const float x = (float)p[0];
        const float z = (float)p[1];
        const float y = (float)p[2];

        V_lifted.emplace_back(x, y, z);
        V_planar.emplace_back(x, 0.0f, z);
    }

    g_env2 = polyscope::registerSurfaceMesh(name2d, V_planar, out.F);
    g_env3 = polyscope::registerSurfaceMesh(name3d, V_lifted, out.F);

    g_env2->setEdgeWidth(1.0f);
    g_env3->setEdgeWidth(1.0f);

    if (out.tri_tag.size() == out.F.size()) {
        std::vector<double> winner(out.F.size());
        for (size_t i = 0; i < out.tri_tag.size(); ++i)
        winner[i] = (double)out.tri_tag[i][0];

        g_env2->addFaceScalarQuantity("winner node", winner);
        g_env3->addFaceScalarQuantity("winner node", winner);
    }

    apply_env_toggles();
}

// ------------------------------------------------------------
// Show all max(I') node meshes
// ------------------------------------------------------------
// build cut node meshes from P1 step history 
static void build_cut_node_meshes_from_P1_history(const std::vector<int>& nodes,const std::string& prefix){
    clear_cut_node_meshes();

    for (int node_id : nodes) {
        // node id not valid 
        if (!(0 <= node_id && node_id < (int)g_P1.nodes.size())) continue;

        // build triangulation from history
        df::Tri2 tri = g_in.tri_poset;
        pst::replay_history_poset(tri, g_P1.nodes[node_id].history, g_in);

        auto ids      = viz::present_ids(tri);
        auto to_local = viz::make_local_index(ids);
        auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

        float S = 20.0f;
        auto V2 = viz_helpers::make_planar_poset_vertices(ids, g_in.points2d, 0.0f, 20.0f, S);

        auto V3 = viz_helpers::make_lifted_poset_vertices(ids, g_in.points2d, 0.0f, 20.0f, S, 0.5f * S);

        // register meshes
        auto* m2 = polyscope::registerSurfaceMesh(prefix + " node " + std::to_string(node_id) + " 2D", V2, faces);
        auto* m3 = polyscope::registerSurfaceMesh(prefix + " node " + std::to_string(node_id) + " 3D", V3, faces);

        m2->setEdgeWidth(1.0f);
        m3->setEdgeWidth(1.0f);

        g_cut_nodes_2d.push_back(m2);
        g_cut_nodes_3d.push_back(m3);
    }

  // apply current group visibility (not necessary)
  for (auto* m : g_cut_nodes_2d) if (m) m->setEnabled(g_show_cut_nodes_2d);
  for (auto* m : g_cut_nodes_3d) if (m) m->setEnabled(g_show_cut_nodes_3d);
}

} // namespace

// ============================================================
// UI
// ============================================================
void dedekind_envelope_ui() {

  ImGui::SetNextItemOpen(false, ImGuiCond_Once);
  if (!ImGui::CollapsingHeader("Envelope", ImGuiTreeNodeFlags_DefaultOpen))
    return;

  const int C = (int)g_completion.cuts.size();

  ImGui::InputText("cut id", g_cut_buf, IM_ARRAYSIZE(g_cut_buf));
  const int cut_id = std::atoi(g_cut_buf);

  // ---- group toggles ----
  bool env2dChanged = ImGui::Checkbox("env 2D", &g_show_env_2d); ImGui::SameLine();
  bool env3dChanged = ImGui::Checkbox("env 3D", &g_show_env_3d);

  if (env2dChanged && g_env2) g_env2->setEnabled(g_show_env_2d);
  if (env3dChanged && g_env3) g_env3->setEnabled(g_show_env_3d);

  bool cut2dChanged = ImGui::Checkbox("cut nodes 2D", &g_show_cut_nodes_2d); ImGui::SameLine();
  bool cut3dChanged = ImGui::Checkbox("cut nodes 3D", &g_show_cut_nodes_3d);

  if (cut2dChanged)
    for (auto* m : g_cut_nodes_2d) if (m) m->setEnabled(g_show_cut_nodes_2d);

  if (cut3dChanged)
    for (auto* m : g_cut_nodes_3d) if (m) m->setEnabled(g_show_cut_nodes_3d);

  // ---- info ----
  if (!(0 <= cut_id && cut_id < C)) {
    ImGui::Text("Valid cut indices: 0..%d", std::max(0, C - 1));
  } else {
    ImGui::Text("max(I') size: %d",
      (int)g_completion.cuts[cut_id].maximal_elements_Iprime.size());
  }

  // ---- compute ----
  if (ImGui::Button("Compute envelope of max(I')")) {

    if (0 <= cut_id && cut_id < C) {
        // get node ids in max (I')    
        const auto& nodes = g_completion.cuts[cut_id].minimal_elements_F;

        // build 2D and 3D meshes from P1 step history 
        build_cut_node_meshes_from_P1_history(nodes, "Cut " + std::to_string(cut_id));

        // collect lifted triangles of lifted triangulations in max(I')
        std::vector<env_max::InputTriangle> tris;
        tris.reserve(4000);

        // get lifted triangles from each node in max(I')
        // collects all lifted triangle faces, them into surface patches for envelope computation
        for (int node_id : nodes)
            append_lifted_triangles_from_P1_history(node_id, node_id, tris);

        // compute upper envelope
        // this is reconstructed as triangle mesh with per triangle tags
        env_max::Mesh out;
        if (env_max::compute_lower_envelope(tris, out)) {
            register_envelope_meshes("envelope 2D"," envelope 3D",out);
        }

    }
  }

  ImGui::SameLine();

  if (ImGui::Button("Clear")) {
    clear_env();
    clear_cut_node_meshes();
  }
}

} 

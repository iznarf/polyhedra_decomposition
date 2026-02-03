#pragma once


#include <string>
#include <vector>
#include <polyscope/surface_mesh.h>
#include <glm/glm.hpp>

namespace df { struct InputData; }
namespace pst { struct Poset1; }

namespace envelope_vis {


struct Options {

  // ---- geometry / lifting ----
  // used to build lifted vertices via viz_helpers::make_lifted_poset_vertices(...)
  float scale    = 1.0f;
  float lift_min = 0.0f;
  float lift_max = 2.0f;

  // ---- visualization ----
  bool show_2d = true;
  bool show_3d = true;

  float edge_width = 1.0f;

  // ---- envelope computation ----
  // just a reserve hint 
  int reserve_tris = 0;

  // ---- quantities ----
  std::string winner_quantity_name = "winner mesh";
};

// Removes currently registered envelope meshes (the internal cached pointers)
void clear();

// Enable/disable the currently registered envelope meshes
void set_enabled(bool show2d, bool show3d);

// Build envelope from the list of P1 node ids and register meshes.
// Returns true on success.
bool compute_and_register_from_P1_nodes(
    const std::string& name2d,
    const std::string& name3d,
    const df::InputData& in,
    const pst::Poset1& P,
    const std::vector<int>& node_ids,
    const Options& opt = Options{}
);

struct MeshPair {
  polyscope::SurfaceMesh* m2 = nullptr;
  polyscope::SurfaceMesh* m3 = nullptr;
};

// Build 2D+3D surface mesh for a poset node (triangulation), using P1 history
MeshPair build_mesh_pair_from_P1_node(
    const std::string& name2d,
    const std::string& name3d,
    const df::InputData& in,
    const pst::Poset1& P1,
    int node_id,
    float lift_min,
    float lift_max,
    float scale,
    float edge_width = 1.0f,
    bool add_global_id_quantity = true
);



} // namespace envelope_vis


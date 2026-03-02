#pragma once

#include "poset_ear_star.h"
#include "ear_star_vis.h"   
#include <string>
#include <vector>

namespace pst_es_viz {

// Registers every triangulation in P as its own surface mesh in Polyscope,
// laid out on a grid by translating the vertex positions.
// - tri_start: the root triangulation used for reconstruction
// - points2d:  array indexed by global vertex id (must include star_id too)
// - star_id, star_point: for reconstruction
// - base_name: mesh names will be base_name + "_" + index
// - spacing: grid spacing in world units 
void register_poset_triangulations_grid(
    const df::Tri2& tri_start,
    const pst_es::FlipPoset& P,
    const std::vector<df::P2>& points2d,
    df::vertex_id star_id,
    const df::P2& star_point,
    const std::string& base_name = "T",
    double spacing = 4.0
);

} // namespace pst_es_viz
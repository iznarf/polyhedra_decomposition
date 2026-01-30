#pragma once

#include "input.h"          
#include <glm/glm.hpp>     

#include <vector>
#include <array>
#include <unordered_map>

namespace polyscope { class SurfaceMesh; }

namespace viz_helpers {

// faces as local index triples for planar triangulation
std::vector<std::array<int,3>>
faces_from_triangles(const df::Tri2& t,
                     const std::unordered_map<df::vertex_id,int>& to_local);

// attach global id scalar quantity to polyscope mesh
void add_global_id_quantity(polyscope::SurfaceMesh* mesh,
                            const std::vector<df::vertex_id>& ids);

// make planar vertex positions for a given triangulation centered at (cx,cz)
// in the X-Z plane, with a small uniform scale
std::vector<glm::vec3>
make_planar_poset_vertices(const std::vector<df::vertex_id>& ids,
                           const std::vector<df::P2>& points2d,
                           float cx, float cz, float scale);

// make lifted vertices above same grid center (cx,cz)
std::vector<glm::vec3>
make_lifted_poset_vertices(const std::vector<df::vertex_id>& ids,
                           const std::vector<df::P2>& points2d,
                           float cx, float cz,
                           float scale_xy, float scale_z);

glm::vec3 bbox_center(const std::vector<glm::vec3>& V);

void recenter(std::vector<glm::vec3>& V, const glm::vec3& c);

// attach node index scalar quantity to polyscope mesh
void add_node_id_quantity(polyscope::SurfaceMesh* mesh, int nodeIndex);

} // namespace viz_helpers



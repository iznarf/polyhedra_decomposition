#include "visualization.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <glm/glm.hpp>

#include <CGAL/Triangulation_2.h>
#include <unordered_map>

namespace viz {

using SurfaceMesh = polyscope::SurfaceMesh;


namespace {


// takes CGAL 2D triangulation and optional height (omega) function
// it produces a vertex list V and face list F for polyscope
template <typename Tri2, typename OmegaOrNull>
void extract_mesh(const Tri2& tri,  
                  const OmegaOrNull* omega, // pointer to omega function or nullptr
                  std::vector<glm::vec3>& V, // vertex positions
                  std::vector<std::array<size_t,3>>& F) // face indices, one face consists of 3 vertex indices
{
  // prepare look up table vidx (vertex_handle, index in vertex list)
  // maps CGAL vertex handles (pointers inside that CGAL object) to index in V (vertex list)
  std::unordered_map<typename Tri2::Vertex_handle, size_t> vidx;
  // allocate space for vertex indices
  vidx.reserve(static_cast<size_t>(std::distance(tri.finite_vertices_begin(),
                                                 tri.finite_vertices_end())));
  V.clear(); F.clear();

  // iterate over all finite vertices to extract positions
  for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit) {
    const auto& p = vit->point();
    float x = static_cast<float>(p.x());
    float z = static_cast<float>(p.y());
    float y = 0.0f;
    // lift z if omega is provided
    if constexpr (!std::is_same<OmegaOrNull, void>::value) {
      if (omega) y = static_cast<float>((*omega)(p));
    } 
    size_t idx = V.size();
    // store position
    V.emplace_back(x, y, z);
    // store index in lookup table
    vidx.emplace(vit, idx);
  }

  // iterate over all finite faces to extract connectivity
  // each face is a triangle 
  // each triangle face has 3 vertex handles
  for (auto fit = tri.finite_faces_begin(); fit != tri.finite_faces_end(); ++fit) {
    size_t i0 = vidx.at(fit->vertex(0));
    size_t i1 = vidx.at(fit->vertex(1));
    size_t i2 = vidx.at(fit->vertex(2));
    F.push_back({i0, i1, i2});
  }
}

// get or register a SurfaceMesh in Polyscope with given name, vertex positions V and faces F
// if the mesh already exists, return the existing one
SurfaceMesh* get_or_register(const std::string& name,
                             const std::vector<glm::vec3>& V,
                             const std::vector<std::array<size_t,3>>& F) {
  SurfaceMesh* M = nullptr;
  if (polyscope::hasSurfaceMesh(name)) {
    M = polyscope::getSurfaceMesh(name);
  } else {
    M = polyscope::registerSurfaceMesh(name, V, F);
    if (M) M->setSmoothShade(false);
  }
  return M;
}


} // anon

// initialize polyscope
void init() {
  polyscope::init();
}

// this is for visualizing the planar triangulation, therefore omega is nullptr
void register_planar_triangulation(const df::Tri2& tri, const std::string& name) {
  std::vector<glm::vec3> V;
  std::vector<std::array<size_t,3>> F;
  extract_mesh<df::Tri2, void>(tri, nullptr, V, F);
  auto* M = get_or_register(name, V, F);
  (void)M;
}

// this is for visualizing the lifted triangulation, therefore omega is provided
void register_lifted_triangulation(const df::Tri2& tri,
                                   const df::OmegaQuad& omega,
                                   const std::string& name) {
  std::vector<glm::vec3> V;
  std::vector<std::array<size_t,3>> F;
  extract_mesh(tri, &omega, V, F);
  auto* M = get_or_register(name, V, F);
  (void)M;
}


// these functions are then for updating existing meshes in Polyscope after we performed a flip i think
void update_planar_triangulation(const df::Tri2& tri, const std::string& name) {
  std::vector<glm::vec3> V;
  std::vector<std::array<size_t,3>> F;
  extract_mesh<df::Tri2, void>(tri, nullptr, V, F);

  if (polyscope::hasSurfaceMesh(name)) {
    polyscope::removeSurfaceMesh(name, /*errorIfAbsent=*/false);
  }
  auto* M = polyscope::registerSurfaceMesh(name, V, F);
  if (M) M->setSmoothShade(false);
}

void update_lifted_triangulation(const df::Tri2& tri,
                                 const df::OmegaQuad& omega,
                                 const std::string& name) {
  std::vector<glm::vec3> V;
  std::vector<std::array<size_t,3>> F;
  extract_mesh(tri, &omega, V, F);

  if (polyscope::hasSurfaceMesh(name)) {
    polyscope::removeSurfaceMesh(name, /*errorIfAbsent=*/false);
  }
  auto* M = polyscope::registerSurfaceMesh(name, V, F);
  if (M) M->setSmoothShade(false);
}


void show() { polyscope::show(); }

} // namespace viz

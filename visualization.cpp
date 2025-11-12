#include "visualization.h"
#include "input.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <array>


using df::vertex_id;
using glm::vec3;

namespace  {


// vertcies for planar triangulation z = 0
std::vector<vec3> points_planar(const std::vector<vertex_id>& ids, const std::vector<df::P2>& P) {
    std::vector<vec3> V; V.reserve(ids.size());
    for (auto id : ids) {
        const auto& p = P[id];
        V.push_back(vec3((float)p.x(), 0.0f, (float)p.y()));
    }
    return V;
}

// vertices for lifted triangulation z = x^2 + y^2 
std::vector<vec3> points_lifted(const std::vector<vertex_id>& ids, const std::vector<df::P2>& P) {
    std::vector<vec3> V; V.reserve(ids.size());
    for (auto id : ids) {
        const auto& p = P[id];
        float z = (float)(p.x()*p.x() + p.y()*p.y());
        V.push_back(vec3((float)p.x(), z, (float)p.y()));
    }
    return V;
}

// faces as local index triples
std::vector<std::array<int,3>> faces_from_triangles(const df::Tri2& t, const std::unordered_map<vertex_id,int>& to_local) {
    std::vector<std::array<int,3>> F;
    F.reserve(t.number_of_faces());
    for (auto f = t.finite_faces_begin(); f != t.finite_faces_end(); ++f) {
        int a = to_local.at(f->vertex(0)->info());
        int b = to_local.at(f->vertex(1)->info());
        int c = to_local.at(f->vertex(2)->info());
        F.push_back({a,b,c});
    }
    return F;
}

// register one triangulation as (2D mesh + lifted mesh)
void register_triangulation_as_mesh(const df::Tri2& tri,
                         const std::vector<df::P2>& points2d,
                         const std::string& name_planar,
                         const std::string& name_lifted) {
    auto ids   = viz::present_ids(tri); // which vertices does the triangulation contain -> global indices
    auto to_local = viz::make_local_index(ids); // map to local indices for polyscope

    auto vertices_2d = points_planar(ids, points2d); // coordinates of points in 2D
    auto faces  = faces_from_triangles(tri, to_local); // faces as local index triples
    polyscope::registerSurfaceMesh(name_planar, vertices_2d, faces); // planar mesh 

    auto vertices_3d = points_lifted(ids, points2d); // coordinates of lifted points 
    auto mesh = polyscope::registerSurfaceMesh(name_lifted, vertices_3d, faces); // lifted mesh
    mesh->setTransparency(0.5); 
}

} // anon

namespace viz {

static inline glm::vec3 lift_paraboloid(const df::P2& p) {
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  return glm::vec3((float)x, (float)(x*x + y*y), (float)y); // y-up
}

void show_flip_tetra(const df::InputData& D, df::vertex_id ia, df::vertex_id ib,  const std::string& label) {

    // map global id -> handle
    auto& tri = D.tri_current;
    auto& index_to_vertex_handle = D.index_to_vertex_handle_current;

    // use the global id -> vertex handle map to get vertex handles
    df::Tri2::Vertex_handle va = index_to_vertex_handle[ia];
    df::Tri2::Vertex_handle vb = index_to_vertex_handle[ib];

    // find opposite vertices c,d for edge (va,vb)
    df::Tri2::Face_handle f; int i = -1;
    if (!tri.is_edge(va, vb, f, i)) return;
    if (tri.is_infinite(f)) return;

    auto g = f->neighbor(i);
    int mi = tri.mirror_index(f, i);

    auto vc = f->vertex(i);
    auto vd = g->vertex(mi);


    const df::P2& A2 = va->point();
    const df::P2& B2 = vb->point();
    const df::P2& C2 = vc->point();
    const df::P2& D2 = vd->point();

    std::vector<glm::vec3> V = {
        lift_paraboloid(A2), lift_paraboloid(B2),
        lift_paraboloid(C2), lift_paraboloid(D2)
    };
    std::vector<std::array<size_t,3>> F = { {0,1,2}, {0,3,1}, {1,3,2}, {0,2,3} };

    auto mesh = polyscope::registerSurfaceMesh(label, V, F);

    mesh->setTransparency(0.5); 


}


// collect the global indices present in a triangulation
std::vector<vertex_id> present_ids(const df::Tri2& t) {
    std::vector<vertex_id> ids;
    ids.reserve(t.number_of_vertices());
    for (auto v = t.finite_vertices_begin(); v != t.finite_vertices_end(); ++v)
        ids.push_back(v->info());
    return ids;
}

// map global index to local compact index [0,..,V-1] for polyscope
std::unordered_map<vertex_id,int> make_local_index(const std::vector<vertex_id>& ids) {
  std::unordered_map<vertex_id,int> m;
  m.reserve(ids.size());
  for (int i = 0; i < (int)ids.size(); ++i) m[ids[i]] = i;
  return m;
}

void show_four_meshes(const df::InputData& D) {

  // lower (all points) as planar & lifted
  register_triangulation_as_mesh(D.tri_lower, D.points2d, "lower 2D", "lower lifted");

  // upper (hull) as planar & lifted
  register_triangulation_as_mesh(D.tri_upper, D.points2d, "upper 2D", "upper lifted");

}

} // namespace viz



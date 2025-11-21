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
    //mesh->setTransparency(0.5); 
}

} // anon

namespace viz {

static inline glm::vec3 lift_paraboloid(const df::P2& p) {
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  return glm::vec3((float)x, (float)(x*x + y*y), (float)y); // y-up
}

void show_flip_tetra(const df::InputData& D, df::vertex_id ia, df::vertex_id ib,  const std::string& label) {

    const auto& tri = D.tri_current;

    // find the edge (ia, ib) in the current triangulation by vertex ids
    df::Tri2::Face_handle f;
    int i = -1;
    bool found = false;

    for (auto e = tri.finite_edges_begin(); e != tri.finite_edges_end(); ++e) {
        auto fe = e->first;
        int  ei = e->second;

        auto va = fe->vertex(tri.cw(ei));
        auto vb = fe->vertex(tri.ccw(ei));

        df::vertex_id ja = va->info();
        df::vertex_id jb = vb->info();

        if ((ja == ia && jb == ib) || (ja == ib && jb == ia)) {
            f     = fe;
            i     = ei;
            found = true;
            break;
        }
    }

        if (!found) {
            std::cout << "[viz] show_flip_tetra: edge (" << ia << "," << ib
                    << ") not found in tri_current\n";
            return;
        }

        // if the incident face or its neighbor is infinite, just bail out
        if (tri.is_infinite(f)) return;
        auto g = f->neighbor(i);
        if (tri.is_infinite(g)) return;

        // endpoints of the edge and the two opposite vertices
        auto va = f->vertex(tri.cw(i));
        auto vb = f->vertex(tri.ccw(i));
        auto vc = f->vertex(i);
        int  mi = tri.mirror_index(f, i);
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
        mesh->setTransparency(0.5f);
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


void show_or_update_current(const df::InputData& D) {
  // build fresh buffers from the current triangulation
  const auto& tri = D.tri_current;

  const auto ids      = present_ids(tri);
  const auto to_local = make_local_index(ids);
  const auto V2       = points_planar(ids,  D.points2d);
  const auto V3       = points_lifted(ids,  D.points2d);
  const auto F        = faces_from_triangles(tri, to_local);

  // remove old meshes only if they exist
  if (polyscope::hasSurfaceMesh("current 2D")) {
    polyscope::getSurfaceMesh("current 2D")->remove();
  }
  if (polyscope::hasSurfaceMesh("current lifted")) {
    polyscope::getSurfaceMesh("current lifted")->remove();
  }

  // (re)register meshes
  polyscope::registerSurfaceMesh("current 2D", V2, F);
  auto* m3 = polyscope::registerSurfaceMesh("current lifted", V3, F);
  m3->setTransparency(0.5f);
}


void show_or_update_replay(const df::InputData& D) {
    const auto& tri = D.tri_replay;

    const auto ids      = present_ids(tri);
    const auto to_local = make_local_index(ids);
    const auto V2       = points_planar(ids,  D.points2d);
    const auto V3       = points_lifted(ids,  D.points2d);
    const auto F        = faces_from_triangles(tri, to_local);

    // remove old replay meshes if they exist
    if (polyscope::hasSurfaceMesh("replay 2D")) {
        polyscope::getSurfaceMesh("replay 2D")->remove();
    }
    if (polyscope::hasSurfaceMesh("replay lifted")) {
        polyscope::getSurfaceMesh("replay lifted")->remove();
    }

    polyscope::registerSurfaceMesh("replay 2D", V2, F);
    auto* m3 = polyscope::registerSurfaceMesh("replay lifted", V3, F);
    m3->setTransparency(0.5f);
}



} // namespace viz



#include "regularity_check.h"

#include <CGAL/enum.h>
#include <unordered_map>
#include <iomanip>
#include <iostream>
#include <vector>

namespace df {
namespace reg {

using K   = df::K;
using P2  = df::P2;
using P3  = df::P3;
using Tri = df::Tri2;

static inline P3 lift_point(const P2& p, const OmegaQuad& omega) {
  return P3(p.x(), p.y(), omega(p));
}

// ensure (a,b,c) CCW in 2D, then evaluate orientation(a',b',c',d')
// returns POSITIVE if d' above plane a'b'c', NEGATIVE if below, COLLINEAR if coplanar
static inline CGAL::Orientation oriented_height_sign(const P2& a2, const P2& b2, const P2& c2,
                                                     const P3& a3, const P3& b3, const P3& c3, const P3& d3) {
  auto o2 = CGAL::orientation(a2, b2, c2);
  if (o2 == CGAL::RIGHT_TURN) return CGAL::orientation(a3, c3, b3, d3);
  if (o2 == CGAL::COLLINEAR)  return CGAL::COLLINEAR;
  return CGAL::orientation(a3, b3, c3, d3);
}

// Quad strictly convex around edge ab: c and d on opposite sides of ab, neither collinear.
static inline bool quad_strictly_convex(const P2& a, const P2& b, const P2& c, const P2& d) {
  auto oc = CGAL::orientation(a, b, c);
  auto od = CGAL::orientation(a, b, d);
  return (oc != CGAL::COLLINEAR) && (od != CGAL::COLLINEAR) && (oc != od);
}

// build CGAL vertex handle -> polyscope vertex index map
static std::unordered_map<Tri::Vertex_handle, size_t>
make_vertex_index_map(const Tri& tri) {
  std::unordered_map<Tri::Vertex_handle, size_t> vidx;
  vidx.reserve(static_cast<size_t>(std::distance(tri.finite_vertices_begin(),
                                                 tri.finite_vertices_end())));
  size_t idx = 0;
  for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit) {
    vidx.emplace(vit, idx++);
  }
  return vidx;
}


std::vector<EdgeRecord> collect_flippable_edges_with_ids(const Tri& tri, const OmegaQuad& omega) {
    std::vector<EdgeRecord> out; 
    out.reserve(64);

    // vertex handle map to polyscope vertex indices
    auto vidx = make_vertex_index_map(tri);
    // iterate over all finite edges 
    for (auto eit = tri.finite_edges_begin(); eit != tri.finite_edges_end(); ++eit) {
        auto f = eit->first; // incident face of the edge
        int  i = eit->second; // index of opposite vertex of edge eit in face f 

    // internal edge?
    auto g = f->neighbor(i); // opposite face across edge
    if (tri.is_infinite(g)) continue; // if g is infinite, it is a boundary edge. we skip it since we can just flip internal edges 

    // extract the quad around the edge  
    auto va = f->vertex(tri.cw(i)); // vertex handle of endpoint of edge eit
    auto vb = f->vertex(tri.ccw(i)); // vertex handle of endpoint of edge eit
    auto vc = f->vertex(i); // vertex handle of opposite vertex of edge eit in face f 
    int  j  = tri.mirror_index(f, i); // index of opposite vertex of edge eit in face g
    auto vd = g->vertex(j); // vertex handle of vertex opposite to edge eit in face g

    // extract 2D points of the quad around the edge 
    const P2 &a2 = va->point(), &b2 = vb->point(), &c2 = vc->point(), &d2 = vd->point();

    // check if the quad is convex
    if (!quad_strictly_convex(a2, b2, c2, d2)) continue;

    // lift points to 3D
    P3 a3 = lift_point(a2, omega);
    P3 b3 = lift_point(b2, omega);
    P3 c3 = lift_point(c2, omega);
    P3 d3 = lift_point(d2, omega);

    CGAL::Orientation s = oriented_height_sign(a2, b2, c2, a3, b3, c3, d3);

    // information about this edge
    EdgeRecord rec;
    rec.ia = vidx.at(va); // polyscope vertex index of edge endpoint
    rec.ib = vidx.at(vb); // polyscope vertex index of edge endpoint
    rec.a2 = a2; rec.b2 = b2; // store 2D coordinates of vertices for debug
    rec.nonreg_down = (s == CGAL::NEGATIVE);  // d' below plane a'b'c' -> down flip
    rec.nonreg_up   = (s == CGAL::POSITIVE); // d' above plane a'b'c' -> up flip
    rec.degenerate  = (s == CGAL::COLLINEAR); // coplanar in lifted test
    
    // list of flippable edges with flags, vertex IDs, and 2D coordinates
    out.push_back(rec);
  }
  return out;
}

// print table of those edges 
void print_flippable_edges_table_ids(const Tri& tri, const OmegaQuad& omega, std::ostream& os) {
  auto edges = collect_flippable_edges_with_ids(tri, omega);

  using std::left;
  using std::setw;

  os << "\n=== locally non-regular internal edges (planar IDs match polyscope mesh) ===\n";
  os << left
     << setw(6)  << "start"
     << setw(6)  << "end"
     << setw(14) << "start_x"
     << setw(14) << "start_y"
     << setw(14) << "end_x"
     << setw(14) << "end_y"
     << setw(8)  << "up"
     << setw(8)  << "down"
     << setw(11) << "deg"
     << "\n";

  os << std::boolalpha << std::fixed << std::setprecision(6);

  for (const auto& e : edges) {
    os << left
       << setw(6)  << e.ia
       << setw(6)  << e.ib
       << setw(14) << e.a2.x()
       << setw(14) << e.a2.y()
       << setw(14) << e.b2.x()
       << setw(14) << e.b2.y()
       << setw(8)  << e.nonreg_up
       << setw(8)  << e.nonreg_down
       << setw(11) << e.degenerate
       << "\n";
  }

  if (edges.empty()) os << "(none)\n";
  os << "=== end ===\n";
}


// Optional: print map and counts to sanity-check you’re on the same triangulation.
void print_triangulation_summary(const Tri& tri) {
  size_t nV = static_cast<size_t>(std::distance(tri.finite_vertices_begin(), tri.finite_vertices_end()));
  size_t nF = static_cast<size_t>(std::distance(tri.finite_faces_begin(), tri.finite_faces_end()));
  std::cout << "Triangulation summary: " << nV << " vertices, " << nF << " faces.\n";
  size_t idx = 0;
  for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit, ++idx) {
    auto p = vit->point();
    std::cout << "  id " << idx << " -> (" << p.x() << ", " << p.y() << ")\n";
  }
}

} // namespace reg
} // namespace df



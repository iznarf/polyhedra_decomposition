#include "check_edges.h"

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


  // function 1
  //takes input2D struct and checks every internal edge in the triangulation object 
  // first it checks if the quad around the edge is strictly convex in 2D
  // then it lifts the four points to 3D using omega
  // then it checks the orientation of the lifted points to determine if the edge is flippable or not
  // it returns a list of edges that are locally non-regular (so the resulting flip would be a down flip) 
  // we return couples of vertex ids (global indices) for each edge
  // since we have a map global id -> vertex handle in the input2D struct
  // then use vertex handles and is_edge(vertex_handle a, vertex_handle b) function of CGAL triangulation to get the face handle of incident face and opposite vertex of the edge contained in face
  // we need face handle and the index for flipping then


  // function 2 
  // we write another function that checks all locally non regular edges if they are already in the target triangulation
  // here we can again use the global indices to get the vertex handles of the target triangulation
  // then we can simply check with is_edge(vertex_handle a, vertex_handle b) function of CGAL triangulation
  // if the edge is already in the target triangulation we can delete it form the list of edges to flip




  // ensure (a,b,c) ccw in 2D, then evaluate orientation(a',b',c',d')
  // returns POSITIVE if d' above plane a'b'c', NEGATIVE if below, COLLINEAR if coplanar
  static inline CGAL::Orientation oriented_height_sign(const P2& a2, const P2& b2, const P2& c2,
                                                      const P3& a3, const P3& b3, const P3& c3, const P3& d3) {
    auto o2 = CGAL::orientation(a2, b2, c2);
    if (o2 == CGAL::RIGHT_TURN) return CGAL::orientation(a3, c3, b3, d3);
    if (o2 == CGAL::COLLINEAR)  return CGAL::COLLINEAR;
    return CGAL::orientation(a3, b3, c3, d3);
  }

  // quad strictly convex around edge ab: c and d on opposite sides of ab, neither collinear
  static inline bool quad_strictly_convex(const P2& a, const P2& b, const P2& c, const P2& d) {
    auto oc = CGAL::orientation(a, b, c);
    auto od = CGAL::orientation(a, b, d);
    return (oc != CGAL::COLLINEAR) && (od != CGAL::COLLINEAR) && (oc != od);
  }

  // function 1
  std::vector<std::array<df::vertex_id, 2>> find_locally_non_regular_edges(const Tri2& tri) {
      std::vector<std::array<df::vertex_id, 2>> non_regular_edges;
      // iterate over all finite edges in the triangulation
      for (auto edge = tri.finite_edges_begin(); edge != tri.finite_edges_end(); ++edge) {
          auto f = edge->first; // incident face of the edge
          int  i = edge->second; // index of opposite vertex of edge in face f 

          // internal edge?
          auto g = f->neighbor(i); // opposite face across edge
          if (tri.is_infinite(g)) continue; // if g is infinite, it is a boundary edge. we skip it since we can just flip internal edges 

          // extract the quad around the edge  
          auto va = f->vertex(tri.cw(i)); // vertex handle of endpoint of edge 
          auto vb = f->vertex(tri.ccw(i)); // vertex handle of endpoint of edge 
          auto vc = f->vertex(i); // vertex handle of opposite vertex of edge in face f 
          int  j  = tri.mirror_index(f, i); // index of opposite vertex of edge in face g
          auto vd = g->vertex(j); // vertex handle of vertex opposite to edge in face g

          // extract 2D points of the quad around the edge 
          const P2 &a2 = va->point(), &b2 = vb->point(), &c2 = vc->point(), &d2 = vd->point();

          // check if the quad is convex
          if (!quad_strictly_convex(a2, b2, c2, d2)) continue;

          // lift points to 3D
          P3 a3(a2.x(), a2.y(), a2.x()*a2.x() + a2.y()*a2.y());
          P3 b3(b2.x(), b2.y(), b2.x()*b2.x() + b2.y()*b2.y());
          P3 c3(c2.x(), c2.y(), c2.x()*c2.x() + c2.y()*c2.y());
          P3 d3(d2.x(), d2.y(), d2.x()*d2.x() + d2.y()*d2.y());


          CGAL::Orientation s = oriented_height_sign(a2, b2, c2, a3, b3, c3, d3);

          // down flip?
          if (s == CGAL::NEGATIVE) {
              df::vertex_id ia = va->info();
              df::vertex_id ib = vb->info();
              non_regular_edges.push_back({ia, ib});
          }
      }
      return non_regular_edges;
  }

// function 2
bool edge_in_target(const df::vertex_id ia, const df::vertex_id ib, const df::InputData& D) {
    // iterate over all finite edges in the target triangulation
    const df::Tri2& target = D.tri_lower;  
    for (auto e = target.finite_edges_begin(); e != target.finite_edges_end(); ++e) {
        auto f = e->first;
        int  i = e->second;

        auto va = f->vertex(target.cw(i));
        auto vb = f->vertex(target.ccw(i));

        df::vertex_id ja = va->info();
        df::vertex_id jb = vb->info();

        // undirected edge comparison
        if ((ja == ia && jb == ib) || (ja == ib && jb == ia)) {
            return true;
        }
    }
    return false;
}
} // namespace reg
} // namespace df



#include "check_edges.h"
#include "geometry_utils.h"

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

    using df::oriented_height_sign;


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
  // function that checks all locally non regular edges if they are already in the target triangulation
  // here we can again use the global indices to get the vertex handles of the target triangulation
  // then we can simply check with is_edge(vertex_handle a, vertex_handle b) function of CGAL triangulation
  // if the edge is already in the target triangulation we can delete it form the list of edges to flip





  // function 1
  std::vector<std::array<df::vertex_id, 2>> find_locally_non_regular_edges(const Tri2& tri) {
        std::vector<std::array<df::vertex_id, 2>> non_regular_edges;

    for (auto edge = tri.finite_edges_begin(); edge != tri.finite_edges_end(); ++edge) {
        auto f = edge->first;  // incident face of the edge
        int  i = edge->second; // index of opposite vertex in face f

        // internal edge? both incident faces must be finite
        if (tri.is_infinite(f)) continue;
        auto g = f->neighbor(i);
        if (tri.is_infinite(g)) continue;

        // extract the quad around the edge
        auto va = f->vertex(tri.cw(i));
        auto vb = f->vertex(tri.ccw(i));
        auto vc = f->vertex(i);
        int  j  = tri.mirror_index(f, i);
        auto vd = g->vertex(j);

        const P2 &a2 = va->point(), &b2 = vb->point(), &c2 = vc->point(), &d2 = vd->point();

        if (!df::quad_strictly_convex(a2, b2, c2, d2)) continue;

        // lift to 3D using lift in geometry_utils
        P3 a3 = lift(a2);
        P3 b3 = lift(b2);
        P3 c3 = lift(c2);
        P3 d3 = lift(d2);

       
        CGAL::Orientation s = oriented_height_sign(a2, b2, c2, a3, b3, c3, d3);

        // if negative then edge (a,b) is locally non-regular -> down flip
        if (s == CGAL::NEGATIVE) {
            df::vertex_id ia = va->info();
            df::vertex_id ib = vb->info();
            non_regular_edges.push_back({ia, ib});
        }

        if (s == CGAL::COPLANAR) {
            // just print that there is a colplanar quad 
            df::vertex_id ia = va->info();
            df::vertex_id ib = vb->info();
            std::cout << "[check_edges] WARNING: edge (" << ia << "," << ib
                      << ") is in a coplanar quad\n";
        }
    }
    return non_regular_edges;
}


// function 2
bool edge_in_target(const df::vertex_id ia, const df::vertex_id ib, const df::InputData& D) {
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



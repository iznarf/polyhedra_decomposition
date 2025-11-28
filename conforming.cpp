#include "conforming.h"

#include <CGAL/Segment_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>
#include <algorithm>
#include <limits>
#include <optional>
#include <variant>
#include <vector>

namespace df { namespace reg {

using K   = df::K;    
using P2  = df::P2;
using P3  = df::P3;
using Tri = df::Tri2;
using Seg3  = CGAL::Segment_3<K>;


static inline P3 lift(const P2& p) {
  return P3(p.x(), p.y(), p.x()*p.x() + p.y()*p.y());
}

// we check if flipping the edge with endpoints ia, ib is conforming to the target triangulation
// conforming to the target triangulation means that the tetrahedron formed by the lifted points a,b,c,d does not intersect the faces of the target triangulation


bool is_flip_conforming(df::vertex_id ia, df::vertex_id ib, const df::InputData& D)
{
    
    // find the edge with global indices (ia, ib) in the current triangulation by vertex ids 
    Tri::Face_handle fh;
    int i = -1;
    bool found = false;
    
    for (auto e = D.tri_current.finite_edges_begin(); e != D.tri_current.finite_edges_end(); ++e)
    {
        auto f  = e->first; // incident face of the edge
        int ei  = e->second; // local index of the edge in face f

        auto va = f->vertex(D.tri_current.cw(ei)); // vertex opposite to index ei in face f
        auto vb = f->vertex(D.tri_current.ccw(ei)); // vertex opposite to index ei in face f

        df::vertex_id ja = va->info(); // global index of vertex a
        df::vertex_id jb = vb->info(); // global index of vertex b

        // we check if indices match, if yes then we found the edge in the current triangulation
        if ((ja == ia && jb == ib) || (ja == ib && jb == ia)) {
            fh    = f;
            i     = ei;
            found = true;
            break;
        }
    }

    if (!found) {
        std::cout << "[conform] ERROR: edge (" << ia << "," << ib
                  << ") not found in tri_current\n";
        return false;
    }

    // neighbor face across the edge
    auto gh = fh->neighbor(i);

    // if we hit a boundary edge (which should not happen), treat as non-conforming
    if (D.tri_current.is_infinite(fh) || D.tri_current.is_infinite(gh)) {
        std::cout << "[conform] WARNING: edge (" << ia << "," << ib
                  << ") is on boundary / infinite, treating as non-conforming\n";
        return false; 
    }

    // opposite vertices c (in fh) and d (in gh)
    auto vc = fh->vertex(i);
    int  j  = D.tri_current.mirror_index(fh, i);
    auto vd = gh->vertex(j);

    df::vertex_id ic = vc->info();
    df::vertex_id id = vd->info();


    // candidate (c,d)
    P2 c2 = vc->point();
    P2 d2 = vd->point();
    const CGAL::Segment_2<K> edge_cd_2d(c2, d2);

    // scan edges in LOWER triangulation
    for (auto e = D.tri_lower.finite_edges_begin(); e != D.tri_lower.finite_edges_end(); ++e) {
        auto f  = e->first;
        int  ei = e->second;

        auto vu = f->vertex(D.tri_lower.cw(ei));
        auto vv = f->vertex(D.tri_lower.ccw(ei));
        df::vertex_id iu = vu->info();
        df::vertex_id iv = vv->info();


        // if (c,d) is already in lower triangulation the flip is conforming
        if (((iu == ic) && (iv == id)) || ((iu == id) && (iv == ic))) {
            std::cout << "[conform] (c,d)=(" << ic << "," << id << ") is already in lower -> flip of edge (" << ia << "," << ib << ") is conforming \n";
            return true;
        }

        // skip if edge shares endpoint with (c,d) ("intersection" at endpoint)
        if (iu == ic || iu == id || iv == ic || iv == id) continue;

        // get coordinates of vertices of edge (u,v) in lower triangulation
        P2 u2 = vu->point();
        P2 v2 = vv->point();

        // 2D intersection
        const CGAL::Segment_2<K> edge_uv_2d(u2, v2);

        // if they intersect in 2D and interesction is not at endpoint, we have to check 3D intersection
        if (CGAL::do_intersect(edge_cd_2d, edge_uv_2d)) {
            // lift points to 3D
            P3 c3 = lift(c2);
            P3 d3 = lift(d2);
            P3 u3 = lift(u2);
            P3 v3 = lift(v2);

            // check if edges intersect in 3D
            // if they intersect in 3D we have to check if the segments lie in one plane
            // if yes, the flip is conforming since this is only boundary contact which is allowed 
            // if not, the flip is non-conforming
            if (CGAL::do_intersect(Seg3(c3, d3), Seg3(u3, v3))) {
                if (CGAL::orientation(c3, d3, u3, v3) == CGAL::COPLANAR) {
                    // segments intersect in 3D but are coplanar -> allowed
                    std::cout << "[conform] PASS: segments intersect in 3D but are coplanar with lower (u,v)=(" << iu << "," << iv
                              << ")\n";
                    continue;
                }
                else {
                    // segments intersect in 3D and are not coplanar -> non-conforming
                    return false;
                }
            }

            // (u,v) intersects in 2D but not in 3D -> check height orientation
            // for height orientation we need a consistent ordering of (u,v) relative to (c,d)
            // we say u is left of (c,d) and v is right of (c,d)
            auto su = CGAL::orientation(c2, d2, u2); // >0: u is LEFT of c->d, <0: RIGHT
            if (su == CGAL::COLLINEAR) continue;
            if (su == CGAL::NEGATIVE) {
                std::swap(u2, v2);
                std::swap(u3, v3);
            }

            // orientation-based height test: we check if v3 lies in direction of normal of triangle (c3,d3,u3)
            const auto o = CGAL::orientation(c3, d3, u3, v3);
            // CGAL::orientation is positive if v3 is in normal direction of plane (c3,d3,u3)
            // it is negative if v3 is not in normal direction
            // block -> (c,d) below (u,v); v3 lies in opposite direction of normal
            if (o == CGAL::POSITIVE) {
                std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
                        << "below lower (u,v)=(" << iu << "," << iv << ")\n";
                return false;
            }
            
            // // (c,d) not below (u,v) -> pass this lower edge; v3 lies in normal direction
            if (o == CGAL::NEGATIVE || o == CGAL::COPLANAR) {
                std::cout << "[conform] PASS: (c,d)=(" << ic << "," << id << ") "
                        << "above lower (u,v)=(" << iu << "," << iv << ")\n";
                continue;
            }
        }
    }

    //std::cout << "[conform] OK: flip (a,b)=(" << ia << "," << ib << ") -> (c,d)=("
    //          << ic << "," << id << ") is conforming\n";
    return true;
}


} } // namespace df::reg

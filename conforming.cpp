#include "conforming.h"
#include "height_test.h"
#include "geometry_utils.h"

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


// we check if flipping the edge with endpoints ia, ib is conforming to the target triangulation
// conforming to the target triangulation means that the tetrahedron formed by the lifted points a,b,c,d does not intersect the faces of the target triangulation


bool is_flip_conforming(df::vertex_id ia, df::vertex_id ib, const df::InputData& D, const df::Tri2& tri_current) {
    
    // print which edge we are checking
    std::cout << "[conform] checking flip of edge (" << ia << "," << ib << ")\n";
    // find the edge with global indices (ia, ib) in the current triangulation by vertex ids 
    Tri::Face_handle fh;
    int i = -1;
    bool found = false;
    
    for (auto e = tri_current.finite_edges_begin(); e != tri_current.finite_edges_end(); ++e) {
        auto f  = e->first; // incident face of the edge
        int ei  = e->second; // local index of the edge in face f

        auto va = f->vertex(tri_current.cw(ei)); // vertex opposite to index ei in face f
        auto vb = f->vertex(tri_current.ccw(ei)); // vertex opposite to index ei in face f

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
    if (tri_current.is_infinite(fh) || tri_current.is_infinite(gh)) {
        std::cout << "[conform] WARNING: edge (" << ia << "," << ib
                  << ") is on boundary / infinite, treating as non-conforming\n";
        return false; 
    }

    // opposite vertices c (in fh) and d (in gh)
    auto vc = fh->vertex(i);
    int  j  = tri_current.mirror_index(fh, i);
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
            P3 c3 = df::lift(c2);
            P3 d3 = df::lift(d2);
            P3 u3 = df::lift(u2);
            P3 v3 = df::lift(v2);

            // check if edges intersect in 3D
            // if they intersect in 3D we have to check if the segments lie in one plane
            // if yes, the flip is conforming since this is only boundary contact which is allowed 
            // if not, the flip is non-conforming
            if (CGAL::do_intersect(Seg3(c3, d3), Seg3(u3, v3))) {
                // print that segments intersect in 3D and that this is not allowed
                std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
                          << "intersects lower (u,v)=(" << iu << "," << iv << ") in 3D\n";
                return false;
            }

           

            auto height = oriented_height_sign(u2, v2, c2, u3, v3, c3, d3);
            int cmp = compare_heights_at_intersection(c2, d2, u2, v2, c3, d3, u3, v3);
           
        
            if (height < 0){
                std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
                        << "below lower (u,v)=(" << iu << "," << iv << ")\n";
                // also print cmp value
                std::cout << "cmp = " << cmp << "\n";
                return false;
            }
            if (height > 0){
                std::cout << "[conform] PASS: (c,d)=(" << ic << "," << id << ") "
                        << "above lower (u,v)=(" << iu << "," << iv << ")\n";
                std::cout << "cmp = " << cmp << "\n";
                continue;
            }
            
        }
    }

    //std::cout << "[conform] OK: flip (a,b)=(" << ia << "," << ib << ") -> (c,d)=("
    //          << ic << "," << id << ") is conforming\n";
    return true;
}


} } // namespace df::reg

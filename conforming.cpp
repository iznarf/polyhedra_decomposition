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


bool is_flip_conforming(df::vertex_id ia,
                        df::vertex_id ib,
                        const df::InputData& D,
                        const std::unordered_map<df::vertex_id,int>& li_current,
                        const std::unordered_map<df::vertex_id,int>& li_lower)
{
    auto idxStr = [](const auto& M, df::vertex_id g)->std::string {
        auto it = M.find(g);
        return (it == M.end()) ? "nc" : std::to_string(it->second);
    };
    auto Lcur = [&](df::vertex_id g){ return idxStr(li_current, g); };
    auto Llow = [&](df::vertex_id g){ return idxStr(li_lower,   g); };

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

    //std::cout << "[conform] check edge (a,b)=(" << ia << "," << ib << ")"
    //          << "  polyscope indices =(" << Lcur(ia) << "," << Lcur(ib) << ")\n";


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

        // the edge (a,b) will always intersect edge (c,d)
        //if ((iu == ia && iv == ib) || (iu == ib && iv == ia)) continue;

        // if (c,d) is already in lower triangulation the flip is conforming
        // maybe we have to leave this out? 
        if (((iu == ic) && (iv == id)) || ((iu == id) && (iv == ic))) {
            std::cout << "[conform] (c,d)=(" << ic << "," << id << ") is already in lower -> flip is conforming \n";
            return true;
        }

        // skip if edge shares endpoint with (c,d) ("intersection" at endpoint)
        if (iu == ic || iu == id || iv == ic || iv == id) continue;

        // get coordinates of vertices of edge (u,v) in lower triangulation
        P2 u2 = vu->point();
        P2 v2 = vv->point();

        // 2D intersection
        const CGAL::Segment_2<K> edge_uv_2d(u2, v2);

        // we want to avoid to exactly compute intersectin objects 
        // so do intersect would be better here and then do not compute the intersection point because we do not need it
        // fix it later (check if it still works correctly after relaxing the intersection test)
        auto intersection_object = CGAL::intersection(edge_cd_2d, edge_uv_2d);
        // if (c,d) and (u,v) do not intersect in 2D, continue
        if (!intersection_object) continue;

        // treat segment result: if non-degenerate, it's overlap -> skip (boundary contact)
        P2 x;

        if (const CGAL::Segment_2<K>* s = std::get_if<CGAL::Segment_2<K>>(&*intersection_object)) {
            // check squared difference of endpoints of segment to make sure it is really a segment
            const P2& x0 = s->source(); const P2& x1 = s->target();
            const double dx = CGAL::to_double(x1.x() - x0.x());
            const double dy = CGAL::to_double(x1.y() - x0.y());
            // if distance is smaller than threshold, treat as point intersection
            // intersection point is then x0 
            if (dx*dx + dy*dy < 1e-24) x = x0; else continue;
            // otherwise treat as overlap -> skip
        } else if (const P2* p = std::get_if<P2>(&*intersection_object)) {
            x = *p;
        } else {
            continue;
        }

        // check for endpoint touch or outside segment
        // we want proper intersection in the interior of both segments
        // we only do this because maybe due to numerical issues the intersection test is slightly off the segment? 
        const bool cd_strict = CGAL::collinear_are_strictly_ordered_along_line(c2, x, d2);
        const bool uv_strict = CGAL::collinear_are_strictly_ordered_along_line(u2, x, v2);
        if (!(cd_strict && uv_strict)) continue;  

        // lift points to 3D
        P3 c3 = lift(c2);
        P3 d3 = lift(d2);
        P3 u3 = lift(u2);
        P3 v3 = lift(v2);

        // check if edges intersect in 3D
        // if they intersect in 3D the flip is non-conforming 
        if (CGAL::do_intersect(Seg3(c3, d3), Seg3(u3, v3))) {
            //std::cout << "[conform] BLOCK: segments intersect in 3D with lower (u,v)=(" << iu << "," << iv
            //          << ") local indices =(" << Llow(iu) << "," << Llow(iv) << ")\n";
            return false;
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
        // (c,d) not below (u,v) -> pass this lower edge; v3 lies in normal direction
        if (o == CGAL::POSITIVE) {
            //std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
            //        << "below lower (u,v)=(" << iu << "," << iv << ")\n";
            return false;
        }
        // block -> (c,d) below (u,v); v3 lies in opposite direction of normal
        if (o == CGAL::NEGATIVE) {
            //std::cout << "[conform] PASS: (c,d)=(" << ic << "," << id << ") "
            //        << "above lower (u,v)=(" << iu << "," << iv << ")\n";
            continue;
        }
    }

    //std::cout << "[conform] OK: flip (a,b)=(" << ia << "," << ib << ") -> (c,d)=("
    //          << ic << "," << id << ") is conforming\n";
    return true;
}


} } // namespace df::reg


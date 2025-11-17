#include "conforming_insertion.h"
#include <CGAL/Segment_3.h>
#include <CGAL/Segment_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Object.h>
#include <algorithm>
#include <limits>
#include <optional>
#include <variant>
#include <vector>
#include <iostream>
#include <unordered_map>

namespace df { namespace reg {

using K   = df::K;
using P2  = df::P2;
using P3  = df::P3;
using Tri = df::Tri2;
using Seg3 = CGAL::Segment_3<K>;
using Seg2 = CGAL::Segment_2<K>;



static inline P3 lift(const P2& p) {
    return P3(p.x(), p.y(), p.x()*p.x() + p.y()*p.y());
}

bool is_insertion_conforming(df::vertex_id id,
                             const df::InputData& D)
{

    const Tri& current = D.tri_current;

    // 2D position of the vertex to be inserted
    const P2& d2 = D.points2d[id];

    // locate in current triangulation
    Tri::Locate_type lt;
    int li; // index used if point lies on an edge
    Tri::Face_handle fh = current.locate(d2, lt, li);

    if (lt != Tri::FACE) {
        // p is not inside a face, we do not handle this case here
        std::cout << "[conform] WARNING: insertion point is not inside a face, cannot check conforming\n";
        return false;
    }

    // face containing the point to be inserted
    auto va = fh->vertex(0);
    auto vb = fh->vertex(1);
    auto vc = fh->vertex(2);

    df::vertex_id ia = va->info();
    df::vertex_id ib = vb->info();
    df::vertex_id ic = vc->info();

    const P2& a2 = va->point();
    const P2& b2 = vb->point();
    const P2& c2 = vc->point();

    // edges that will be created by inserting d
    Seg2 edge_ad_2d(a2, d2); // (id, ia)
    Seg2 edge_bd_2d(b2, d2); // (id, ib)
    Seg2 edge_cd_2d(c2, d2); // (id, ic)



    // iterate over all edges of the lower triangulation
    for (auto e = D.tri_lower.finite_edges_begin();
         e != D.tri_lower.finite_edges_end(); ++e) {

        auto f  = e->first;
        int  ei = e->second;

        auto vu = f->vertex(D.tri_lower.cw(ei));
        auto vv = f->vertex(D.tri_lower.ccw(ei));

        df::vertex_id u = vu->info();
        df::vertex_id v = vv->info();

        // skip neighbor edges (sharing vertices with a,b,c)
        if (u == ia || u == ib || u == ic || u == id ||
            v == ia || v == ib || v == ic || v == id) {
            continue;
        }

        // **Important change**: make *copies* so we can swap them
        P2 u2 = vu->point();
        P2 v2 = vv->point();

        Seg2 edge_uv_2d(u2, v2);

        bool intersection_ad = CGAL::do_intersect(edge_ad_2d, edge_uv_2d);
        bool intersection_bd = CGAL::do_intersect(edge_bd_2d, edge_uv_2d);
        bool intersection_cd = CGAL::do_intersect(edge_cd_2d, edge_uv_2d);

        // if there is no intersection, nothing to check for this lower edge
        if (!intersection_ad && !intersection_bd && !intersection_cd) {
            continue;
        }

        // if any of the segments in 3D intersect, the insertion is non-conforming
        // lift all points to 3D and make segments
        P3 a3 = lift(a2);
        P3 b3 = lift(b2);
        P3 c3 = lift(c2);
        P3 d3 = lift(d2);
        P3 u3 = lift(u2);
        P3 v3 = lift(v2);

        // create 3D segments
        Seg3 seg_ad_3d(a3, d3);
        Seg3 seg_bd_3d(b3, d3);
        Seg3 seg_cd_3d(c3, d3);
        Seg3 seg_uv_3d(u3, v3);

        if (CGAL::do_intersect(seg_ad_3d, seg_uv_3d)) {
            std::cout << "[conform] WARNING: inserting vertex with global id "
                      << id << " is non-conforming (edge (a,d))\n";
            return false;
        }
        if (CGAL::do_intersect(seg_bd_3d, seg_uv_3d)) {
            std::cout << "[conform] WARNING: inserting vertex with global id "
                      << id << " is non-conforming (edge (b,d))\n";
            return false;
        }
        if (CGAL::do_intersect(seg_cd_3d, seg_uv_3d)) {
            std::cout << "[conform] WARNING: inserting vertex with global id "
                      << id << " is non-conforming (edge (c,d))\n";
            return false;
        }



        // For each intersecting edge, check the height condition

        // (a,d) vs (u,v)
        if (intersection_ad) {
            // we need the copies for consistent orientation tests
            // make copies again
            P2 u2 = vu->point();
            P2 v2 = vv->point();
            P3 u3 = lift(u2);
            P3 v3 = lift(v2);

            auto su = CGAL::orientation(a2, d2, u2); // >0: u left of a->d, <0: right
            // if collinear, we just skip orientation-based test
            if (su == CGAL::NEGATIVE) {
                std::swap(u2, v2);
                std::swap(u3, v3);            
            }

            const auto o = CGAL::orientation(a3, d3, u3, v3);
            // interpret: if (a,d) is "below" (u,v), insertion is non-conforming
            if (o == CGAL::NEGATIVE) {
                std::cout << "[conform] WARNING: inserting vertex with global id "
                          << id << " is non-conforming (edge (a,d))\n";
                return false;
            }
        }

        // (b,d) vs (u,v)
        if (intersection_bd) {
            P2 u2 = vu->point();
            P2 v2 = vv->point();
            P3 u3 = lift(u2);
            P3 v3 = lift(v2);

            auto su = CGAL::orientation(b2, d2, u2);
            if (su == CGAL::NEGATIVE) {
                std::swap(u2, v2);
                std::swap(u3, v3);            
            }

            const auto o = CGAL::orientation(b3, d3, u3, v3);
            if (o == CGAL::NEGATIVE) {
                std::cout << "[conform] WARNING: inserting vertex with global id "
                          << id << " is non-conforming (edge (b,d))\n";
                return false;
            }
        }

        // (c,d) vs (u,v)
        if (intersection_cd) {
            P2 u2 = vu->point();
            P2 v2 = vv->point();
            P3 u3 = lift(u2);
            P3 v3 = lift(v2);
           
            auto su = CGAL::orientation(c2, d2, u2);
            if (su == CGAL::NEGATIVE) {
                std::swap(u2, v2);
                std::swap(u3, v3);
            }

            const auto o = CGAL::orientation(c3, d3, u3, v3);
            if (o == CGAL::NEGATIVE) {
                std::cout << "[conform] WARNING: inserting vertex with global id "
                          << id << " is non-conforming (edge (c,d))\n";
                return false;
            }
        }
    }

    // if we never found a blocking lower edge, insertion is conforming
    return true;
}

}} // namespace df::reg






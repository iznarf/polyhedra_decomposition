
#include "conforming_insertion.h"
#include "input.h"
#include "height_test.h"
#include "geometry_utils.h"
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




// fix: one for loop not three
bool is_insertion_conforming(df::vertex_id id, const df::InputData& D, const df::Tri2& current) {
  
    const Tri& lower   = D.tri_lower;

    // 2D position of the vertex to be inserted
    const P2& d2 = D.points2d[id];
    const P3& d3 = df::lift(d2);

    // locate in current triangulation
    Tri::Locate_type lt;
    int li; // index used if point lies on an edge -> will not happen
    Tri::Face_handle fh = current.locate(d2, lt, li);

    if (lt != Tri::FACE) {
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

    // fast path: are (a,d), (b,d), (c,d) already in the lower triangulation 
    // if yes, we flip is conforming
    bool ad_in_lower = false;
    bool bd_in_lower = false;
    bool cd_in_lower = false;


    for (auto e = lower.finite_edges_begin(); e != lower.finite_edges_end(); ++e) {
        auto f  = e->first;
        int  ei = e->second;

        auto vu = f->vertex(lower.cw(ei));
        auto vv = f->vertex(lower.ccw(ei));

        df::vertex_id u = vu->info();
        df::vertex_id v = vv->info();

        if ((u == ia && v == id) || (u == id && v == ia))
            ad_in_lower = true;
        if ((u == ib && v == id) || (u == id && v == ib))
            bd_in_lower = true;
        if ((u == ic && v == id) || (u == id && v == ic))
            cd_in_lower = true;

    }

    if (ad_in_lower && bd_in_lower && cd_in_lower) {
        std::cout << "[conform] all three edges (" << ia << "," << id << "), ("
                  << ib << "," << id << "), (" << ic << "," << id
                  << ") are already in lower triangulation -> insertion is conforming\n";
        return true;
    }

    // lift the triangle vertices
    P3 a3 = df::lift(a2);
    P3 b3 = df::lift(b2);
    P3 c3 = df::lift(c2);

    Seg2 edge_ad_2d(a2, d2);
    Seg2 edge_bd_2d(b2, d2);
    Seg2 edge_cd_2d(c2, d2);

    // edge (a,d)
    if (ad_in_lower == false) {

    
        for (auto e = lower.finite_edges_begin(); e != lower.finite_edges_end(); ++e) {
            auto f  = e->first;
            int  ei = e->second;

            auto vu = f->vertex(lower.cw(ei));
            auto vv = f->vertex(lower.ccw(ei));

            df::vertex_id u = vu->info();
            df::vertex_id v = vv->info();

            // skip edges that are incident to d or incident to a ("intersection" at endpoint)
            if (u == ia || u == id || v == ia || v == id) continue;

            P2 u2 = vu->point();
            P2 v2 = vv->point();

            Seg2 edge_uv_2d(u2, v2);

            // cheap 2D filter: if segments don't even intersect in 2D, ignore
            if (CGAL::do_intersect(edge_ad_2d, edge_uv_2d) == true){ 

                // full 3D test
                P3 u3 = df::lift(u2);
                P3 v3 = df::lift(v2);

                Seg3 seg_uv_3d(u3, v3);
                Seg3 seg_ad_3d(a3, d3);

                // if they intersect in 3D, flip is not conforming
                if (CGAL::do_intersect(seg_ad_3d, seg_uv_3d)) {
                    std::cout << "[conform] BLOCK: (a,d)=(" << ia << "," << id << ") "
                            << "intersects lower (u,v)=(" << u << "," << v << ") in 3D\n";

                    return false;
                }

                int cmp_ad = compare_heights_at_intersection(a2, d2, u2, v2, a3, d3, u3, v3);

                auto height_ad = oriented_height_sign(u2, v2, a2, u3, v3, a3, d3);

       
                if (height_ad < 0){
                    std::cout << "[conform] BLOCK: (a,d)=(" << ia << "," << id << ") "
                            << "below lower (u,v)=(" << u << "," << v << ")\n";
                    std::cout << "cmp_ad = " << cmp_ad << "\n";
                    return false;
                }
                
                if (height_ad > 0){
                    std::cout << "[conform] PASS: (a,d)=(" << ia << "," << id << ") "
                            << "above lower (u,v)=(" << u << "," << v << ")\n";
                    std::cout << "cmp_ad = " << cmp_ad << "\n";
                }
                
        
            }
        }
    }

    // edge (b,d)
    if (bd_in_lower == false) {
   
        for (auto e = lower.finite_edges_begin(); e != lower.finite_edges_end(); ++e) {
            auto f  = e->first;
            int  ei = e->second;

            auto vu = f->vertex(lower.cw(ei));
            auto vv = f->vertex(lower.ccw(ei));

            df::vertex_id u = vu->info();
            df::vertex_id v = vv->info();


            // skip edges that are incident to d or incident to b
            if (u == ib || u == id || v == ib || v == id) continue;

            P2 u2 = vu->point();
            P2 v2 = vv->point();

            Seg2 edge_uv_2d(u2, v2);

            if (CGAL::do_intersect(edge_bd_2d, edge_uv_2d) == true){ 

                P3 u3 = df::lift(u2);
                P3 v3 = df::lift(v2);

                Seg3 seg_uv_3d(u3, v3);
                Seg3 seg_bd_3d(b3, d3);

                if (CGAL::do_intersect(seg_bd_3d, seg_uv_3d)) {
                    // print that segments intersect in 3D and that this is not allowed
                    std::cout << "[conform] BLOCK: (b,d)=(" << ib << "," << id << ") "
                            << "intersects lower (u,v)=(" << u << "," << v << ") in 3D\n";
                    return false;
                }

                int cmp_bd = compare_heights_at_intersection(b2, d2, u2, v2, b3, d3, u3, v3);
                auto height_bd = oriented_height_sign(u2, v2, b2, u3, v3, b3, d3);
            
                if (height_bd < 0){
                    std::cout << "[conform] BLOCK: (b,d)=(" << ib << "," << id << ") "
                            << "below lower (u,v)=(" << u << "," << v << ")\n";
                    std::cout << "cmp_bd = " << cmp_bd << "\n";
                    return false;
                }
                
                if (height_bd > 0){
                    std::cout << "[conform] PASS: (b,d)=(" << ib << "," << id << ") "
                            << "above lower (u,v)=(" << u << "," << v << ")\n";
                    std::cout << "cmp_bd = " << cmp_bd << "\n";
                }
                
            }
        }
    }
    // edge (c,d)
    if (cd_in_lower == false) {
        for (auto e = lower.finite_edges_begin(); e != lower.finite_edges_end(); ++e) {
            auto f  = e->first;
            int  ei = e->second;

            auto vu = f->vertex(lower.cw(ei));
            auto vv = f->vertex(lower.ccw(ei));

            df::vertex_id u = vu->info();
            df::vertex_id v = vv->info();

            // skip edges that are incident to d or incident to c
            if (u == ic || u == id || v == ic || v == id) continue;

            P2 u2 = vu->point();
            P2 v2 = vv->point();

            Seg2 edge_uv_2d(u2, v2);

            if (CGAL::do_intersect(edge_cd_2d, edge_uv_2d) == true){

                P3 u3 = lift(u2);
                P3 v3 = lift(v2);

                Seg3 seg_uv_3d(u3, v3);
                Seg3 seg_cd_3d(c3, d3);

                if (CGAL::do_intersect(seg_cd_3d, seg_uv_3d)) {
                    std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
                            << "intersects lower (u,v)=(" << u << "," << v << ") in 3D\n";
                    return false;
                }

                int cmp_cd = compare_heights_at_intersection(c2, d2, u2, v2, c3, d3, u3, v3);
                auto height_cd = oriented_height_sign(u2, v2, c2, u3, v3, c3, d3);
              
                if (height_cd < 0){
                    std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
                            << "below lower (u,v)=(" << u << "," << v << ")\n";
                    std::cout << "cmp_cd = " << cmp_cd << "\n";
                    return false;
                }
                if (height_cd > 0){
                    std::cout << "[conform] PASS: (c,d)=(" << ic << "," << id << ") "
                            << "above lower (u,v)=(" << u << "," << v << ")\n";
                    std::cout << "cmp_cd = " << cmp_cd << "\n";
                }
            }  
        }
    }

    std::cout << "[conform] insertion of vertex with global id " << id
              << " is conforming\n";
    return true;
}
}}
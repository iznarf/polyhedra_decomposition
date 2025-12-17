
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

    // vertex handles of face in ccw order in current triangulation containing the point to be inserted
    auto va = fh->vertex(0);
    auto vb = fh->vertex(1);
    auto vc = fh->vertex(2);

    // global ids of these vertices
    df::vertex_id ia = va->info();
    df::vertex_id ib = vb->info();
    df::vertex_id ic = vc->info();

    // 2D points of these vertices
    const P2& a2 = va->point();
    const P2& b2 = vb->point();
    const P2& c2 = vc->point();

    // lift the triangle vertices
    P3 a3 = df::lift(a2);
    P3 b3 = df::lift(b2);
    P3 c3 = df::lift(c2);

    CGAL::Segment_2<K> edge_ad_2d(a2, d2);
    CGAL::Segment_2<K> edge_bd_2d(b2, d2);
    CGAL::Segment_2<K> edge_cd_2d(c2, d2);

    // fast path: are (a,d), (b,d), (c,d) already in the lower triangulation 
    // if yes, we flip is conforming
    bool ad_in_lower = false;
    bool bd_in_lower = false;
    bool cd_in_lower = false;

    // bool which tells us if we need to check this edge or can skip it
    bool check_this_edge = true;

    // scan all edges in lower triangulation
    for (auto e = lower.finite_edges_begin(); e != lower.finite_edges_end(); ++e) {
        auto f  = e->first;
        int  ei = e->second;

        // vertex handles of edge in lower triangulation
        auto vu = f->vertex(lower.cw(ei));
        auto vv = f->vertex(lower.ccw(ei));

        // global ids of these vertices
        df::vertex_id u = vu->info();
        df::vertex_id v = vv->info();

        if ((u == ia && v == id) || (u == id && v == ia))
            ad_in_lower = true;
        if ((u == ib && v == id) || (u == id && v == ib))
            bd_in_lower = true;
        if ((u == ic && v == id) || (u == id && v == ic))
            cd_in_lower = true;

        // early exit if all three edges (a,d), (b,d), (c,d) found
        if (ad_in_lower && bd_in_lower && cd_in_lower) {
            std::cout << "[conform] all three edges (" << ia << "," << id << "), ("
                    << ib << "," << id << "), (" << ic << "," << id
                    << ") are already in lower triangulation -> insertion is conforming\n";
            return true;
        }


        P2 u2 = vu->point();
        P2 v2 = vv->point();

        CGAL::Segment_2<K> edge_uv_2d(u2, v2);

        P3 u3 = df::lift(u2);
        P3 v3 = df::lift(v2);

        CGAL::Segment_3<K> seg_uv_3d(u3, v3);
        CGAL::Segment_3<K> seg_ad_3d(a3, d3);
        CGAL::Segment_3<K> seg_bd_3d(b3, d3);
        CGAL::Segment_3<K> seg_cd_3d(c3, d3);

        check_this_edge = true; // reset for each edge

        // edge (a,d)
        if (ad_in_lower == false) {

            // skip edges that are incident to d or incident to a ("intersection" at endpoint)
            if ((u == ia || u == id || v == ia || v == id)) {
                check_this_edge = false;
            } 
            
            // and also skip the edge (b,c)
            if ((u == ib && v == ic) || (u == ic && v == ib)) {
                check_this_edge = false;
            } 
            if (check_this_edge == true){
                // cheap 2D filter: if segments don't even intersect in 2D, ignore
                if (CGAL::do_intersect(edge_ad_2d, edge_uv_2d) == true){ 
                    

                    // if they intersect in 3D, flip is not conforming
                    if (CGAL::do_intersect(seg_ad_3d, seg_uv_3d)) {
                        std::cout << "[conform] BLOCK: (a,d)=(" << ia << "," << id << ") "
                                << "intersects lower (u,v)=(" << u << "," << v << ") in 3D\n";
                        return false;
                    }

                    // if they intersect in 2D but not in 3D, we need to check heights
                    // compute height of (a,d) at intersection with (u,v)
                    int cmp_ad = compare_heights_at_intersection(a2, d2, u2, v2, a3, d3, u3, v3);
                    auto orientation_ad = oriented_height_sign(a2, d2, u2, v2, a3, d3, u3, v3);

                    if (orientation_ad == CGAL::NEGATIVE){
                        std::cout << "[conform] PASS: (a,d)=(" << ia << "," << id << ") "
                                << "above lower (u,v)=(" << u << "," << v << ")\n";
                        std::cout << "cmp_ad = " << cmp_ad << "\n";  
                    }
                    if (orientation_ad == CGAL::POSITIVE){
                        std::cout << "[conform] BLOCK: (a,d)=(" << ia << "," << id << ") "
                                << "below lower (u,v)=(" << u << "," << v << ")\n";
                        std::cout << "cmp_ad = " << cmp_ad << "\n";
                        return false;
                    }
                    

                }
            }
        }
        check_this_edge = true; // reset for next edge

        // edge (b,d)
        if (bd_in_lower == false) {
   
            // skip edges that are incident to d or incident to b
            if (u == ib || u == id || v == ib || v == id) {
                check_this_edge = false;
            } 

            // and also skip the edge (a,c)
            if ((u == ia && v == ic) || (u == ic && v == ia)){ 
                check_this_edge = false;
            }
           
            if (check_this_edge == true){
                // cheap 2D filter: if segments don't even intersect in 2D, ignore
                if (CGAL::do_intersect(edge_bd_2d, edge_uv_2d) == true){ 


                    if (CGAL::do_intersect(seg_bd_3d, seg_uv_3d)) {
                        std::cout << "[conform] BLOCK: (b,d)=(" << ib << "," << id << ") "
                                << "intersects lower (u,v)=(" << u << "," << v << ") in 3D\n";
                        return false;
                    }

                    int cmp_bd = compare_heights_at_intersection(b2, d2, u2, v2, b3, d3, u3, v3);
                    auto orientation_bd = oriented_height_sign(b2, d2, u2, v2, b3, d3, u3, v3);
                    
                    if (orientation_bd == CGAL::NEGATIVE){
                        std::cout << "[conform] PASS: (b,d)=(" << ib << "," << id << ") "
                                << "above lower (u,v)=(" << u << "," << v << ")\n";
                        std::cout << "cmp_bd = " << cmp_bd << "\n";
                    }
                
                    if (orientation_bd == CGAL::POSITIVE){
                        std::cout << "[conform] BLOCK: (b,d)=(" << ib << "," << id << ") "
                                << "below lower (u,v)=(" << u << "," << v << ")\n";
                        std::cout << "cmp_bd = " << cmp_bd << "\n";
                        return false;
                    }
                    
                }
            }
        }
        check_this_edge = true; // reset for next edge
        // edge (c,d)
        if (cd_in_lower == false) {
        
            // skip edges that are incident to d or incident to c
            if (u == ic || u == id || v == ic || v == id) {
                check_this_edge = false;
            }

            // and also skip the edge (a,b)
            if ((u == ia && v == ib) || (u == ib && v == ia)) {
                check_this_edge = false;
            }
            

            if (check_this_edge == true){
                // cheap 2D filter: if segments don't even intersect in 2D, ignore
                if (CGAL::do_intersect(edge_cd_2d, edge_uv_2d) == true){


                    if (CGAL::do_intersect(seg_cd_3d, seg_uv_3d)) {
                        std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
                                << "intersects lower (u,v)=(" << u << "," << v << ") in 3D\n";
                        return false;
                    }

                    int cmp_cd = compare_heights_at_intersection(c2, d2, u2, v2, c3, d3, u3, v3);
                    auto orientation_cd = oriented_height_sign(c2, d2, u2, v2, c3, d3, u3, v3);
                
                    
                    if (orientation_cd == CGAL::NEGATIVE){
                        std::cout << "[conform] PASS: (c,d)=(" << ic << "," << id << ") "
                                << "above lower (u,v)=(" << u << "," << v << ")\n";
                        std::cout << "cmp_cd = " << cmp_cd << "\n";
                    }
                    if (orientation_cd == CGAL::POSITIVE){
                        std::cout << "[conform] BLOCK: (c,d)=(" << ic << "," << id << ") "
                                << "below lower (u,v)=(" << u << "," << v << ")\n";
                        std::cout << "cmp_cd = " << cmp_cd << "\n";
                        return false;
                    }

                }  
            }
        }
    } // end for all edges in lower triangulation

    std::cout << "[conform] insertion of vertex with global id " << id
              << " is conforming\n";
    return true;
}
}}
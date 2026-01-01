#include "compare.h" 
#include "height_test.h"
#include "geometry_utils.h"
#include "poset.h"   

#include <CGAL/Segment_2.h>
#include <CGAL/Segment_3.h>
#include <CGAL/intersections.h>
#include <iostream>
#include <unordered_set>
#include <CGAL/Triangulation_2.h>



namespace pst2 {

using K   = df::K;
using P2  = df::P2;
using P3  = df::P3;
using Tri = df::Tri2;
using Seg2 = CGAL::Segment_2<K>;
using Seg3 = CGAL::Segment_3<K>;


static df::Tri2 build_from_history(const df::InputData& D, const pst::Node& node) {
    Tri T = D.tri_upper;              // root triangulation
    pst::replay_history_poset(T, node.history, D);
    return T;
}


static bool edge_in_triangulation(df::vertex_id a, df::vertex_id b, const Tri& T){
    for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
        auto f  = e->first;
        int  ei = e->second;

        auto va = f->vertex(T.cw(ei));
        auto vb = f->vertex(T.ccw(ei));

        df::vertex_id ia = va->info();
        df::vertex_id ib = vb->info();

        if ((ia == a && ib == b) || (ia == b && ib == a)) return true;
    }
    return false;
}

static bool all_edges_in(const df::Tri2& T1, const df::Tri2& T2){
    for (auto e = T1.finite_edges_begin(); e != T1.finite_edges_end(); ++e) {
        auto f  = e->first;
        int  i  = e->second;
        df::vertex_id a = f->vertex(T1.cw(i))->info();
        df::vertex_id b = f->vertex(T1.ccw(i))->info();

        bool found = false;
        for (auto eb = T2.finite_edges_begin(); eb != T2.finite_edges_end(); ++eb) {
            auto fb = eb->first;
            int  ib = eb->second;
            df::vertex_id c = fb->vertex(T2.cw(ib))->info();
            df::vertex_id d = fb->vertex(T2.ccw(ib))->info();
            if ((a == c && b == d) || (a == d && b == c)) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}


static bool check_edges(const df::Tri2& T1, const df::Tri2& T2){
    for (auto ea = T1.finite_edges_begin(); ea != T1.finite_edges_end(); ++ea) {
        auto fa = ea->first;
        int  ia = ea->second;

        auto va = fa->vertex(T1.cw(ia));
        auto vb = fa->vertex(T1.ccw(ia));
        df::vertex_id a = va->info();
        df::vertex_id b = vb->info();

        // skip if edge already in T2
        bool shared = edge_in_triangulation(a,b,T2);
        if (shared) continue;

        Seg2 seg1(va->point(), vb->point());

        for (auto eb = T2.finite_edges_begin(); eb != T2.finite_edges_end(); ++eb) {
            auto fb = eb->first;
            int  ib = eb->second;

            auto vc = fb->vertex(T2.cw(ib));
            auto vd = fb->vertex(T2.ccw(ib));
            df::vertex_id c = vc->info();
            df::vertex_id d = vd->info();

            if (c == a || c == b || d == a || d == b) continue;

            Seg2 seg2(vc->point(), vd->point());
            if (!CGAL::do_intersect(seg1, seg2)) continue;

            auto s = df::oriented_height_sign(
                va->point(), vb->point(),
                vc->point(), vd->point(),
                df::lift(va->point()), df::lift(vb->point()),
                df::lift(vc->point()), df::lift(vd->point())
            );

            // NEGATIVE: T1 edge ABOVE T2 edge ->  T1 not <=2 T2
            if (s == CGAL::NEGATIVE) {
                return false;
            }
        }
    }
    return true;
}


static bool check_extra_vertices_T2_vs_T1(
    const Tri& T1,
    const Tri& T2)
{
    // collect vertex ids of T1
    std::unordered_set<df::vertex_id> ids1;
    ids1.reserve((size_t)T1.number_of_vertices());

    for (auto v = T1.finite_vertices_begin(); v != T1.finite_vertices_end(); ++v) {
        ids1.insert(v->info());
    }

    // loop over vertices of T2 not in T1
    for (auto v2 = T2.finite_vertices_begin(); v2 != T2.finite_vertices_end(); ++v2) {

        df::vertex_id id = v2->info();
        if (ids1.find(id) != ids1.end()) continue; // shared vertex

        const P2& p = v2->point();

        // locate p in T1
        auto fh = T1.locate(p);
        if (fh == nullptr) {
            // should not happen
            return false;
        }

        // face vertices in T1
        const P2& a = fh->vertex(0)->point();
        const P2& b = fh->vertex(1)->point();
        const P2& c = fh->vertex(2)->point();

        const P3& pa = df::lift(a);
        const P3& pb = df::lift(b);
        const P3& pc = df::lift(c);
        const P3& pp = df::lift(p);

        auto s = df::oriented_height_sign(a, b, c, p, pa, pb, pc, pp);

        // NEGATIVE means: lifted p is BELOW the face -> violation
        if (s == CGAL::NEGATIVE) {
            return false;
        }
    }

    return true;
}



// return true iff t1 <=_2 t2  (T1 surface never above T2)
bool compare(int t1, int t2, const df::InputData& D, const std::vector<pst::Node>& nodes)
{
    Tri T1 = build_from_history(D, nodes[t1]);
    Tri T2 = build_from_history(D, nodes[t2]);

    // check extra vertices of T2 against T1: locate vertex of T2 in T1, check height
    if (!check_extra_vertices_T2_vs_T1(T1, T2)) {
        return false;
    }

    // check all edges of T1 against T2
    if (all_edges_in(T1, T2)) {
        // if (T2 <=2 T1) then reject (T1 <=2 T2)
        return !check_edges(T2, T1);
    }

    // check all edges of T1 against T2: height/ intersection test
    return check_edges(T1, T2);
}



} // namespace pst2

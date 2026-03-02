#include "poset_ear_star.h"

#include <CGAL/enum.h>
#include <algorithm>
#include <iostream>

namespace pst_es {

// ------------------------------------------------------------

static inline std::array<df::vertex_id,2>
norm_edge(df::vertex_id a, df::vertex_id b) {
    if (a > b) std::swap(a,b);
    return {a,b};
}

// ------------------------------------------------------------
// signature
// ------------------------------------------------------------

static TriSignature make_signature(const df::Tri2& tri) {
    TriSignature sig;

    for (auto e = tri.finite_edges_begin();
         e != tri.finite_edges_end(); ++e) {

        auto f = e->first;
        int  i = e->second;

        auto va = f->vertex(tri.cw(i));
        auto vb = f->vertex(tri.ccw(i));

        sig.edges.push_back(
            norm_edge(va->info(), vb->info())
        );
    }

    std::sort(sig.edges.begin(), sig.edges.end());
    sig.edges.erase(
        std::unique(sig.edges.begin(), sig.edges.end()),
        sig.edges.end()
    );

    return sig;
}

std::size_t TriSignatureHash::operator()(
    TriSignature const& s) const noexcept {

    std::size_t h = 0;
    for (auto const& e : s.edges) {
        h ^= std::hash<std::size_t>{}(e[0]) + 0x9e3779b9 + (h<<6) + (h>>2);
        h ^= std::hash<std::size_t>{}(e[1]) + 0x9e3779b9 + (h<<6) + (h>>2);
    }
    return h;
}

// ------------------------------------------------------------
// convex quad test
// ------------------------------------------------------------

static bool quad_convex(
    const df::P2& A,
    const df::P2& B,
    const df::P2& C,
    const df::P2& D)
{
    auto o1 = CGAL::orientation(A,B,C);
    auto o2 = CGAL::orientation(A,B,D);
    if (o1 == CGAL::COLLINEAR || o2 == CGAL::COLLINEAR) return false;
    if (o1 == o2) return false;

    auto o3 = CGAL::orientation(C,D,A);
    auto o4 = CGAL::orientation(C,D,B);
    if (o3 == CGAL::COLLINEAR || o4 == CGAL::COLLINEAR) return false;
    if (o3 == o4) return false;

    return true;
}

// ------------------------------------------------------------
// find flippable edges
// ------------------------------------------------------------

static std::vector<std::array<df::vertex_id,2>>
find_flips(const df::Tri2& tri)
{
    std::vector<std::array<df::vertex_id,2>> out;

    for (auto e = tri.finite_edges_begin();
         e != tri.finite_edges_end(); ++e) {

        auto f = e->first;
        int  i = e->second;
        auto g = f->neighbor(i);

        if (tri.is_infinite(g)) continue;

        auto va = f->vertex(tri.cw(i));
        auto vb = f->vertex(tri.ccw(i));
        auto vc = f->vertex(i);
        int j = tri.mirror_index(f,i);
        auto vd = g->vertex(j);

        if (!quad_convex(
                va->point(),
                vb->point(),
                vc->point(),
                vd->point())) continue;

        out.push_back(norm_edge(
            va->info(), vb->info()));
    }

    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());

    return out;
}

// ------------------------------------------------------------
// apply flip
// ------------------------------------------------------------

static bool apply_flip(
    df::Tri2& tri,
    df::vertex_id a,
    df::vertex_id b)
{
    for (auto e = tri.finite_edges_begin();
         e != tri.finite_edges_end(); ++e) {

        auto f = e->first;
        int  i = e->second;

        auto va = f->vertex(tri.cw(i));
        auto vb = f->vertex(tri.ccw(i));

        if ((va->info()==a && vb->info()==b) ||
            (va->info()==b && vb->info()==a))
        {
            if (tri.is_infinite(f->neighbor(i))) return false;
            tri.flip(f,i);
            return true;
        }
    }
    return false;
}

// ------------------------------------------------------------
// build full poset
// ------------------------------------------------------------

FlipPoset build_flip_poset(
    const df::Tri2& tri_start,
    df::vertex_id star_id,
    const df::P2& star_point)
{
    FlipPoset P;

    Node root;
    root.sig = make_signature(tri_start);
    P.nodes.push_back(root);
    P.sig_to_idx[root.sig] = 0;

    for (int idx=0; idx<(int)P.nodes.size(); ++idx) {

        df::Tri2 tri = tri_start;
        reconstruct_triangulation(
            tri_start, P, idx,
            star_id, star_point, tri);

        // flips
        for (auto const& e : find_flips(tri)) {

            df::Tri2 tri2 = tri;
            apply_flip(tri2, e[0], e[1]);

            TriSignature sig2 = make_signature(tri2);
            if (P.sig_to_idx.count(sig2)) continue;

            int child = P.nodes.size();
            Node n;
            n.sig = sig2;
            n.parent = idx;
            n.move_from_parent = {MoveKind::Flip, e[0], e[1]};

            P.nodes.push_back(n);
            P.sig_to_idx[sig2] = child;
            P.nodes[idx].children.push_back(child);
        }

        // insertion of star
        bool star_present = false;
        for (auto v = tri.finite_vertices_begin();
             v != tri.finite_vertices_end(); ++v)
            if (v->info()==star_id)
                star_present = true;

        if (!star_present) {
            df::Tri2 tri2 = tri;
            auto fh = tri2.locate(star_point);
            auto vh = tri2.insert(star_point, fh);
            vh->info() = star_id;

            TriSignature sig2 = make_signature(tri2);
            if (!P.sig_to_idx.count(sig2)) {
                int child = P.nodes.size();
                Node n;
                n.sig = sig2;
                n.parent = idx;
                n.move_from_parent = {MoveKind::InsertStar,0,0};
                P.nodes.push_back(n);
                P.sig_to_idx[sig2] = child;
                P.nodes[idx].children.push_back(child);
            }
        }
        // deletion of star (only if star is present and has degree 3)
        {
            bool star_present_deg3 = false;
            for (auto v = tri.finite_vertices_begin(); v != tri.finite_vertices_end(); ++v) {
                if (v->info() == star_id && tri.degree(v) == 3) {
                    star_present_deg3 = true;
                    break;
                }
            }

            if (star_present_deg3) {
                df::Tri2 tri2 = tri;

                bool found = false;
                df::Tri2::Vertex_handle vh_star; // do NOT set to nullptr; it's not a pointer in CGAL

                for (auto vit = tri2.finite_vertices_begin(); vit != tri2.finite_vertices_end(); ++vit) {
                    if (vit->info() == star_id) {
                        vh_star = vit;
                        found = true;
                        break;
                    }
                }

                if (found && tri2.degree(vh_star) == 3) {
                    tri2.remove_degree_3(vh_star);

                    TriSignature sig2 = make_signature(tri2);
                    if (!P.sig_to_idx.count(sig2)) {
                        int child = (int)P.nodes.size();
                        Node n;
                        n.sig = std::move(sig2);
                        n.parent = idx;
                        n.move_from_parent = {MoveKind::DeleteStar, 0, 0};
                        P.nodes.push_back(std::move(n));
                        P.sig_to_idx[P.nodes.back().sig] = child;
                        P.nodes[idx].children.push_back(child);
                    }
                }
            }
        }
    }

    std::cout << "[ear_star] nodes: "
              << P.nodes.size() << "\n";

    return P;
}

// ------------------------------------------------------------
// reconstruction
// ------------------------------------------------------------

bool reconstruct_triangulation(
    const df::Tri2& tri_start,
    const FlipPoset& P,
    int idx,
    df::vertex_id star_id,
    const df::P2& star_point,
    df::Tri2& out_tri)
{
    out_tri = tri_start;

    std::vector<Move> chain;
    for (int cur=idx;
         P.nodes[cur].parent!=-1;
         cur=P.nodes[cur].parent)
        chain.push_back(P.nodes[cur].move_from_parent);

    std::reverse(chain.begin(), chain.end());

    for (auto const& m : chain) {

        if (m.kind==MoveKind::Flip)
            apply_flip(out_tri, m.a, m.b);

        else if (m.kind==MoveKind::InsertStar) {
            auto fh = out_tri.locate(star_point);
            auto vh = out_tri.insert(star_point, fh);
            vh->info() = star_id;
        }

        else if (m.kind==MoveKind::DeleteStar) {
            for (auto v = out_tri.finite_vertices_begin();
                 v != out_tri.finite_vertices_end(); ++v)
                if (v->info()==star_id) {
                    out_tri.remove_degree_3(v);
                    break;
                }
        }
    }

    return true;
}

} // namespace pst_es
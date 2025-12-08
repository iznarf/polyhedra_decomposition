#include "debug.h"
#include "input.h"
#include "visualization.h"
#include "geometry_utils.h"
#include "check_edges.h"
#include "insertion.h"
#include <iostream>
#include <set>
#include <vector>
#include <CGAL/Segment_2.h>
#include <CGAL/intersections.h>

namespace df {
   

// print mapping from polyscope local indices to global vertex ids in current triangulation
static std::vector<df::vertex_id>
collect_vertex_ids_in_order(const df::Tri2& tri) {
    std::vector<df::vertex_id> ids;
    ids.reserve(tri.number_of_vertices());
    for (auto v = tri.finite_vertices_begin(); v != tri.finite_vertices_end(); ++v) {
        ids.push_back(v->info());
    }
    return ids;
}

// prints all edges in the current triangulation
void debug_print_edge_list(const InputData& D) {
    const auto& tri = D.tri_current;

    std::vector<std::pair<vertex_id, vertex_id>> edges;

    for (auto e = tri.finite_edges_begin(); e != tri.finite_edges_end(); ++e) {
        auto f  = e->first;
        int ei  = e->second;

        auto va = f->vertex(tri.cw(ei));
        auto vb = f->vertex(tri.ccw(ei));

        vertex_id a = va->info();
        vertex_id b = vb->info();

        // canonicalize (min,max)
        if (b < a) std::swap(a, b);

        edges.emplace_back(a, b);
    }

    std::cout << "[debug] edges in current triangulation (global ids):\n";

    for (size_t i = 0; i < edges.size(); i++) {
        std::cout << "(" << edges[i].first << "," << edges[i].second << ")";
        if (i + 1 < edges.size()) std::cout << ", ";
    }
    std::cout << "\n";
}


    // returns all undirected edges of a triangulation as normalized (min,max) pairs
    static std::vector<std::pair<vertex_id, vertex_id>> collect_edges(const Tri2& T) {
        std::vector<std::pair<vertex_id, vertex_id>> edges;

        for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
            auto f = e->first;
            int i  = e->second;

            auto va = f->vertex(T.cw(i));
            auto vb = f->vertex(T.ccw(i));

            vertex_id a = va->info();
            vertex_id b = vb->info();
            if (a > b) std::swap(a, b);

            edges.emplace_back(a, b);
        }

        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        return edges;
    }


    bool triangulations_equal(const Tri2& A, const Tri2& B) {
        auto EA = collect_edges(A);
        auto EB = collect_edges(B);
        return EA == EB;
    }

void debug_print_local_to_global_map(
    const df::InputData& D,
    TriKind which,
    const std::vector<int>& local_indices)
{
    const df::Tri2& tri = (which == TriKind::Current)
                        ? D.tri_current
                        : D.tri_lower;

    // get global ids in polyscope order for this triangulation
    auto ids = collect_vertex_ids_in_order(tri);

    const char* label = (which == TriKind::Current) ? "current" : "lower";

    std::cout << "\n[debug] " << label
              << " polyscope local -> global mapping:\n";

    // if no subset given: dump full map
    if (local_indices.empty()) {
        for (int li = 0; li < static_cast<int>(ids.size()); ++li) {
            df::vertex_id gid = ids[li];
            std::cout << "  local " << li << " -> global " << gid << "\n";
        }
        return;
    }

    // otherwise: only print the requested local indices
    for (int li : local_indices) {
        if (li < 0 || li >= static_cast<int>(ids.size())) {
            std::cout << "  local " << li << " : out of range (0.."
                      << (ids.size() - 1) << ")\n";
            continue;
        }
        df::vertex_id gid = ids[li];
        std::cout << "  local " << li << " -> global " << gid << "\n";
    }
}

void print_step_history(const df::InputData& D) {
    std::cout << "\n==== Step history (" 
              << D.step_history.size() << " steps) ====\n";

    int i = 0;
    for (const auto& s : D.step_history) {
        std::cout << "Step " << i++ << ": ";

        if (s.kind == df::StepKind::EdgeFlip_down) {
            std::cout << "EdgeFlip  (a=" << s.a 
                      << ", b=" << s.b
                      << ", c=" << s.c
                      << ", d=" << s.d << ")";
        } else {
            std::cout << "VertexInsertion_down  (face=(" 
                      << s.a << "," << s.b << "," << s.c 
                      << "), new=" << s.d << ")";
        }

        std::cout << "\n";
    }

    std::cout << "========================================\n\n";
}

// helper: given an edge (ia, ib), find the 4 vertices of the quad around it
// in current triangulation (a,b,c,d) with edge (a,b) and opposite vertices (c,d)
static bool
compute_quad_for_edge(const Tri2& tri,
                      vertex_id ia, vertex_id ib,
                      std::array<vertex_id,4>& out)
{
    // find the edge in the triangulation
    Tri2::Face_handle f;
    int edge_idx = -1;
    bool found = false;

    for (auto e = tri.finite_edges_begin(); e != tri.finite_edges_end(); ++e) {
        auto fh = e->first;
        int  ei = e->second;

        auto va = fh->vertex(tri.cw(ei));
        auto vb = fh->vertex(tri.ccw(ei));

        vertex_id a = va->info();
        vertex_id b = vb->info();

        if ((a == ia && b == ib) || (a == ib && b == ia)) {
            f        = fh;
            edge_idx = ei;
            found    = true;
            break;
        }
    }

    if (!found) {
        std::cout << "[debug] compute_quad_for_edge: edge (" << ia << "," << ib
                  << ") not found in tri_current\n";
        return false;
    }

    if (tri.is_infinite(f)) {
        std::cout << "[debug] compute_quad_for_edge: incident face is infinite\n";
        return false;
    }

    auto g = f->neighbor(edge_idx);
    if (tri.is_infinite(g)) {
        std::cout << "[debug] compute_quad_for_edge: neighbor face is infinite\n";
        return false;
    }

    auto va = f->vertex(tri.cw(edge_idx));
    auto vb = f->vertex(tri.ccw(edge_idx));
    auto vc = f->vertex(edge_idx);
    int  mi = tri.mirror_index(f, edge_idx);
    auto vd = g->vertex(mi);

    out = { va->info(), vb->info(), vc->info(), vd->info() };
    return true;
}


// this function collects the tetrahedra that correspond to edge flips and vertex insertions if the algorithm did not terminate correctly
std::vector<DebugTetrahedron>
collect_debug_tetrahedra(const InputData& D)
{
    std::vector<DebugTetrahedron> result;

    const Tri2& current = D.tri_current;
    const Tri2& target  = D.tri_lower;

    // collecting locally non-regular edges for 2-2 flip tetrahedra 
    auto non_regular_edges = df::reg::find_locally_non_regular_edges(current);

    for (const auto& e : non_regular_edges) {
        vertex_id ia = e[0];
        vertex_id ib = e[1];

        // we should only collect downflippable edges that are not in the target triangulation 
        // right now we collect all non regular edges 

        bool in_target = df::reg::edge_in_target(ia, ib, D);
       

        std::array<vertex_id,4> quad_ids;
        if (!compute_quad_for_edge(current, ia, ib, quad_ids)) {
            // if we can't build the quad, just skip this one
            continue;
        }

        if(in_target == true){
            // edge is already in target triangulation, skip it
            continue;
        }

            DebugTetrahedron tet;
            tet.verts = quad_ids;
            tet.kind  = DebugTetKind::EdgeFlip;
            result.push_back(tet);
        }
    

    // collect missing vertices for insertion tetrahedra 
    auto missing = df::find_missing_vertices(current, target);

    for (vertex_id id : missing) {
        const P2& p = D.points2d[id];

        Tri2::Locate_type lt;
        int li;
        Tri2::Face_handle fh = current.locate(p, lt, li);

        if (lt != Tri2::FACE || current.is_infinite(fh)) {
            std::cout << "[debug] collect_debug_tetrahedra: "
                      << "missing vertex " << id
                      << " not in a finite face (lt=" << lt << ")\n";
            continue;
        }

        vertex_id a = fh->vertex(0)->info();
        vertex_id b = fh->vertex(1)->info();
        vertex_id c = fh->vertex(2)->info();

        DebugTetrahedron tet;
        tet.verts = { a, b, c, id };
        tet.kind  = DebugTetKind::Insertion;
        result.push_back(tet);
        // tetrahedron with vertices 0,1,2,3 has faces 
        // {0,1,2}, {0,3,1}, {1,3,2}, {0,2,3}

    }

    return result;
}
}
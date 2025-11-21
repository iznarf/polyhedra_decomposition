#include "debug.h"
#include "input.h"
#include "visualization.h"
#include "geometry_utils.h"
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




 // helper function to collect edges in current but not in target
std::vector<std::pair<vertex_id, vertex_id>>
edges_in_current_not_in_target(const Tri2& current, const Tri2& target)
{
    auto canon = [](vertex_id a, vertex_id b) {
        if (b < a) std::swap(a, b);
        return std::make_pair(a, b);
    };

    // Store all edges of target triangulation
    std::set<std::pair<vertex_id, vertex_id>> target_edges;

    for (auto e = target.finite_edges_begin(); e != target.finite_edges_end(); ++e) {
        auto f  = e->first;
        int ei  = e->second;

        vertex_id a = f->vertex(target.cw(ei))->info();
        vertex_id b = f->vertex(target.ccw(ei))->info();

        target_edges.insert(canon(a, b));
    }

    std::vector<std::pair<vertex_id, vertex_id>> diff;
        for (auto e = current.finite_edges_begin(); e != current.finite_edges_end(); ++e) {
            auto f  = e->first;
            int ei  = e->second;

            vertex_id a = f->vertex(current.cw(ei))->info();
            vertex_id b = f->vertex(current.ccw(ei))->info();

            auto ce = canon(a, b);
            if (target_edges.find(ce) == target_edges.end()) { 
                diff.push_back(ce);
            }
        }

        return diff;
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

// checks whether there are edges in current triangulation that are not in lower triangulation
// if returns true, there are differences
bool edge_diff_with_lower(const df::InputData& D) {
    auto diff = edges_in_current_not_in_target(D.tri_current, D.tri_lower);
    auto ids      = viz::present_ids(D.tri_current);
    auto to_local = viz::make_local_index(ids);
    

    std::cout << "\n=== Edges in current but NOT in lower (global + local ids in current) ===\n";
    if (diff.empty()) {
        std::cout << "(none)\n";
        return false;
    }

    for (auto [a, b] : diff) {
        int la = to_local.at(a);
        int lb = to_local.at(b);
        std::cout << "  global (" << a << "," << b << ")"
                  << "  local (" << la << "," << lb << ")\n";
    }
    // if difference was found, return true 
    return true;
}

// check if an edge (i,j) exists in a triangulation T
static bool edge_in_triangulation(const Tri2& T,
                                  vertex_id i,
                                  vertex_id j)
{
    if (i == j) return false;
    if (i > j) std::swap(i, j);

    for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
        auto  f  = e->first;
        int   ei = e->second;

        auto  va = f->vertex(T.cw(ei));
        auto  vb = f->vertex(T.ccw(ei));

        vertex_id a = va->info();
        vertex_id b = vb->info();
        if (a > b) std::swap(a, b);

        if (a == i && b == j) {
            return true;
        }
    }
    return false;
}

// print all edges created by insertion of vertex id
void debug_edges_created_by_insertion(vertex_id id, const InputData& D)
{
    const Tri2& cur   = D.tri_current;
    const Tri2& lower = D.tri_lower;

    std::cout << "\n[debug-insert] edges incident to inserted vertex " << id << ":\n";

    // 1) find the vertex handle for the newly inserted vertex
    Tri2::Vertex_handle vd = Tri2::Vertex_handle();
    for (auto v = cur.finite_vertices_begin(); v != cur.finite_vertices_end(); ++v) {
        if (v->info() == id) {
            vd = v;
            break;
        }
    }

    if (vd == Tri2::Vertex_handle()) {
        std::cout << "[debug-insert] ERROR: could not find vertex " << id
                  << " in current triangulation\n";
        return;
    }

    // 2) iterate over all incident vertices -> edges (id, n)
    Tri2::Vertex_circulator vc = cur.incident_vertices(vd);
    Tri2::Vertex_circulator done = vc;

    
    if (vc == 0) {
        std::cout << "[debug-insert] vertex has no incident vertices (unexpected)\n";
        return;
    }

    do {
        if (cur.is_infinite(vc)) {
            ++vc;
            continue;
        }

        vertex_id n = vc->info();
        vertex_id a = id;
        vertex_id b = n;
        if (a > b) std::swap(a, b);

        bool in_lower = edge_in_triangulation(lower, a, b);

        std::cout << "  edge (" << a << "," << b << ")"
                  << (in_lower ? "  [in lower]\n" : "  [NOT in lower]\n");

        ++vc;
    } while (vc != done);
}

// check edge (i,j) against lower triangulation for intersections
void debug_check_edge_against_lower(df::vertex_id i, df::vertex_id j, const df::InputData& D) {
    const auto& pts = D.points2d;
    CGAL::Segment_2<K> e2d(pts[i], pts[j]);

    std::cout << "[post-insert-debug] checking edge (" << i << "," << j << ") against lower triangulation\n";

    for (auto e = D.tri_lower.finite_edges_begin();
         e != D.tri_lower.finite_edges_end(); ++e) {

        auto f  = e->first;
        int  ei = e->second;

        auto vu = f->vertex(D.tri_lower.cw(ei));
        auto vv = f->vertex(D.tri_lower.ccw(ei));

        df::vertex_id u = vu->info();
        df::vertex_id v = vv->info();

        CGAL::Segment_2<K> uv2d(vu->point(), vv->point());

        if (CGAL::do_intersect(e2d, uv2d)) {
            std::cout << "  --> intersects lower edge (" << u << "," << v << ")\n";
        }
    }
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

        if (s.kind == df::StepKind::EdgeFlip) {
            std::cout << "EdgeFlip  (a=" << s.a 
                      << ", b=" << s.b
                      << ", c=" << s.c
                      << ", d=" << s.d << ")";
        } else {
            std::cout << "VertexInsertion  (face=(" 
                      << s.a << "," << s.b << "," << s.c 
                      << "), new=" << s.d << ")";
        }

        std::cout << "\n";
    }

    std::cout << "========================================\n\n";
}









}
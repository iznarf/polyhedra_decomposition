#pragma once
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <random>

namespace df {

    //using K  = CGAL::Simple_cartesian<double>; // cartesian kernel
    using K  = CGAL::Exact_predicates_inexact_constructions_kernel;
    using P2 = K::Point_2;
    using P3 = K::Point_3;
    using vertex_id = std::size_t;  // global vertex id
    using vertex_base = CGAL::Triangulation_vertex_base_with_info_2<vertex_id, K>;
    using face_base  = CGAL::Triangulation_face_base_2<K>;
    using Tri_ds = CGAL::Triangulation_data_structure_2<vertex_base, face_base>;
    using Tri2 = CGAL::Triangulation_2<K, Tri_ds>;



    inline constexpr double PI = 3.14159265358979323846;

    enum class StepKind { EdgeFlip, VertexInsertion };

    struct StepRecord {
    StepKind kind;

    // For both flip + insertion, we can store 4 ids:
    // Edge flip: (a,b) edge, (c,d) opposite vertices before flip
    // Insertion: (a,b,c) is the face containing d, and d is the new vertex
    vertex_id a, b, c, d;
    };


    struct FlipRecord {
        vertex_id a, b;  // edge BEFORE flip  (a,b)
        vertex_id c, d;  // edge AFTER  flip (c,d)

        FlipRecord() = default;
        FlipRecord(vertex_id aa, vertex_id bb, vertex_id cc, vertex_id dd)
            : a(aa), b(bb), c(cc), d(dd) {}
    };

    struct InputData {
        std::vector<P2> points2d;   // all random points of planar point set A
        Tri2            tri_upper;  // source: triangulation using hull only (no interior vertices)
        Tri2            tri_lower;  // target: triangulation using all points (includes interior)
        Tri2           tri_current;  // triangulation of the current state in the algorithm
        Tri2          tri_replay;  // triangulation for applying the recorded flips

        std::vector<df::FlipRecord> flip_history; // edge flip history 
        std::vector<std::pair<vertex_id, vertex_id>> hull_edges;
        std::vector<StepRecord> step_history;  // sequence of all flips: edge flips and insertions


  
    };

    // check if edge (u,v) is a hull edge
    bool is_hull_edge(const InputData& D, vertex_id u, vertex_id v);

    
    // create random input polyhedron (not necessarily valid)
    InputData make_random_input(int n_points, unsigned seed);

    // create random VALID input polyhedron
    InputData make_random_valid_input(int n_points, unsigned seed_start);


} // namespace df


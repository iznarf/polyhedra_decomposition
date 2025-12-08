#pragma once
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>



#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Weighted_point_2.h>


#include <CGAL/Triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>
#include <random>

namespace df {

    // simple cartesian kernel with double would be better 
    //using K  = CGAL::Simple_cartesian<double>; // cartesian kernel
    using K  = CGAL::Exact_predicates_inexact_constructions_kernel;
    using P2 = K::Point_2;
    using P3 = K::Point_3;
    using vertex_id = std::size_t;  // global vertex id
    using vertex_base = CGAL::Triangulation_vertex_base_with_info_2<vertex_id, K>;
    using face_base  = CGAL::Triangulation_face_base_2<K>;
    using Tri_ds = CGAL::Triangulation_data_structure_2<vertex_base, face_base>;
    using Tri2 = CGAL::Triangulation_2<K, Tri_ds>;
    using Polyhedron_3 = CGAL::Polyhedron_3<K>;

    // regular triangulation types: new TDS with info
    using RT_Vb_base  = CGAL::Regular_triangulation_vertex_base_2<K>;
    using RT_Vb       = CGAL::Triangulation_vertex_base_with_info_2<vertex_id, K, RT_Vb_base>;
    using RT_Fb       = CGAL::Regular_triangulation_face_base_2<K>;
    using RT_Tds      = CGAL::Triangulation_data_structure_2<RT_Vb, RT_Fb>;
    using Tri2Regular = CGAL::Regular_triangulation_2<K, RT_Tds>;
    using P2_weighted = CGAL::Weighted_point_2<K>;


    inline constexpr double PI = 3.14159265358979323846;

    enum class StepKind { EdgeFlip, VertexInsertion, VertexDeletion };

    struct StepRecord {
        StepKind kind;

        // For both flip + insertion, we can store 4 ids:
        // Edge flip: (a,b) edge, (c,d) opposite vertices before flip
        // Insertion: (a,b,c) is the face containing d, and d is the new vertex
        vertex_id a, b, c, d;
    };


    struct InputData {
        std::vector<P2> points2d;   // all random points of planar point set A
        std::vector<P2_weighted> points2d_weighted; // weighted points for regular triangulation of planar point set A
        Tri2            tri_upper;  // source: triangulation using hull only (no interior vertices)
        Tri2            tri_lower;  // target: triangulation using all points 
        Tri2           tri_current;  // triangulation of the current state in the algorithm
        Tri2          tri_replay;  // triangulation for applying the recorded flips
        Tri2            tri_poset;  // triangulation for poset computation
        Tri2Regular tri_regular; // regular triangulation of point set A

        std::vector<std::pair<vertex_id, vertex_id>> hull_edges; // this is just for checking if input is valid
        std::vector<StepRecord> step_history;  // sequence of all flips: edge flips and insertions
    };

    // check if edge (u,v) is a hull edge
    bool is_hull_edge(const InputData& D, vertex_id u, vertex_id v);

    
    // create random input polyhedron (not necessarily valid)
    InputData make_random_input(int n_points, unsigned seed);

    // create random VALID input polyhedron
    InputData make_random_valid_input(int n_points, unsigned seed_start);

    // check if all lifted points are vertices of their convex hull
    bool vertices_in_convex_position(const InputData& D);


} // namespace df


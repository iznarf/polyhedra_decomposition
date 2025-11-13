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

    struct InputData {
        std::vector<P2> points2d;   // all random points of planar point set A
        Tri2            tri_upper;  // source: triangulation using hull only (no interior vertices)
        Tri2            tri_lower;  // target: triangulation using all points (includes interior)
        Tri2           tri_current;  // triangulation of the current state in the algorithm
  
    };
  
  // create random input points, hull triangulation, full triangulation
  InputData make_random_input(int n_points = 10, unsigned seed = 42);

} // namespace df


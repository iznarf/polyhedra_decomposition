#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <random>

namespace df {

  using K  = CGAL::Exact_predicates_inexact_constructions_kernel; // change to cartesian kernel
  using P2 = K::Point_2;
  using P3 = K::Point_3;
  using Tri2 = CGAL::Triangulation_2<K>;

  inline constexpr double PI = 3.14159265358979323846;

  struct OmegaQuad {
    // omega(x,y) = a x^2 + b x y + c y^2 + d x + e y + f
    double a, b, c, d, e, f;
    double operator()(const P2& p) const {
      const double x = p.x(), y = p.y();
      return a*x*x + b*x*y + c*y*y + d*x + e*y + f;
    }
  };


    struct InputData {
    std::vector<P2> points2d;   // all random points of planar point set A
    std::vector<P2> hull2d;     // convex hull of A (CCW)
    OmegaQuad       omega;      // height
    Tri2            tri_lower;  // source: triangulation using hull only (no interior vertices)
    Tri2            tri_upper;  // target: triangulation using all points (includes interior)
    
    std::vector<P3> lifted() const; // (you can keep this if you use it elsewhere)
    };

  // Create random input: points, hull triangulation, full triangulation
  InputData make_random_input(int n_points = 10, unsigned seed = 42);

} // namespace df


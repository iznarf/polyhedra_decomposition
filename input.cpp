#include "input.h"
#include <CGAL/convex_hull_2.h>
#include <algorithm>
#include <random>
#include <cmath>
#include <cassert>

namespace df {

// randomly generate n points inside a disk of radius R, uniformly over the area (not just uniformly in radius) 
// we pick random angle and random radius for polar coordinates and then convert to cartesian coordinates
// we use uniform area to spread points evenly
// return 2D points (P2 for CGAL)
static std::vector<P2> sample_points_in_disk(int n, double R, std::mt19937& rng) {
  std::uniform_real_distribution<double> Uang(0.0, 2.0 * PI); // random angle
  std::uniform_real_distribution<double> Urad(0.0, 1.0); // random radius
  std::vector<P2> pts; pts.reserve(n);
  for (int i = 0; i < n; ++i) {
    double ang = Uang(rng); // random angle between 0 and 2pi
    double r   = R * std::sqrt(Urad(rng));     // sqrt of random radius
    pts.emplace_back(r * std::cos(ang), r * std::sin(ang)); // convert to cartesian
  }
  return pts;
}

// return lifted 3D points according to height function omega (P3 for CGAL)
std::vector<P3> InputData::lifted() const {
  std::vector<P3> out; out.reserve(points2d.size());
  for (auto& p : points2d) out.emplace_back(p.x(), p.y(), omega(p));
  return out;
}

// create random input: points, hull triangulation, full triangulation
// n_points: number of vertices 
// seed: starting value for the random number generator -> if seed stays the same, we can get the same random points
InputData make_random_input(int n_points, unsigned seed) {
  assert(n_points >= 3);
  std::mt19937 rng(seed);

  InputData D;

  // 1) random points with interior points
  D.points2d = sample_points_in_disk(n_points, /*R=*/1.0, rng);

  // 2) convex hull (CCW) of all points
  std::vector<P2> hull;
  hull.reserve(n_points);
  CGAL::convex_hull_2(D.points2d.begin(), D.points2d.end(), std::back_inserter(hull));
  D.hull2d = hull;

  // 3) triangulations
  // source triangulation: only hull vertices (no interior vertices)
  D.tri_upper.clear();
  D.tri_upper.insert(hull.begin(), hull.end());

  // target triangulation: all points (includes interior vertices)
  // (shuffle insertion order to get a "random" triangulation with Triangulation_2)
  std::vector<P2> pts_shuf = D.points2d;
  std::shuffle(pts_shuf.begin(), pts_shuf.end(), rng);
  D.tri_lower.clear();
  D.tri_lower.insert(pts_shuf.begin(), pts_shuf.end());

  // 4) convex height function omega(x,y) = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f
  D.omega = OmegaQuad{/*a=*/1.0, /*b=*/0.0, /*c=*/1.0,
                      /*d=*/0.0, /*e=*/0.0, /*f=*/0.0};

  return D;
}

} // namespace df



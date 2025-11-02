#include <polyscope/polyscope.h>
#include <CGAL/Simple_cartesian.h>
#include <portable-file-dialogs.h>

typedef CGAL::Simple_cartesian<double> Kernel;

#include "input.h"             // generate input
#include "visualization.h"     // Polyscope visualization
#include "callback.h"          // callback interface
#include "regularity_check.h"  // regularity checks for edges

// NEW: target-edge set + candidate filtering
#include "target_edges.h"
#include "candidates.h"

int main() {
  // build input polyhedron with 22 vertices
  df::InputData in = df::make_random_input(22, 50);

  // start from upper triangultion (source triangulation) which is the convex hull of input points
  df::Tri2& curr = in.tri_upper;

  // target triangulation is the lower one (all points)
  const df::Tri2& Tv = in.tri_lower;

  // initialize Polyscope
  viz::init();

  // polyhedron is union of lifted triangulations at height omega
  // visualize polyhedron
  viz::register_lifted_triangulation(in.tri_lower, in.omega, "lower");
  viz::register_lifted_triangulation(in.tri_upper, in.omega, "upper");

  // show planar triangulations of upper and lower boundary 
  viz::register_planar_triangulation(in.tri_upper, "planar, upper boundary");
  viz::register_planar_triangulation(in.tri_lower, "planar, lower boundary");

  
  //df::reg::print_triangulation_summary(curr);
  // print all internal edges for local regularity
  df::reg::print_flippable_edges_table_ids(curr, in.omega);

  // 
  const auto target_edges = df::target::build_edge_set_from_triangulation(Tv);

  // --- Filter: keep only DOWN, non-degenerate, and NOT already in Tv
  auto cand = df::cand::down_not_in_target(curr, in.omega, target_edges);

  // Print concise shortlist for cross-checking in Polyscope
  df::cand::print_candidates(cand);

  // NOTE: We are NOT flipping yet. After you confirm the shortlist matches
  // what you expect visually, we’ll call the flip and then recompute/print again.

  viz::show();
}


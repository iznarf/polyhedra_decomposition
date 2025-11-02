#pragma once
#include <vector>
#include <iosfwd>
#include "input.h"

namespace df {
namespace reg {

// Record for one flippable internal edge 
struct EdgeRecord {
  size_t ia = 0, ib = 0;     // Polyscope-style vertex indices (planar mesh)
  P2 a2, b2;                 // 2D coordinates of endpoints (handy for debug)
  bool nonreg_up   = false;  // d' above plane a'b'c'
  bool nonreg_down = false;  // d' below plane a'b'c'
  bool degenerate  = false;  // coplanar in lifted test
};

// Collect flippable internal edges (strictly convex quad around the edge),
// annotate local regularity flags under omega, return with vertex IDs
std::vector<EdgeRecord>
collect_flippable_edges_with_ids(const Tri2& tri, const OmegaQuad& omega);

// Print table of those edges (IDs match the planar mesh
void print_flippable_edges_table_ids(const Tri2& tri, const OmegaQuad& omega,
                                     std::ostream& out = std::cout);

// (Optional) Quick sanity check: print counts and the ID -> (x,y) map for the
// triangulation you’re inspecting (IDs match PLANAR Polyscope mesh).
void print_triangulation_summary(const Tri2& tri);


} // namespace reg
} // namespace df

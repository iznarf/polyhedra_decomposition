#include "candidates.h"
#include <iomanip>
#include <iostream>

namespace df { namespace cand {

using Tri = df::Tri2;

// list of CGAL vertex handles
static std::vector<Tri::Vertex_handle>
make_id_to_vertex_handle(const Tri& tri) {
  std::vector<Tri::Vertex_handle> id2vh;
  id2vh.reserve(static_cast<size_t>(
      std::distance(tri.finite_vertices_begin(), tri.finite_vertices_end())));
  for (auto v = tri.finite_vertices_begin(); v != tri.finite_vertices_end(); ++v)
    id2vh.push_back(v);
  return id2vh;
}

// checks if locally non-regular edges (down fliapabble) are not already in the target triangulation
std::vector<df::reg::EdgeRecord>
down_not_in_target(const Tri& tri,
                   const df::OmegaQuad& omega,
                   const df::target::EdgeSet& target_edges)
{

  // collect all edges w/ flags & polyscope IDs 
  auto rows = df::reg::collect_flippable_edges_with_ids(tri, omega);

  // Map planar ID -> CGAL vertex (to get coordinates in this triangulation)
  auto id2vh = make_id_to_vertex_handle(tri);

  std::vector<df::reg::EdgeRecord> out;
  out.reserve(rows.size());

  for (const auto& e : rows) {
    if (!(e.nonreg_down && !e.degenerate)) continue;
    // create edge key from the two endpoints
    const df::P2& A = id2vh[e.ia]->point(); // vertex handle of first vertex (polyscope ID map to vertex handle)
    const df::P2& B = id2vh[e.ib]->point(); // vertex handle of second vertex (polyscope ID map to vertex handle)
    auto key = df::target::make_edge_key(A, B); // create ordered edge key

    // if it's already in target triangulation, drop it 
    if (target_edges.find(key) != target_edges.end()) continue;

    out.push_back(e);
  }
  return out;
}

void print_candidates(const std::vector<df::reg::EdgeRecord>& rows) {
  using std::left;
  using std::setw;

  std::cout << "=== down-flip candidates not in target triangulation ===\n";
  std::cout << left
            << setw(6)  << "start"
            << setw(6)  << "end"
            << setw(14) << "start_x"
            << setw(14) << "start_y"
            << setw(14) << "end_x"
            << setw(14) << "end_y"
            << "\n";

  for (const auto& e : rows) {
    std::cout << left
              << setw(6)  << e.ia
              << setw(6)  << e.ib
              << setw(14) << e.a2.x()
              << setw(14) << e.a2.y()
              << setw(14) << e.b2.x()
              << setw(14) << e.b2.y()
              << "\n";
  }

  if (rows.empty()) std::cout << "(none)\n";
  std::cout << "=== end ===\n";
}

}} // namespace


#include "target_edges.h"
namespace df { namespace target {

// create an edge key (ordered pair of P2) from two points u and v
EdgeKey make_edge_key(const df::P2& u, const df::P2& v) {
  static const P2Less less{};
  return less(u, v) ? EdgeKey{u, v} : EdgeKey{v, u};
}

// build edge set from all edges in target triangulation 
EdgeSet build_edge_set_from_triangulation(const df::Tri2& Tv) {
  EdgeSet S;
  // loop over all finite faces to extract edges
  for (auto f = Tv.finite_faces_begin(); f != Tv.finite_faces_end(); ++f) {
    const df::P2& a = f->vertex(0)->point();
    const df::P2& b = f->vertex(1)->point();
    const df::P2& c = f->vertex(2)->point();
    // insert edges into set -> duplicates automatically handled by set
    S.insert(make_edge_key(a,b));
    S.insert(make_edge_key(b,c));
    S.insert(make_edge_key(c,a));
  }
  return S;
}

}}


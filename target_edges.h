#pragma once
#include "input.h"
#include <set>

namespace df { namespace target {

// comparison operator for P2 points
struct P2Less {
  bool operator()(const df::P2& a, const df::P2& b) const {
    auto cx = CGAL::compare(a.x(), b.x());
    if (cx != CGAL::EQUAL) return cx == CGAL::SMALLER;
    auto cy = CGAL::compare(a.y(), b.y());
    return cy == CGAL::SMALLER;
  }
};

// Edge key: ordered pair of P2 points representing an edge
using EdgeKey = std::pair<df::P2, df::P2>;

// Comparison operator for EdgeKey to be used in std::set
struct EdgeKeyLess {
  bool operator()(const EdgeKey& e1, const EdgeKey& e2) const {
    static const P2Less less{};
    if (less(e1.first, e2.first))  return true;
    if (less(e2.first, e1.first))  return false;
    return less(e1.second, e2.second);
  }
};

// Edge set: set of EdgeKeys
using EdgeSet = std::set<EdgeKey, EdgeKeyLess>;

// build edge set from all edges in target triangulation
EdgeSet build_edge_set_from_triangulation(const df::Tri2& Tv);
// create an edge key (ordered pair of P2) from two points u and v
EdgeKey make_edge_key(const df::P2& u, const df::P2& v);

}}

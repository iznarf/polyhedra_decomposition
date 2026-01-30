#pragma once

#include "poset_view.h"
#include "poset_utils.h"

#include <vector>
#include <utility>

namespace pst_interval {

// ------------------------------------------------------------
// helpers
// ------------------------------------------------------------

// Build a bitset mask from a node list
inline pst::Bitset make_mask_bitset(int n, const std::vector<int>& nodes) {
  pst::Bitset in(n);
  for (int v : nodes) {
    if (0 <= v && v < n) in.set(v);
  }
  return in;
}

// ------------------------------------------------------------
// order-derived sets (bitset versions, optimized via reachability)
// ------------------------------------------------------------

// Upper set ↑x = { z | x <= z }
inline pst::Bitset upper_set(const PosetView& P, int x) {
  pst::Bitset out(P.n);
  if (x < 0 || x >= P.n) return out;

  if (P.reachability_up) {
    out = (*P.reachability_up)[x];
    return out;
  }

  // fallback: compute using leq
  for (int z = 0; z < P.n; ++z)
    if (leq(P, x, z)) out.set(z);
  return out;
}

// Lower set ↓y = { z | z <= y }
inline pst::Bitset lower_set(const PosetView& P, int y) {
  pst::Bitset out(P.n);
  if (y < 0 || y >= P.n) return out;

  if (P.reachability_down) {
    out = (*P.reachability_down)[y];
    return out;
  }

  // fallback: compute using leq
  for (int z = 0; z < P.n; ++z)
    if (leq(P, z, y)) out.set(z);
  return out;
}


// Interval [x,y] = { z | x <= z <= y } = ↑x ∩ ↓y
// Empty if x !<= y
inline pst::Bitset interval_xy(const PosetView& P, int x, int y) {
    pst::Bitset out(P.n);

    if (x < 0 || x >= P.n || y < 0 || y >= P.n) return out;
    if (!leq(P, x, y)) return out;

    auto up = upper_set(P, x);   // ↑x
    auto dn = lower_set(P, y);   // ↓y
    up &= dn;                    // ↑x ∩ ↓y
    return up;
}

// ------------------------------------------------------------
// bitset -> node list
// ------------------------------------------------------------

inline std::vector<int> nodes_from_bitset(const pst::Bitset& bs) {
    std::vector<int> nodes;
    for (int i = 0; i < bs.size(); ++i) {
        if (bs.test(i)) nodes.push_back(i);
    }
    return nodes;
}

// ------------------------------------------------------------
// cover edges in induced subgraph
// ------------------------------------------------------------

// Uses cover_up, restricted to nodes in mask
inline std::vector<std::pair<int,int>>
induced_cover_up_edges(const PosetView& P, const pst::Bitset& mask) {
    std::vector<std::pair<int,int>> edges;
    if (!P.cover_up) return edges;

    for (int u = 0; u < P.n; ++u) {
        if (!mask.test(u)) continue;

        for (int v : (*P.cover_up)[u]) {
        if (mask.test(v)) {
            edges.emplace_back(u, v);
        }
        }
    }
    return edges;
}

} 

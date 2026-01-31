#include "dedekind_cut.h"

#include <algorithm>
#include <functional>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_set>

namespace dm_completion {

// -------------------------
// hashing for vector<int>
// -------------------------
struct VecHash {
  size_t operator()(const std::vector<int>& v) const noexcept {
    size_t h = 1469598103934665603ull;
    for (int x : v) {
      h ^= (size_t)x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    }
    return h;
  }
};

// -------------------------
// printing
// -------------------------
void print_cut(int idx, const DedekindCut& C) {
  std::cout << "  [" << idx << "] gens={";
  for (int i = 0; i < (int)C.generators.size(); ++i) {
    std::cout << C.generators[i] << (i + 1 < (int)C.generators.size() ? "," : "");
  }
  std::cout << "}  minF={";
  for (int i = 0; i < (int)C.minimal_elements_F.size(); ++i) {
    std::cout << C.minimal_elements_F[i] << (i + 1 < (int)C.minimal_elements_F.size() ? "," : "");
  }
  std::cout << "}  maxI'={";
  for (int i = 0; i < (int)C.maximal_elements_Iprime.size(); ++i) {
    std::cout << C.maximal_elements_Iprime[i] << (i + 1 < (int)C.maximal_elements_Iprime.size() ? "," : "");
  }
  std::cout << "}  |F|=" << C.F.size() << " |I'|=" << C.Iprime.size() << "\n";
}

// -------------------------
// core helpers 
// -------------------------

// Minimal generators: for an antichain, this is already itself
// right now this functon is not necessary
static std::vector<int> normalize_generators_antichain(std::vector<int> gens) {
  std::sort(gens.begin(), gens.end());
  gens.erase(std::unique(gens.begin(), gens.end()), gens.end());
  return gens;
}


// Build the full cut from a generator antichain
DedekindCut build_cut(const PosetView& P, const std::vector<int>& gens_in) {
  DedekindCut C;
  C.generators = normalize_generators_antichain(gens_in);

  // Filter F = common upper bounds of generators
  pst::Bitset Fmask = meet_join::common_upper_bounds(P, C.generators);
  C.F = pst_interval::nodes_from_bitset(Fmask);

  // min(F)
  pst::Bitset minFmask = meet_join::minimal_elements(P, Fmask);
  C.minimal_elements_F = pst_interval::nodes_from_bitset(minFmask);

  // I' = common lower bounds of min(F)
  pst::Bitset Ipmask = meet_join::common_lower_bounds(P, C.minimal_elements_F);
  C.Iprime = pst_interval::nodes_from_bitset(Ipmask);

  // max(I')
  pst::Bitset maxIpmask = meet_join::maximal_elements(P, Ipmask);
  C.maximal_elements_Iprime = pst_interval::nodes_from_bitset(maxIpmask);

  return C;
}

// Fast key: min(F) sorted
static std::vector<int> cut_key_minF(const PosetView& P, const std::vector<int>& gens_in) {
  std::vector<int> gens = normalize_generators_antichain(gens_in);

  pst::Bitset Fmask = meet_join::common_upper_bounds(P, gens);
  pst::Bitset minFmask = meet_join::minimal_elements(P, Fmask);

  auto key = pst_interval::nodes_from_bitset(minFmask);
  std::sort(key.begin(), key.end());
  return key;
}

// -------------------------
// width via Dilworth (Hopcroft–Karp)
// -------------------------

int width_longest_antichain(const PosetView& P) {
  const int N = P.n;

  // Build bipartite graph: edge u->v iff u < v
  std::vector<std::vector<int>> adjL(N);

  // If we have reachability_down[v] = { u | u <= v }, we can build edges efficiently
  if (P.reachability_down) {
    for (int v = 0; v < N; ++v) {
      pst::Bitset bs = (*P.reachability_down)[v]; // u <= v
      if (v >= 0 && v < (int)bs.size()) bs.reset(v); // strict: u < v

      for (int u = bs.find_first();
           u != pst::Bitset::npos;
           u = bs.find_next(u))
      {
        adjL[(int)u].push_back(v);
      }
    }
  } else {
    // fallback O(n^2): add u->v if u < v
    for (int u = 0; u < N; ++u) {
      for (int v = 0; v < N; ++v) {
        if (u == v) continue;
        if (leq(P, u, v) && !leq(P, v, u)) adjL[u].push_back(v);
      }
    }
  }

  std::vector<int> pairU(N, -1), pairV(N, -1), dist(N);

  auto bfs = [&]() {
    std::queue<int> q;
    for (int u = 0; u < N; ++u) {
      if (pairU[u] == -1) { dist[u] = 0; q.push(u); }
      else dist[u] = -1;
    }
    bool found_free = false;
    while (!q.empty()) {
      int u = q.front(); q.pop();
      for (int v : adjL[u]) {
        int u2 = pairV[v];
        if (u2 == -1) {
          found_free = true;
        } else if (dist[u2] == -1) {
          dist[u2] = dist[u] + 1;
          q.push(u2);
        }
      }
    }
    return found_free;
  };

  std::function<bool(int)> dfs = [&](int u) {
    for (int v : adjL[u]) {
      int u2 = pairV[v];
      if (u2 == -1 || (dist[u2] == dist[u] + 1 && dfs(u2))) {
        pairU[u] = v;
        pairV[v] = u;
        return true;
      }
    }
    dist[u] = -1;
    return false;
  };

  int matching = 0;
  while (bfs()) {
    for (int u = 0; u < N; ++u) {
      if (pairU[u] == -1 && dfs(u)) matching++;
    }
  }

  // Dilworth: width = N - |matching|
  return N - matching;
}

// -------------------------
// enumerate all cuts
// -------------------------

std::vector<DedekindCut> compute_all_cuts(const PosetView& P) {
  const int N = P.n;
  const int w = width_longest_antichain(P);
  std::cout << "[dedekind_cut] poset width (longest antichain size) = " << w << "\n";

  std::vector<DedekindCut> cuts;

  // Key cuts by min(F), which uniquely identifies the DM cut.
  std::unordered_set<std::vector<int>, VecHash> seen;

  std::vector<int> cur;
  int new_this_k = 0;

  auto incomparable_to_cur = [&](int v) -> bool {
    for (int a : cur) {
      if (leq(P, a, v) || leq(P, v, a)) return false;
    }
    return true;
  };

  std::function<void(int,int)> backtrack_k = [&](int start, int k) {
    if ((int)cur.size() == k) {
      // FAST: compute key only
      std::vector<int> key = cut_key_minF(P, cur);

      if (!seen.insert(key).second) return;

      // SLOW: build full cut only if new
      DedekindCut C = build_cut(P, cur);
      print_cut((int)cuts.size(), C);
      cuts.push_back(std::move(C));
      new_this_k++;
      return;
    }

    for (int v = start; v < N; ++v) {
      if (!incomparable_to_cur(v)) continue;
      cur.push_back(v);
      backtrack_k(v + 1, k);
      cur.pop_back();
    }
  };

  for (int k = 1; k <= w; ++k) {
    std::cout << "\n=== k = " << k << " generators (antichains) ===\n";
    new_this_k = 0;
    cur.clear();
    backtrack_k(0, k);
    std::cout << "New cuts for k=" << k << ": " << new_this_k << "\n";

    // early stop if no new cuts found
    if (new_this_k == 0) break;
  }

  std::cout << "[dedekind_cut] total Dedekind cuts found: " << cuts.size() << "\n";
  return cuts;
}

} // namespace dm


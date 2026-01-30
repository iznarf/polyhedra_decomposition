#pragma once
#include "poset_view.h"
#include "poset_utils.h"
#include "poset_interval.h"

#include <vector>
#include <optional>

namespace meet_join {

struct Result {
  // candidates: the set of possible meets/joins (unique iff size==1)
  std::vector<int> candidates;

  bool exists() const { return !candidates.empty(); }
  bool unique() const { return candidates.size() == 1; }
  std::optional<int> value() const {
    if (candidates.size() == 1) return candidates[0];
    return std::nullopt;
  }
};

// ---------------------------------------------
// Helpers
// ---------------------------------------------

inline std::vector<int> to_vector(const pst::Bitset& S) {
  return pst_interval::nodes_from_bitset(S);
}

// Common upper bounds of a set A:  ⋂_{a in A} ↑a
// elems: list of node ids (can be size 0, 1, 2, ...)
inline pst::Bitset common_upper_bounds(const PosetView& P, const std::vector<int>& elems) {
  pst::Bitset U(P.n);

  if (elems.empty()) {
    // Convention: upper bounds of empty set = all elements
    U.set();
    return U;
  }

  // start with ↑elems[0]
  int a0 = elems[0];
  if (a0 < 0 || a0 >= P.n) return pst::Bitset(P.n); // empty
  U = pst_interval::upper_set(P, a0);

  for (int i = 1; i < (int)elems.size(); ++i) {
    int a = elems[i];
    if (a < 0 || a >= P.n) return pst::Bitset(P.n); // empty
    U &= pst_interval::upper_set(P, a);
  }
  return U;
}

// Common lower bounds of a set A:  ⋂_{a in A} ↓a
inline pst::Bitset common_lower_bounds(const PosetView& P, const std::vector<int>& elems) {
  pst::Bitset D(P.n);

  if (elems.empty()) {
    // Convention: lower bounds of empty set = all elements
    D.set();
    return D;
  }

  int a0 = elems[0];
  if (a0 < 0 || a0 >= P.n) return pst::Bitset(P.n); // empty
  D = pst_interval::lower_set(P, a0);

  for (int i = 1; i < (int)elems.size(); ++i) {
    int a = elems[i];
    if (a < 0 || a >= P.n) return pst::Bitset(P.n); // empty
    D &= pst_interval::lower_set(P, a);
  }
  return D;
}

// Minimal elements of S (no strictly smaller element in S)
inline pst::Bitset minimal_elements(const PosetView& P, const pst::Bitset& S) {
  pst::Bitset out(P.n);

  // Fast path: use reachability_down[a] = { b | b <= a }
  if (P.reachability_down) {
    for (int a = 0; a < P.n; ++a) {
      if (!S.test(a)) continue;

      pst::Bitset smaller = (*P.reachability_down)[a] & S;
      smaller.reset(a);                 // ignore itself (works even if diagonal absent)
      if (!smaller.any()) out.set(a);   // no smaller element in S => minimal
    }
    return out;
  }

  // Fallback O(n^2) via leq
  for (int a = 0; a < P.n; ++a) {
    if (!S.test(a)) continue;

    bool hasSmaller = false;
    for (int b = 0; b < P.n; ++b) {
      if (a == b || !S.test(b)) continue;
      if (leq(P, b, a) && !leq(P, a, b)) { hasSmaller = true; break; } // b < a
    }
    if (!hasSmaller) out.set(a);
  }
  return out;
}

// Maximal elements of S (no strictly larger element in S)
inline pst::Bitset maximal_elements(const PosetView& P, const pst::Bitset& S) {
  pst::Bitset out(P.n);

  // Fast path: use reachability_up[a] = { b | a <= b }
  if (P.reachability_up) {
    for (int a = 0; a < P.n; ++a) {
      if (!S.test(a)) continue;

      pst::Bitset bigger = (*P.reachability_up)[a] & S;
      bigger.reset(a);
      if (!bigger.any()) out.set(a);
    }
    return out;
  }

  // Fallback O(n^2) via leq
  for (int a = 0; a < P.n; ++a) {
    if (!S.test(a)) continue;

    bool dominated = false;
    for (int b = 0; b < P.n; ++b) {
      if (a == b || !S.test(b)) continue;
      if (leq(P, a, b) && !leq(P, b, a)) { dominated = true; break; } // a < b
    }
    if (!dominated) out.set(a);
  }
  return out;
}

// ---------------------------------------------
// Meet / Join for ANY number of inputs
// ---------------------------------------------

// join(A) = minimal elements of common_upper_bounds(A)
// returns candidates (unique iff lattice join exists for A)
inline Result join(const PosetView& P, const std::vector<int>& elems) {
  Result r;
  pst::Bitset U = common_upper_bounds(P, elems);
  pst::Bitset mins = minimal_elements(P, U);
  r.candidates = to_vector(mins);
  return r;
}

// meet(A) = maximal elements of common_lower_bounds(A)
inline Result meet(const PosetView& P, const std::vector<int>& elems) {
  Result r;
  pst::Bitset D = common_lower_bounds(P, elems);
  pst::Bitset maxs = maximal_elements(P, D);
  r.candidates = to_vector(maxs);
  return r;
}

// Convenience overloads for 2 inputs (nice for UI)
inline Result join(const PosetView& P, int a, int b) {
  return join(P, std::vector<int>{a, b});
}
inline Result meet(const PosetView& P, int a, int b) {
  return meet(P, std::vector<int>{a, b});
}

} // namespace meet_join


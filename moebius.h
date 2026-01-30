#pragma once

#include <cstdint>
#include <vector>

struct PosetView;

namespace pst {

// Compute μ(global_min, x) for all x.
//
// Requirements on PosetView:
//   - cover_down          (to detect global minimum)
//   - topo_up             (bottom -> top order of cover_up)
//   - reachability_down   (bitset: down[b] contains all a with a <= b)
//
// Returns:
//   mu[x] = μ(global_min, x)
//   If there is no UNIQUE global minimum, returns all zeros.
std::vector<std::int64_t>
mobius_from_global_min(const PosetView& P);

} // namespace pst




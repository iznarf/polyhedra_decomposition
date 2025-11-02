#pragma once
#include "input.h"
#include "regularity_check.h"
#include "target_edges.h"
#include <vector>

namespace df { namespace cand {

// Return only edges that are:
//  - locally non-regular DOWN
//  - non-degenerate
//  - NOT already an edge of the target triangulation Tv
std::vector<df::reg::EdgeRecord>
down_not_in_target(const df::Tri2& tri,
                   const df::OmegaQuad& omega,
                   const df::target::EdgeSet& target_edges);

// Convenience: print a compact table of the filtered candidates
void print_candidates(const std::vector<df::reg::EdgeRecord>& rows);

}} // namespace


#pragma once

#include "poset_view.h"
#include "poset_utils.h"
#include "poset_interval.h"
#include "meet_join.h"

#include <vector>

namespace dm_completion {

struct DedekindCut {
  // the antichain we used to generate it (your backtracking "cur")
  std::vector<int> generators;

  // filter and I' as node lists (global ids); not really necessary to store 
  std::vector<int> F;
  std::vector<int> Iprime;

  // canonical boundary antichains
  std::vector<int> minimal_elements_F;       // min(F)
  std::vector<int> maximal_elements_Iprime;  // max(I')
};

DedekindCut build_cut(const PosetView& P, const std::vector<int>& gens_in);

void print_cut(int idx, const DedekindCut& C);

// width(P) = size of longest antichain
int width_longest_antichain(const PosetView& P);

// Enumerate all Dedekind–MacNeille cuts by backtracking over antichains
std::vector<DedekindCut> compute_all_cuts(const PosetView& P);




} // namespace dm


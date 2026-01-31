#pragma once
#include "poset.h"
#include "poset_utils.h"
#include "input.h"
#include <vector>


namespace pst2 {


struct Poset2 {
    std::vector<std::vector<int>> cover_down;  // covering edges u -> v where u < v and no w with u < w < v
    std::vector<std::vector<int>> cover_up;    // covering edges v -> u where u < v and no w with u < w < v
    std::vector<int> topo_up;                 // topo_sort: node indices in topological order
    std::vector<int> topo_down;               // topo_sort: node indices in topological order
    std::vector<pst::Bitset> reachability_down;  // reachability bitsets in down direction for <=2; // reachability_down[u][v]=1 iff u ->* v , i.e. v <= u
    std::vector<pst::Bitset> reachability_up;    // reachability bitsets in up direction for <=2;   // reachability_up[u][v]=1 iff v ->* u , i.e. u <= v
    std::vector<int> levels;                   // level[u] = length of longest chain from minimal element to u
};

Poset2 build_poset2(const df::InputData& D, const pst::Poset1& P1);

bool leq(const Poset2& P, int a, int b);

void debug_compare_poset1(const df::InputData& D, const pst::Poset1& P1);

} 

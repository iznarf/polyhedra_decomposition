#pragma once
#include "poset.h"
#include "poset_utils.h"
#include "input.h"
#include <vector>


namespace pst2 {


struct Poset2 {
    std::vector<std::vector<int>> cover_down;  // covering edges u -> v where u < v and no w with u < w < v
    std::vector<std::vector<int>> cover_up;    // covering edges v -> u where u < v and no w with u < w < v
    std::vector<pst::Bitset> reachability_down;  // reachability bitsets in down direction for <=2; // reachability_down[u][v]=1 iff u ->* v (DOWN), i.e. v <= u
};

Poset2 build_poset2(const df::InputData& D, const pst::Poset1& P1);

bool leq(const Poset2& P, int a, int b);

void debug_compare_poset1(const df::InputData& D, const pst::Poset1& P1);

std::vector<char> upper_set(const Poset2& P, int x);
std::vector<char> lower_set(const Poset2& P, int x);

std::vector<int> interval_xy(const Poset2& P, int x, int y);
std::vector<int> interval_min_x(const Poset2& P, int x);
std::vector<int> interval_x_max(const Poset2& P, int x);
std::vector<int> meet_candidates_xy(const Poset2& P, int x, int y);
std::vector<int> join_candidates_xy(const Poset2& P, int x, int y);

std::vector<int> minima(const Poset2& P);
std::vector<int> maxima(const Poset2& P);
int unique_min_or_minus1(const Poset2& P);
int unique_max_or_minus1(const Poset2& P);

} 

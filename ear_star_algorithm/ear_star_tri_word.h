#pragma once

#include "input.h"
#include "poset_ear_star.h"

#include <string>
#include <vector>

namespace earstar {

// compute the ear-star word for a triangulation on boundary vertices 0..n-1
// with fixed base edge (0, n-1) and interior vertex star_id.
//
// returns a token list like {"3","1","*","6","2"} (numbers as strings + "*").
//
std::vector<std::string> triangulation_to_word(
    const df::Tri2& tri,
    int n_boundary,
    df::vertex_id star_id
);

// join tokens by spaces
std::string word_to_string(const std::vector<std::string>& w);

// for every node in the flip poset, reconstruct triangulation and print its word
void print_words_for_poset(
    const df::Tri2& tri_start,
    const pst_es::FlipPoset& P,
    int n_boundary,
    df::vertex_id star_id,
    const df::P2& star_point
);

} // namespace earstar
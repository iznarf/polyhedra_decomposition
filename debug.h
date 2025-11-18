#pragma once

#include "input.h"   
#include <vector>
#include <utility>

namespace df {

    // Return all global edges that are in current but NOT in target.
    // Edges are treated as undirected (a,b) == (b,a), and stored with (min,max).
    std::vector<std::pair<vertex_id, vertex_id>>
    edges_in_current_not_in_target(const Tri2& current, const Tri2& target);

    // compare D.tri_current vs D.tri_lower and print the differences.
    void debug_print_edge_diff_current_vs_lower(const InputData& D);

    void debug_print_edge_list(const df::InputData& D);


    void debug_print_current_vertex_ids(const df::InputData& D,
                                     const std::vector<int>& local_indices);



    void debug_print_edge_diff_with_local(const df::InputData& D);

    void debug_edges_created_by_insertion(df::vertex_id id, const df::InputData& D);

    void debug_check_edge_against_lower(df::vertex_id i,
                                           df::vertex_id j,
                                           const df::InputData& D);



} // namespace df

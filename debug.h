#pragma once

#include "input.h"   
#include <vector>
#include <array>
#include <utility>


namespace df {

    enum class DebugTetKind {
        EdgeFlip,   
        Insertion   
    };

    struct DebugTetrahedron {
        std::array<vertex_id, 4> verts; // global ids of the four vertices
        DebugTetKind kind;
    };

    // function:
    //  - search for locally non regular edges in current triangulation
    //  - take the four involved vertices and store them
    //  - take the missing vertices, find the triangle they are in, and
    //    store the triangle vertices + missing vertex
    //  - output: list of tetrahedra + kind (2-2 flip vs insertion)
    std::vector<DebugTetrahedron>
    collect_debug_tetrahedra(const InputData& D);



    enum class TriKind {
        Current,
        Lower
    };

    // returns all global edges that are in current but NOT in target 
    std::vector<std::pair<vertex_id, vertex_id>>
    edges_in_current_not_in_target(const Tri2& current, const Tri2& target);

    void debug_print_edge_list(const df::InputData& D);


    void debug_print_current_vertex_ids(const df::InputData& D,
                                     const std::vector<int>& local_indices);



    bool edge_diff_with_lower(const df::InputData& D);

    void debug_edges_created_by_insertion(df::vertex_id id, const df::InputData& D);

    void debug_check_edge_against_lower(df::vertex_id i,
                                           df::vertex_id j,
                                           const df::InputData& D);

    void debug_print_local_to_global_map(const df::InputData& D, TriKind which, const std::vector<int>& local_indices = {});

    void print_step_history(const df::InputData& D);



} // namespace df

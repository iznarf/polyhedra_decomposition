#include "input.h"
#include "check_edges.h"
#include "conforming.h"
#include "flip.h"
#include "visualization.h"
#include <iostream>
#include <unordered_map>

namespace df { 

void perform_all_conforming_down_flips(df::InputData& in) {
    const df::Tri2& target_triangulation  = in.tri_lower;
    df::Tri2&       current_triangulation = in.tri_current;

    // global vertex_id -> local polyscope index of target triangulation
    std::unordered_map<df::vertex_id, int> id2local_lower;
    id2local_lower.reserve(target_triangulation.number_of_vertices());
    int loc = 0;
    for (auto v = target_triangulation.finite_vertices_begin();
         v != target_triangulation.finite_vertices_end(); ++v) {
        id2local_lower[v->info()] = loc++;
    }

    // global vertex_id -> local polyscope index of current triangulation
    std::unordered_map<df::vertex_id, int> id2local_current;

    int iteration = 0;
    while (true) {

        // rebuild id2local_current for the *current* triangulation
        id2local_current.clear();
        id2local_current.reserve(current_triangulation.number_of_vertices());
        loc = 0;
        for (auto v = current_triangulation.finite_vertices_begin();
             v != current_triangulation.finite_vertices_end(); ++v) {
            id2local_current[v->info()] = loc++;
        }

        // collect locally non-regular (down-flippable) edges
        std::vector<std::array<df::vertex_id,2>> down_flip_edges = df::reg::find_locally_non_regular_edges(current_triangulation);

        // for storing flippable edges which are NOT in the target triangulation and conforming 
        std::vector<std::pair<df::vertex_id, df::vertex_id>> flippable_edges;
        flippable_edges.reserve(down_flip_edges.size());

        std::cout << "\n=== Flip iteration " << iteration
                  << " : locally non-regular (down-flippable) edges ===\n";
        std::cout << "count: " << down_flip_edges.size() << "\n";
        std::cout << "global indices -> local indices  in target  conforming \n";

        // iterate over all down-flippable edges
        for (size_t k = 0; k < down_flip_edges.size(); ++k) {
            const auto [ia, ib] = down_flip_edges[k];

            // local indices for polyscope
            int local_index_a = id2local_current.at(ia);
            int local_index_b = id2local_current.at(ib);
            std::ostringstream oss;
            oss << "(" << local_index_a << "," << local_index_b << ")";
            std::string pl = oss.str();

            // compute flags
            bool in_target = df::reg::edge_in_target(ia, ib, in);
            bool conf      = df::reg::is_flip_conforming(ia, ib, in,
                                                         id2local_current,
                                                         id2local_lower);

            // right now also edges which are already in target are checked for conformity, this is unnecessary    
            // in_target is false and conf is true -> we can flip this edge
            if (!in_target && conf) {
                flippable_edges.push_back({ia, ib});
            }

            // print one line per edge
            std::cout << std::setw(3) << k << " (" << ia << "," << ib << ")->"
                      << std::setw(8) << pl << "  "
                      << (in_target ? "T" : "F") << "         "
                      << (conf ? "T" : "F") << "\n";
        }

        // if no conforming flippable edges left -> stop
        if (flippable_edges.empty()) {
            std::cout << "\nNo conforming flippable edges left in this phase.\n";
            break;
        }

        // apply an edge flip to the first element in the list of flippable edges
        std::cout << "\n=== Applying conforming flip (iteration " << iteration << ") ===\n";

        const auto [ia, ib] = flippable_edges.front();
        int local_index_a = id2local_current.at(ia);
        int local_index_b = id2local_current.at(ib);
        std::cout << "flipping edge (local " << local_index_a << "," << local_index_b
                  << ")  (global " << ia << "," << ib << ")\n";

        df::apply_edge_flip(ia, ib, in);

        // update the visualization 
        viz::show_or_update_current(in);

        viz::debug_print_edge_list(in);

        ++iteration;
    }
}

} 

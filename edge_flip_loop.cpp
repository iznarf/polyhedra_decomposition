#include "input.h"
#include "check_edges.h"
#include "conforming.h"
#include "flip.h"
#include "visualization.h"
#include "debug.h"
#include <iostream>
#include <unordered_map>

namespace df { 

void perform_all_conforming_down_flips(df::InputData& in) {
    const df::Tri2& target_triangulation  = in.tri_lower;
    df::Tri2&       current_triangulation = in.tri_current;

    int iteration = 0;
    while (true) {

        // collect locally non-regular (down-flippable) edges
        std::vector<std::array<df::vertex_id,2>> down_flip_edges = df::reg::find_locally_non_regular_edges(current_triangulation);

        // for storing flippable edges which are NOT in the target triangulation and conforming 
        std::vector<std::pair<df::vertex_id, df::vertex_id>> flippable_edges;
        flippable_edges.reserve(down_flip_edges.size());


        // iterate over all down-flippable edges
        for (size_t k = 0; k < down_flip_edges.size(); ++k) {
            const auto [ia, ib] = down_flip_edges[k];

            // compute flags
            bool in_target = df::reg::edge_in_target(ia, ib, in);
            bool conf      = false; 

            if (in_target == false) {
                conf = df::reg::is_flip_conforming(ia, ib, in, current_triangulation);
                if (conf == true) {
                    flippable_edges.push_back({ia, ib});

                }
            }

        }

        // if no conforming flippable edges left -> stop
        if (flippable_edges.empty()) {
            std::cout << "\nNo conforming flippable edges left in this phase.\n";
            break;
        }

        // print flippable edges
        std::cout << "\nConforming flippable edges (global ids): ";
        for (const auto& [ia, ib] : flippable_edges) {
            std::cout << "(" << ia << "," << ib << ") ";
        }
        std::cout << "\n";


        // apply an edge flip to the first element in the list of flippable edges
        std::cout << "\n=== Applying conforming flip (iteration " << iteration << ") ===\n";

        const auto [ia, ib] = flippable_edges.front();
        std::cout << "flipping edge (" << ia << "," << ib << ")\n";

        df::apply_edge_flip(ia, ib, in);

        // update the visualization 
        viz::show_or_update_current(in);

        //df::debug_print_edge_list(in);

        ++iteration;
    }
}

} 

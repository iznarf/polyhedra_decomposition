#include <polyscope/polyscope.h>

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <cstdint>
#include <sstream>

#include "input.h"
#include "visualization.h"
#include "edge_flip_loop.h"
#include "insertion.h"
#include "flip.h"
#include "conforming_insertion.h"
#include "debug.h"
#include "geometry_utils.h"
#include "replay.h"
#include "ui_callbacks.h"


int main() {
    polyscope::init();

    // number of vertices in triangulation
    int n_points = 21;
    // random seed to start point generation
    unsigned seed0 = 44;

    df::InputData in = df::make_random_valid_input(n_points, seed0);

    // current setup: we do not sort the to be inserted vertices by height 
    // we do a conforming insertion check 

    // these are valid inputs where the algorithm works:
    //df::InputData in = daf::make_random_input(37, 4458);
    //df::InputData in = df::make_random_input(19, 42);
    //df::InputData in = df::make_random_input(35, 218);
    //df::InputData in = df::make_random_input(40, 40);
    //df::InputData in = df::make_random_input(50, 40);


    // these are valid inputs where the algorithm with complicated conforming insertion check stucks BUT works if we do simple conforming check only:
    //df::InputData in = df::make_random_input(21, 44); 

    // -> the order of vertex insertions matters for the algorithm to succeed
    // -> sorting vertices by descending height is not correct 
    // what is the 'correct' order of vertex insertions?

    viz::register_triangulation_as_mesh(in.tri_lower, in.points2d, "lower 2D", "lower lifted");
    viz::register_triangulation_as_mesh(in.tri_upper, in.points2d, "upper 2D", "upper lifted");

    viz::show_or_update_current(in);


    while (true) {

        // perform all conforming down-flips for the current triangulation until flip list is empty
        df::perform_all_conforming_down_flips(in);

        // see which vertices of the lower triangulation are still missing
        auto missing = df::find_missing_vertices(in.tri_current, in.tri_lower);


        if (missing.empty()) { 
            std::cout << "\n=== no missing vertices left ===\n";
            break; 
        }

        // print missing vertex ids in current triangulation
        std::cout << "\nMissing vertices (global ids): ";
        for (auto id : missing) std::cout << id << " ";
        std::cout << "\n";

        // pick candidates sorted by height (highest first)
        //std::vector<df::vertex_id> insertion_vertex_list =
        //df::sorted_insertion_vertices(missing, in);  // global ids sorted by height

        std::vector<df::vertex_id> insertion_vertex_list = missing; // unsorted version


        df::vertex_id insertion_vertex = 0; // will only be used if we find a conforming one
        bool found_conforming = false;

        // try candidates in order until one is BOTH a downflip and conforming
        for (auto id : insertion_vertex_list) {

            // 1) check downflip condition
            if (!df::is_insertion_downflip(id, in)) {
                std::cout << "[main] skipping vertex " << id
                        << " (insertion is not a down-flip)\n";
                continue;
            }

        
            // 2) check global conformance w.r.t. lower triangulation
            if (df::reg::is_insertion_conforming(id, in)) {
                insertion_vertex = id;
                found_conforming = true;
                break;
            } else {
                std::cout << "[main] WARNING: insertion of vertex " << id
                        << " is non-conforming, trying next candidate\n";
            }
        }

        // if no conforming insertion exists → polyhedron is non-decomposable
        if (!found_conforming) {
            std::cout << "\n[main] ERROR: all candidate vertex insertions are non-conforming.\n"
                    << "[main] The polyhedron appears to be non-decomposable.\n";
            break;  // break out of the main while(true) loop
        }

        std::cout << "Inserting vertex (global id) " << insertion_vertex << "\n";
        df::apply_vertex_insertion(insertion_vertex, in);

        viz::show_or_update_current(in);
        df::debug_print_edge_list(in);
    }
    
    // compare triangulations now 
    if (df::edge_diff_with_lower(in) == true){
        std::cout << "\n[main] ERROR: after all flips and insertions, current triangulation differs from lower triangulation!\n";
        
    } else {
        std::cout << "\n[main] SUCCESS: current triangulation matches lower triangulation!\n";
    }

    //df::debug_print_local_to_global_map(in, df::TriKind::Lower);
    //df::debug_print_local_to_global_map(in, df::TriKind::Current);


    df::print_step_history(in);

    std::vector<df::DebugTetrahedron> debug_tets = df::collect_debug_tetrahedra(in);

    // visualize debug tetrahedra
    viz::load_debug_tetrahedra(in, debug_tets);
    

    // initialize replay data
    df::init_replay(in);
   
    polyscope::state::userCallback = combined_ui_callback;


    // to do: visualize the tetrahedra corresponding to the steps taken by the algorithm



    polyscope::show();
    return 0;
}







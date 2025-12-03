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
#include "poset.h"
#include "flip.h"



int main() {
    polyscope::init();

    // number of vertices in triangulation
    int n_points = 31;
    // random seed to start point generation
    unsigned seed0 = 1312;

    df::InputData in = df::make_random_valid_input(n_points, seed0);


    // these are valid inputs where the algorithm works:
    //df::InputData in = df::make_random_input(21, 44); 
    //df::InputData in = df::make_random_input(19, 42);
    //df::InputData in = df::make_random_input(15,43);
    //df::InputData in = df::make_random_input(18, 23);
    //df::InputData in = df::make_random_input(15, 23);
    //df::InputData in = df::make_random_input(25, 84);
    //df::InputData in = df::make_random_input(29, 4234);
    //df::InputData in = df::make_random_input(31, 4234);
    //df::InputData in = df::make_random_input(34, 4234);
    //df::InputData in = df::make_random_input(45, 4234);



    // these are valid inputs where the algorithm fails:
    //df::InputData in = df::make_random_input(31, 1312); // valid input
    //df::InputData in = daf::make_random_input(37, 4458); // upper and lower intersect in more than boundary edges
    //df::InputData in = df::make_random_input(35, 218); // valid input but non-decomposable
    //df::InputData in = df::make_random_input(40, 40);  // upper and lower intersect in more than boundary edges
    //df::InputData in = df::make_random_input(50, 40); // valid input but non-decomposable
    //df::InputData in = df::make_random_input(55,42);  // valid input but non-decomposable


  
  

    viz::register_triangulation_as_mesh(in.tri_lower, in.points2d, "lower 2D", "lower lifted");
    viz::register_triangulation_as_mesh(in.tri_upper, in.points2d, "upper 2D", "upper lifted");
    viz::register_regular_triangulation_as_mesh(in.tri_regular, in.points2d_weighted, "regular 2D", "regular lifted");

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


        std::vector<df::vertex_id> insertion_vertex_list = missing; // unsorted version

        // sort missing vertices by descending lifted height
        //std::vector<df::vertex_id> insertion_vertex_list = df::sorted_insertion_vertices(missing, in);
        df::vertex_id insertion_vertex = 0; // will only be used if we find a conforming one
        bool found_conforming = false;

        
        // try candidates in order until one is BOTH a downflip and conforming
        for (auto id : insertion_vertex_list) {

            // 1) check downflip condition
            if (!df::is_insertion_downflip(id, in, in.tri_current)) {
                std::cout << "[main] skipping vertex " << id
                        << " (insertion is not a down-flip)\n";
                continue;
            }

            // 2) check global conformance w.r.t. lower triangulation
            if (df::reg::is_insertion_conforming(id, in, in.tri_current)) {
                insertion_vertex = id;
                found_conforming = true;
                break;
            } else {
                std::cout << "[main] WARNING: insertion of vertex " << id
                        << " is non-conforming, trying next candidate\n";
            }
        }

        // if no conforming insertion exists -> polyhedron is non-decomposable
        if (!found_conforming) {
            std::cout << "\n[main] ERROR: all candidate vertex insertions are non-conforming.\n"
                    << "[main] The polyhedron appears to be non-decomposable.\n";
            break;  // break out of the main while(true) loop
        }
        
        
        // apply the insertion
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

    /*
    df::apply_edge_flip(5, 0, in);
    df::apply_edge_flip(4, 2, in);
    df::apply_edge_flip(3, 1, in);
    */

    //df::debug_print_local_to_global_map(in, df::TriKind::Lower);
    //df::debug_print_local_to_global_map(in, df::TriKind::Current);


    df::print_step_history(in);

    // build poset from upper to lower triangulation
    pst::build_poset(in);

    std::vector<df::DebugTetrahedron> debug_tets = df::collect_debug_tetrahedra(in);

    // visualize debug tetrahedra
    viz::load_debug_tetrahedra(in, debug_tets);

    // visualize flip decomposition
    viz::init_flip_decomposition(in);


    // initialize replay data
    df::init_replay(in);
   
    polyscope::state::userCallback = combined_ui_callback;


    polyscope::show();
    return 0;
}






















#include "decomp_algo.h"

#include <iostream>
#include <vector>

#include "edge_flip_loop.h"
#include "conforming_insertion.h"
#include "insertion.h"
#include "geometry_utils.h"  
#include "conforming.h"
#include "debug.h"
#include "flip.h"
#include "visualization.h"




namespace df {

bool run_flip_algorithm(df::InputData& in, bool enable_vis) {

    while (true) {

        // perform all conforming down-flips
        df::perform_all_conforming_down_flips(in);

        if (enable_vis)
            viz::show_or_update_current(in);

        // always show debug info
        df::debug_print_edge_list(in);

        // find missing vertices
        auto missing = df::find_missing_vertices(in.tri_current, in.tri_lower);

        if (missing.empty()) {
            std::cout << "\n=== no missing vertices left ===\n";
            break;
        }

        std::cout << "\nMissing vertices (global ids): ";
        for (auto id : missing) std::cout << id << " ";
        std::cout << "\n";

        df::vertex_id insertion_vertex = 0;
        bool found_conforming = false;

        // try candidates
        for (auto id : missing) {

            if (!df::is_insertion_downflip(id, in, in.tri_current)) {
                std::cout << "[flip_algo] skipping vertex " << id
                          << " (not a down-flip)\n";
                continue;
            }

            if (df::reg::is_insertion_conforming(id, in, in.tri_current)) {
                insertion_vertex = id;
                found_conforming = true;
                break;
            } else {
                std::cout << "[flip_algo] WARNING: insertion of vertex " << id
                          << " is non-conforming\n";
            }
        }

        if (!found_conforming) {
            std::cout << "\n[flip_algo] ERROR: no conforming insertion exists.\n";
            break;
        }

        std::cout << "Inserting vertex (global id) " << insertion_vertex << "\n";
        df::apply_vertex_insertion(insertion_vertex, in);

        if (enable_vis)
            viz::show_or_update_current(in);

        //df::debug_print_edge_list(in);
    }

    const bool ok = df::triangulations_equal(in.tri_current, in.tri_lower);
    if (!ok) {
        std::cout << "\n[flip_algo] ERROR: final triangulation differs from lower!\n";
    } else {
        std::cout << "\n[flip_algo] SUCCESS: current triangulation matches lower.\n";
    }

    return ok;
}

} // namespace df

    // -----------------------------------------------------------------------------------------------------

    // these are valid inputs where the algorithm works:
    //df::InputData in = df::make_random_input(7, 495934895);
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
    //df::InputData in = df::make_random_input(8, 495934895);
    //df::InputData in = df::make_random_input(9, 495934895);
    //df::InputData in = df::make_random_input(10, 495934895);
    //df::InputData in = df::make_random_input(11, 495934895); 
    //df::InputData in = df::make_random_valid_input(15, 495934895);
    //df::InputData in = df::make_random_input(25, 44); 
    //df::InputData in = df::make_random_input(31, 1312); 
    //df::InputData in = daf::make_random_input(37, 4458); 
    //df::InputData in = df::make_random_input(35, 218); 
    //df::InputData in = df::make_random_input(40, 40);  
    //df::InputData in = df::make_random_input(50, 40); 
    //df::InputData in = df::make_random_input(55,42);  


    //these are inputs where the algorithm succeeds and the upper triangulation is the farthees point triangulation
    //df::InputData in = df::make_random_valid_input(8, 495934895);
    //df::InputData in = df::make_random_valid_input(9, 495934895);
    //df::InputData in = df::make_random_valid_input(10, 495934895);


    //these are inputs where the algorithm fails and the upper triangulation is the farthees point triangulation
    // lower triangulation has then to be a local maximum in the poset
    //df::InputData in = df::make_random_valid_input(37, 4458);
    //df::InputData in = df::make_random_valid_input(55, 42);
    //df::InputData in = df::make_random_valid_input(31, 1312);
    //df::InputData in = df::make_random_valid_input(15, 495934895);
    //df::InputData in = df::make_random_valid_input(11, 495934895);


    // THE EXAMPLE to show that decomposition alogrithm is not correct
    //df::InputData in = df::make_random_valid_input(8, 495934895);
    //df::InputData in = df::make_random_valid_input(7, 1516235435);
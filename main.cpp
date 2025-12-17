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
#include "vis_poset.h"



int main() {
    polyscope::init();

    // number of vertices in triangulation
    int n_points = 8;
    // random seed to start point generation
    unsigned seed0 = 495934895;
    df::InputData in = df::make_random_valid_input(n_points, seed0);

   
    /**
    df::apply_edge_flip(5, 0, in, in.tri_lower);
    df::apply_edge_flip(4, 2, in, in.tri_lower);
    df::apply_edge_flip(3, 1, in, in.tri_lower);
    */

    

    // clear steop history to start fresh
    in.step_history.clear();
    


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
    //df::InputData in = df::make_random_valid_input(40, 40);




    //these are inputs where the algorithm fails and the upper triangulation is the farthees point triangulation
    // lower triangulation has then to be a local maximum in the poset
    //df::InputData in = df::make_random_valid_input(37, 4458);
    //df::InputData in = df::make_random_valid_input(55, 42);
    //df::InputData in = df::make_random_valid_input(31, 1312);
    //df::InputData in = df::make_random_valid_input(15, 495934895);
    //df::InputData in = df::make_random_valid_input(11, 495934895);





  

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


        std::vector<df::vertex_id> insertion_vertex_list = missing; 
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
    if (df::triangulations_equal(in.tri_current, in.tri_lower) == false){
        std::cout << "\n[main] ERROR: after all flips and insertions, current triangulation differs from lower triangulation!\n";
        
    } else {
        std::cout << "\n[main] SUCCESS: current triangulation matches lower triangulation!\n";
    }

    // print all global indices of current vertices in current triangulation as a list
    std::cout << "Current triangulation vertex global ids: ";
    for (auto vit = in.tri_current.finite_vertices_begin(); vit != in.tri_current.finite_vertices_end(); ++vit) {
        std::cout << vit->info() << " ";
    }
    std::cout << "\n";

    // print all global indices of current vertices in lower triangulation as a list
    std::cout << "Lower triangulation vertex global ids: ";
    for (auto vit = in.tri_lower.finite_vertices_begin(); vit != in.tri_lower.finite_vertices_end(); ++vit) {
        std::cout << vit->info() << " ";
    }
    std::cout << "\n";


    //df::debug_print_local_to_global_map(in, df::TriKind::Lower);
    //df::debug_print_local_to_global_map(in, df::TriKind::Current);


    df::print_step_history(in);
    std::vector<pst::Node> poset_nodes;



    
    // build the whole down flip poset from upper triangulation
    pst::build_poset(in, poset_nodes);
    viz_poset::register_poset(in, poset_nodes);

    // find minimal nodes in the poset (no incoming down-flips)
    auto mins = pst::nodes_with_no_incoming_down_flips(poset_nodes);
    std::cout << "Minimal nodes: ";
    for (int u : mins) std::cout << u << " ";
    std::cout << "\n";

    // gives us indices of special triangulations in the poset
    auto idx = pst::find_special_triangulations_in_poset(in, poset_nodes);
    std::cout << "poset index of upper (tri_upper):  " << idx.upper  << "\n";
    std::cout << "poset index of current (tri_current): " << idx.current << "\n";
    std::cout << "poset index of lower (tri_lower):   " << idx.lower   << "\n";


    // finds path from upper to lower triangulation in the poset (not necesessarily conforming)
    bool ok = pst::exists_path_via_children(poset_nodes, idx.upper, idx.lower);
    std::cout << "Path upper -> lower exists? " << std::boolalpha << ok << "\n";
    
    std::vector<int> node_path;
    std::vector<df::StepRecord> step_path;

    // finds conforming path in the poset from upper to lower triangulation
    bool ok_2 = pst::find_conforming_down_path_in_global_poset(
        in, poset_nodes, idx.upper, idx.lower, node_path, step_path);

    std::cout << "conforming down path exists? " << std::boolalpha << ok_2 << "\n";
    if (ok_2) {
        std::cout << "global mesh indices: ";
        for (int u : node_path) std::cout << u << " ";
        std::cout << "\n";
    }

    // print the steps along the conforming path
    std::cout << "\nConforming down path:\n";
    for (std::size_t i = 0; i < step_path.size(); ++i) {
        std::cout << "  " << node_path[i]
                << " -> " << node_path[i+1]
                << " : ";
        df::print_step_record(step_path[i]);
        std::cout << "\n";
    }

    




    
    /* this is for building the local poset around a chosen root node

    // empty history = upper triangulation as root
    std::vector<df::StepRecord> empty_history;

    pst::build_poset_local_down_from_history(in, empty_history, 1, poset_nodes, 1000);

    viz_poset::register_poset(in, poset_nodes);

    pst::debug_print_local_poset_histories(
    poset_nodes,
    empty_history.size()  // center_history_len
    );
    */




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






















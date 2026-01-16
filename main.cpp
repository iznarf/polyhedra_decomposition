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
#include "compare.h"
#include "poset2.h"



int main() {
    polyscope::init();

    // number of vertices in triangulation
    int n_points = 7;
    // random seed to start point generation
    unsigned seed0 = 1516235435;
    df::InputData in = df::make_random_valid_input(n_points, seed0);

   
    /*
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


    // THE EXAMPLE to show that alogrithm is not correct
    //df::InputData in = df::make_random_valid_input(8, 495934895);
    //df::InputData in = df::make_random_valid_input(7, 1516235435);

    
    // example to check <=_2 relation
    //df::InputData in = df::make_random_valid_input(5, 495934895);
    //df::InputData in = df::make_random_valid_input(6, 495934895);

    // example to show that unique meet/join do not exist
    //df::InputData in = df::make_random_valid_input(7, 1516235435);


    //viz::register_triangulation_as_mesh(in.tri_lower, in.points2d, "lower 2D", "lower lifted");
    //viz::register_triangulation_as_mesh(in.tri_upper, in.points2d, "upper 2D", "upper lifted");
    //viz::register_regular_triangulation_as_mesh(in.tri_regular, in.points2d_weighted, "regular 2D", "regular lifted");

    //viz::show_or_update_current(in);


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

        //viz::show_or_update_current(in);
        //df::debug_print_edge_list(in);

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
    int down_edge_count = viz_poset::register_poset(in, poset_nodes);

    // find minimal nodes in the poset (no incoming down-flips)
    auto mins = pst::nodes_with_no_incoming_down_flips(poset_nodes);
    std::cout << "minimal nodes: ";
    for (int u : mins) std::cout << u << " ";
    std::cout << "\n";

    // gives us indices of special triangulations in the poset
    auto idx = pst::find_special_triangulations_in_poset(in, poset_nodes);
    std::cout << "poset index of upper:  " << idx.upper  << "\n";
    std::cout << "poset index of current: " << idx.current << "\n";
    std::cout << "poset index of lower:   " << idx.lower   << "\n";


    // finds path from upper to lower triangulation in the poset (not necesessarily conforming)
    bool ok = pst::exists_path_via_children(poset_nodes, idx.upper, idx.lower);
    std::cout << "path upper -> lower exists? " << std::boolalpha << ok << "\n";
    
    std::vector<int> node_path;
    std::vector<df::StepRecord> step_path;

    // finds conforming path in the poset from upper to lower triangulation
    bool ok_2 = pst::find_conforming_down_path_in_global_poset(
        in, poset_nodes, idx.upper, idx.lower, node_path, step_path);

    std::cout << "conforming down path exists? " << std::boolalpha << ok_2 << "\n";

    /*
    if (ok_2) {
        std::cout << "global mesh indices: ";
        for (int u : node_path) std::cout << u << " ";
        std::cout << "\n";
    }
    */

    // print the steps along the conforming path
    std::cout << "\nconforming down path:\n";
    for (std::size_t i = 0; i < step_path.size(); ++i) {
        std::cout << "  " << node_path[i]
                << " -> " << node_path[i+1]
                << " : ";
        df::print_step_record(step_path[i]);
        std::cout << "\n";
    }

    // debug: test the compare function on all <=_1 edges in the poset and compare nodes which have no direct <=_1 relation

    if (poset_nodes.size() <= 0) {
        std::cout << "\n=== comparator test on <=1 edges ===\n";
        int fail = 0;
        int total = 0;

        for (int u = 0; u < (int)poset_nodes.size(); ++u) {
            for (int v : poset_nodes[u].children) {
                ++total;
                if (!pst2::compare(v, u, in, poset_nodes)) {
                    ++fail;
                    std::cout << "FAIL: v<=1u but compare(v,u)=false: "
                    << v << " <=1 " << u
                    << "   (edge stored as " << u << " -> " << v << ")\n";
                }
            }
        }

        std::cout << "checked " << total << " <=1 edges, failures = " << fail << "\n";

        std::cout << "\n=== full compare table (<= 0 nodes) ===\n";
        for (int u = 0; u < (int)poset_nodes.size(); ++u) {
            for (int v = u + 1; v < (int)poset_nodes.size(); ++v) {
                bool uv = pst2::compare(u, v, in, poset_nodes);
                bool vu = pst2::compare(v, u, in, poset_nodes);

                if (uv && vu) {
                    std::cout << u << " == " << v << "\n";
                } else if (uv) {
                    std::cout << u << " <2 " << v << "\n";
                } else if (vu) {
                    std::cout << v << " <2 " << u << "\n";
                } else {
                    std::cout << u << " || " << v << "\n";
                }
            }
        }
    } else {
        std::cout << "poset too large for full table: " << poset_nodes.size() << "\n";
    }

    pst2::Poset2 poset_2 = pst2::build_poset2(in, poset_nodes, true);

    int poset2_cover_edge_count = pst2::print_cover_relations(poset_2);

    viz_poset::register_poset2_cover_edges(poset_2.cover_out);

    //print number of cover edges for <=_1 poset
    std::cout << "[main] poset <=1 cover edges: " << down_edge_count << "\n";
    //print number of cover edges for <=_2 poset
    std::cout << "[main] poset <=2 cover edges: " << poset2_cover_edge_count << "\n";


    //----------------------------------------------------------------
    // test interval_xy function

    /*
    std::vector<int> interval = interval_xy(poset_2, 49, 0);
    if (interval.empty()) {
        std::cout << "[main] interval is empty\n";
    }
    else {
        for (int v : interval) std::cout << v << " ";
        std::cout << "\n";
    }
    */
    //----------------------------------------------------------------
   


    
    //-----------------------------------------------------------------
    // visualization of flip alogrithm

    //std::vector<df::DebugTetrahedron> debug_tets = df::collect_debug_tetrahedra(in);

    // visualize debug tetrahedra
    //viz::load_debug_tetrahedra(in, debug_tets);

    // visualize flip decomposition
    //viz::init_flip_decomposition(in);


    // initialize replay data
    //df::init_replay(in);

    //-----------------------------------------------------------------
   
    polyscope::state::userCallback = combined_ui_callback;


    polyscope::show();
    return 0;
}






















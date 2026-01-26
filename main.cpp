#include <polyscope/polyscope.h>

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <cstdint>
#include <sstream>
#include <GLFW/glfw3.h>

#include "input.h"
#include "decomp_algo.h"
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
#include "dedekind_cut.h"


static void glfw_error_silencer(int error, const char* description) {
    // do nothing
}

int main() {
    polyscope::options::verbosity = 0; // reduces Polyscope std::cout spam

    
    polyscope::init();
    glfwSetErrorCallback(glfw_error_silencer);

    // INPUT DATA GENERATION ------------------------------------------------------------

    // number of vertices in triangulation
    int n_points = 7;
    // random seed to start point generation
    unsigned seed0 = 1516235435;
    df::InputData in = df::make_random_valid_input(n_points, seed0);

    // DECOMPOSITION ALGORITHM -----------------------------------------------------------

    // enable visualization for decomposition algorithm
    bool enable_vis = false;
    if (enable_vis) {
        // setup visualization for the flip algorithm
        viz::setup_algo_visualization(in);
    }

    // clear step history to start fresh
    in.step_history.clear();

    // FLIP ALGORITHM LOOP

    // run flip algorithm with visualization
    df::run_flip_algorithm(in, enable_vis);

    // print step history
    df::print_step_history(in);

    if (enable_vis) {
        // setup debug tools for the flip algorithm 
        viz::setup_algo_debug_tools(in);
    }

    // --------------------------------------------------------------------------------



    // POSET INPUT EXAMPLES -----------------------------------------------------------------

    // example to check <=_2 relation
    //df::InputData in = df::make_random_valid_input(5, 495934895);
    //df::InputData in = df::make_random_valid_input(6, 495934895);

    // example to show that unique meet/join do not exist
    //df::InputData in = df::make_random_valid_input(7, 1516235435);

    // --------------------------------------------------------------------------------------


    // --------------------------------------------------------------------------------------
    // POSET1 COMPUTATION

    std::vector<pst::Node> poset_nodes;
    
    // build first the down flip poset, consiting of node list with children
    pst::build_poset(in, poset_nodes);

    // for later use we want the cover relations, cover_up edges and cover_down edges
    // so build Poset1 structure from node list
    pst::Poset1 P1 = pst::build_poset1(poset_nodes);

    // register poset for visualization
    viz_poset::register_poset(in, P1);

    // bool for debugging poset1
    bool debug_poset1 = false; 
    if (debug_poset1) {
        std::cout << "\n=== poset1 debug info ===\n";
    
        // find minimal nodes in the poset (no incoming down-flips)
        // just prints info if there exist more than one minimal node
        pst::nodes_with_no_incoming_down_edges(P1);

        // gives us indices of special triangulations in the poset
        auto idx = pst::find_special_triangulations_in_poset(in, P1);

        // finds path from upper to lower triangulation in the poset 
        pst::exists_path_via_children(P1, idx.upper, idx.lower);

        // checks if there exists a conforming down path from upper to lower triangulation in poset1
        std::vector<int> node_path;
        std::vector<df::StepRecord> step_path;

        // finds conforming path in the poset from upper to lower triangulation
        bool conf_path_exists = pst::find_conforming_down_path_in_global_poset(in, poset_nodes, idx.upper, idx.lower, node_path, step_path);

        if (conf_path_exists) {
            std::cout << "[poset] found conforming down path from upper triangulation "
                    << idx.upper << " to lower triangulation " << idx.lower << "\n";
            std::cout << "\nconforming down path:\n";
            for (std::size_t i = 0; i < step_path.size(); ++i) {
                std::cout << "  " << node_path[i]
                        << " -> " << node_path[i+1]
                        << " : ";
                df::print_step_record(step_path[i]);
                std::cout << "\n";
            }   
        }
    }

    // -----------------------------------------------------------------------------------------------------
    
    // POSET2 COMPUTATION

    // builds the poset2 structure from the input data and poset1
    pst2::Poset2 P2 = pst2::build_poset2(in, P1);

    bool debug_poset2 = false;
    if (debug_poset2) {
        // debug results of poset1 using comparator (geometric checks)
        pst2::debug_compare_poset1(in, P1);
    }

    // register poset2 edge network for visualization
    viz_poset::register_poset2_cover_edges(P2);

    
   

    polyscope::state::userCallback = combined_ui_callback;

    polyscope::show();
    return 0;
}






















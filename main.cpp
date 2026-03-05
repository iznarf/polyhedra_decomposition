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
#include "compare.h"
#include "poset2.h"
#include "dedekind_cut.h"
#include "regularity_check.h"

#include "input_ear_star.h"
#include "ear_star_vis.h"
#include "poset_ear_star.h"
#include "poset_ear_star_vis.h"
#include "ear_star_tri_word.h"


static void glfw_error_silencer(int error, const char* description) {
    // do nothing
}

pst::Poset1  g_P1;
pst::Poset_just_flips g_P1_flips;
pst2::Poset2 g_P2;
bool g_has_P1 = false;
bool g_has_P1_flips = false;
bool g_has_P2 = false;


df::InputData g_in;   // global input for UI + visualization



int main() {
    polyscope::options::verbosity = 0; // reduces polyscope errors in terminal 

    
    polyscope::init();
    glfwSetErrorCallback(glfw_error_silencer);

    // INPUT DATA GENERATION ------------------------------------------------------------

    // number of vertices in triangulation
    int n_points = 6;
    // random seed to start point generation
    unsigned seed0 = 1516235435;
    g_in = df::make_random_valid_input(n_points, seed0);
    df::InputData& in = g_in; // local ref for easier access

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

    g_P1 = pst::build_poset1(poset_nodes);
    g_has_P1 = true;

    
    // check regularity of all nodes in P1 by calling M2 via WSL and print summary table

    //regcheck::check_poset_regularity_wsl(in, g_P1, 1e6, "wsl", true);


    // register poset for visualization
    //viz_poset::register_poset(in, P1);

    // bool for debugging poset1
    bool debug_poset1 = false; 
    if (debug_poset1) {
        std::cout << "\n=== poset1 debug info ===\n";
    
        // find minimal nodes in the poset (no incoming down-flips)
        // just prints info if there exist more than one minimal node
        pst::nodes_with_no_incoming_down_edges(g_P1);

        // gives us indices of special triangulations in the poset
        auto idx = pst::find_special_triangulations_in_poset(in, g_P1);

        // finds path from upper to lower triangulation in the poset 
        pst::exists_path_via_children(g_P1, idx.upper, idx.lower);

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
    
    // JUST FLIP POSET 


    std::vector<pst::Node> poset_nodes_flips;
    pst::build_poset_just_flips(in, poset_nodes_flips);

    g_P1_flips = pst::build_poset_just_flips1(std::move(poset_nodes_flips));
    g_has_P1_flips = true;

    // --------------------------------------------------------------------------------------
    // POSET2 COMPUTATION

    // builds the poset2 structure from the input data and poset1

    g_P2 = pst2::build_poset2(in, g_P1);
    g_has_P2 = true;

    bool debug_poset2 = false;
    if (debug_poset2) {
        // debug results of poset1 using comparator (geometric checks)
        pst2::debug_compare_poset1(in, g_P1);
    }

    // register poset2 edge network for visualization
    //viz_poset::register_poset2_cover_edges(P2);

    // ------------------------------------------------------------------------------------------------------

    /*
    // EAR STAR TRIANGULATION -> WORD TEST

    // make input struct for ear star: convex n-gon with one point inside
    
    int n_boundary = 5; 

    InputData_ear_star in_ear_star = make_input_ear_star(n_boundary, 1.0);

    // build points2d array indexed by global id (0..n plus star=n)
    df::vertex_id max_id = 0;
    for (auto const& lp : in_ear_star.points) {
        max_id = std::max(max_id, lp.id);
    }

    std::vector<df::P2> points2d(max_id + 1);
    for (auto const& lp : in_ear_star.points) {
        points2d.at(lp.id) = lp.p;
    }

    // visualize triangulation
    register_triangulation_as_mesh_ear_star(in_ear_star.tri_start, points2d, "triangulation");


    // build ear star poset

    df::vertex_id star_id = in_ear_star.points.back().id;
    df::P2 star_point     = in_ear_star.points.back().p;

    // Build flip poset
    pst_es::FlipPoset P =
        pst_es::build_flip_poset(
            in_ear_star.tri_start,
            star_id,
            star_point
        );

    std::cout << "-----------------------------------\n";
    std::cout << "Number of triangulations found: "
            << P.nodes.size() << "\n";
    std::cout << "-----------------------------------\n";

    // Print all signatures
    for (std::size_t i = 0; i < P.nodes.size(); ++i) {
        std::cout << "Triangulation " << i << ":\n";

        for (auto const& e : P.nodes[i].sig.edges) {
            std::cout << "  (" << e[0] << "," << e[1] << ")\n";
        }

        std::cout << "\n";
    }

    // poset ear star visualization: register every triangulation in the poset as a separate mesh in polyscope, laid out on a grid by translating the vertex positions
    pst_es_viz::register_poset_triangulations_grid(
        in_ear_star.tri_start,
        P,
        points2d,
        star_id,
        star_point,
        "earstar",
        2.0  // spacing 
    );

    earstar::print_words_for_poset(
        in_ear_star.tri_start,
        P,
        n_boundary,
        star_id,
        star_point
    );

    */

    // -----------------------------------------------------------------------------------------------------
    polyscope::state::userCallback = combined_ui_callback;

    polyscope::show();
    return 0;
}






















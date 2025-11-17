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


int main() {
    polyscope::init();

    //df::InputData in = df::make_random_input(21, 48);
    df::InputData in = df::make_random_input(19, 40);

    viz::show_four_meshes(in);
    viz::show_or_update_current(in);

    {
        std::vector<int> locals = {6, 4, 18, 17}; // the polyscope vertex indices you see
        viz::debug_print_lower_vertex_ids(in, locals);
    }



    while (true) {

        // perform all conforming down-flips for the current triangulation until flip list is empty
        df::perform_all_conforming_down_flips(in);

        // see which vertices of the lower triangulation are still missing
        auto missing = df::find_missing_vertices(in.tri_current, in.tri_lower);

        // if vertex set of current triangulation matches that of lower triangulation, we are almost done
        // fix: compare CGAL triangulations. same vertex set does not guarantee same triangulation!
        if (missing.empty()) { 
            std::cout << "\n=== Finished: no missing vertices and no down-flips left ===\n";
            break; 
        }

        // print missing vertex ids in current triangulation
        std::cout << "\nMissing vertices (global ids): ";
        for (auto id : missing) std::cout << id << " ";
        std::cout << "\n";

        // pick candidates sorted by height (highest first)
        std::vector<df::vertex_id> insertion_vertex_list = df::sorted_insertion_vertices(missing, in);  // returns list of global ids sorted by height 

        df::vertex_id insertion_vertex = 0; // will only be used if we find a conforming one
        bool found_conforming = false;

        // try candidates in order until one is conforming
        for (auto id : insertion_vertex_list) {
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
                    << "[main] The polyhedron appears to be non-decomposable by conforming flips.\n";
            break;  // break out of the main while(true) loop
        }

        std::cout << "Inserting vertex (global id) " << insertion_vertex << "\n";
        df::apply_vertex_insertion(insertion_vertex, in);



        viz::show_or_update_current(in);
        viz::debug_print_edge_list(in);

        }
    
    std::vector<int> locals = {17, 5, 12, 0};
    viz::debug_print_lower_vertex_ids(in, locals);

    df::print_flip_history(in);


    polyscope::show();
    return 0;
}




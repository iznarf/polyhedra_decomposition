#include <polyscope/polyscope.h>
#include <CGAL/Simple_cartesian.h>

#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <cstdint>
#include <sstream>

typedef CGAL::Simple_cartesian<double> Kernel;

#include "input.h"
#include "visualization.h"
#include "check_edges.h"
#include "conforming.h"
#include "flip.h"


int main() {
    polyscope::init();

    // build input
    df::InputData in = df::make_random_input(21, 48);

    // triangulations 
    const df::Tri2& source_triangulation  = in.tri_upper;
    const df::Tri2& target_triangulation  = in.tri_lower;
    df::Tri2&       current_triangulation = in.tri_current;

    // visualize (keeps Polyscope indices consistent)
    viz::show_four_meshes(in);

    // global vertex_id -> local polyscope index of current triangulation
    std::unordered_map<df::vertex_id, int> id2local_current;
    id2local_current.reserve(current_triangulation.number_of_vertices());
    int loc = 0;
    for (auto v = current_triangulation.finite_vertices_begin();
        v != current_triangulation.finite_vertices_end(); ++v) {
        id2local_current[v->info()] = loc++;
    }

    // global vertex_id -> local polyscope index of target triangulation
    std::unordered_map<df::vertex_id, int> id2local_lower;
    id2local_lower.reserve(target_triangulation.number_of_vertices());
    loc = 0;
    for (auto v = target_triangulation.finite_vertices_begin();
        v != target_triangulation.finite_vertices_end(); ++v) {
        id2local_lower[v->info()] = loc++;
    }

    // dump the mapping once
    #ifdef DF_DEBUG_CONFORMING
    std::cerr << "[main] id2local map (global -> local):\n";
    for (const auto& kv : id2local_current) {
        std::cerr << "  " << kv.first << " -> " << kv.second << "\n";
    }
    #endif

    
    // collect locally non-regular (down-flippable) edges from the current triangulation
    std::vector<std::array<df::vertex_id,2>> down_flip_edges = df::reg::find_locally_non_regular_edges(current_triangulation);

    // for storing flippable edges which are NOT in the target triangulation and conforming 

    std::unordered_set<std::pair<df::vertex_id, df::vertex_id>, boost::hash<std::pair<df::vertex_id, df::vertex_id>>> flippable_edges;
    flippable_edges.reserve(down_flip_edges.size() * 2 + 1);

    std::cout << "\n=== Locally non-regular (down-flippable) edges ===\n";
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

        bool conf = df::reg::is_flip_conforming(ia, ib, in, id2local_current, id2local_lower);

        // in_target is false and conf is true -> we can flip this edge -> add it to the set 
        if (!in_target && conf) {
            flippable_edges.insert({ia, ib});
        }

        // visualize the flip tetrahedra 
        viz::show_flip_tetra(in, ia, ib, std::to_string(k));

        // print one line per edge
        std::cout << std::setw(3) << k << " (" << ia << "," << ib << ")->"
                << std::setw(8) << pl << "  "
                << (in_target ? "T" : "F") << "         "
                << (conf ? "T" : "F") << "\n";
    }

    // now apply an edge flip to the first element in the list of flippable edges
    std::cout << "\n=== Applying conforming flip ===\n";

    if (!flippable_edges.empty()) {
        const auto [ia, ib] = *flippable_edges.begin();
        // we have to print the local polyscope indices here
        int local_index_a = id2local_current.at(ia);
        int local_index_b = id2local_current.at(ib);
        std::cout << "flipping edge (" << local_index_a << "," << local_index_b << ")\n";
        df::apply_edge_flip(ia, ib, in);
        // update the visualization
    }

  polyscope::show();
  return 0;
}



#pragma once
#include "input.h"
#include <vector>
#include <unordered_map>

namespace viz {

// Registers 4 meshes 
void show_four_meshes(const df::InputData& D);

// Collect global vertex IDs in the same order as used for Polyscope registration
std::vector<df::vertex_id> present_ids(const df::Tri2& t);

// Map global vertex ID -> local compact index [0..V-1] in Polyscope arrays
std::unordered_map<df::vertex_id,int> make_local_index(const std::vector<df::vertex_id>& ids);

// visualizes one tetrahedron
void show_flip_tetra(const df::InputData& D, df::vertex_id ia, df::vertex_id ib, const std::string& label);

void show_or_update_current(const df::InputData& D);

void debug_print_lower_vertex_ids(const df::InputData& D,
                                    const std::vector<int>& local_indices);


void debug_print_edge_list(const df::InputData& D);

} // namespace viz







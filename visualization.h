#pragma once
#include "input.h"
#include "debug.h"
#include <vector>
#include <unordered_map>

namespace viz {

// registers a triangulation as lifted and planar meshes in polyscope
void register_triangulation_as_mesh(const df::Tri2& tri,
                                    const std::vector<df::P2>& points2d,
                                    const std::string& name_planar,
                                    const std::string& name_lifted);


// Collect global vertex IDs in the same order as used for Polyscope registration
std::vector<df::vertex_id> present_ids(const df::Tri2& t);

// Map global vertex ID -> local compact index [0..V-1] in Polyscope arrays
std::unordered_map<df::vertex_id,int> make_local_index(const std::vector<df::vertex_id>& ids);

void show_or_update_current(const df::InputData& D);

// this is for replaying the recorded steps, needs fixing 
void show_or_update_replay(const df::InputData& D);

// load a set of debug tetrahedra (indices + type) and
void load_debug_tetrahedra(const df::InputData& D,
                           const std::vector<df::DebugTetrahedron>& tets);

// buttons / slider to step through the debug tets
void debug_tet_ui();




} // namespace viz







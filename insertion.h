#pragma once

#include <vector>
#include <array>
#include <cstdint>

#include "input.h"   

namespace df {


    // collect vertex ids that are in target but not in current
    std::vector<df::vertex_id>
    find_missing_vertices(const Tri2& current, const Tri2& target);

    // pick next vertex to be inserted from the missing set
    df::vertex_id
    pick_next_insertion_vertex(const std::vector<df::vertex_id>& missing, const df::InputData& D);

    // apply vertex insertion into the current triangulation
    void apply_vertex_insertion(const df::vertex_id id, df::InputData& D);

} // namespace df

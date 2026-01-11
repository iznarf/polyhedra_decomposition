#pragma once

#include "input.h"
#include "poset.h"
#include "poset2.h"
#include <vector>

namespace viz_poset {

    
    // draws all 2D + 3D meshes for the given poset nodes
    int register_poset(const df::InputData& D, const std::vector<pst::Node>& nodes);

    // draws the ImGui controls for the poset
    void poset_ui();

    // registers / update a second edge network that visualizes <=2 cover edges
    // cover_out convention: u -> v means u <=2 v (bottom->top)
    // we reverse for drawing (top->down)
    void register_poset2_cover_edges(const std::vector<std::vector<int>>& cover_out);

  
    const std::vector<pst::Node>& get_poset1_nodes();
    const pst2::Poset2& get_poset2();
    

} // namespace viz_poset

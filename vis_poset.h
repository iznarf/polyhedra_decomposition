#pragma once

#include "input.h"
#include "poset.h"
#include <vector>

namespace viz_poset {

    // Call this once after you have built the poset and have the nodes.
    // It stores pointers / copies internally and registers all meshes in Polyscope.
    void register_poset(const df::InputData& D,
                        const std::vector<pst::Node>& nodes);

    // Draws the ImGui controls (checkboxes, sliders, etc.) for the poset.
    void poset_ui();

    // Register / update a second edge network that visualizes <=2 cover edges
    // cover_out convention: u -> v means u <=2 v (bottom->top)
    // we reverse for drawing (top->down)
    void register_poset2_cover_edges(const std::vector<std::vector<int>>& cover_out);


} // namespace viz_poset

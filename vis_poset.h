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

} // namespace viz_poset

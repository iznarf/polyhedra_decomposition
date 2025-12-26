#pragma once

#include "input.h"
#include "poset.h"
#include "poset2.h"
#include <vector>

namespace viz_poset2 {

    // register all node meshes and the <=2 cover edge network 
    // P2.cover_out: <=2 covering edges: u->v means u <=2 v 
    void register_poset2(const df::InputData& D,
                         const std::vector<pst::Node>& nodes,
                         const pst2::Poset2& P2);

    // ImGui controls for toggling visibility
    void poset2_ui();

} 

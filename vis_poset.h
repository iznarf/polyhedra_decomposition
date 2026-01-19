#pragma once

#include "input.h"
#include "poset.h"
#include "poset2.h"
#include <vector>
#include <glm/glm.hpp>


namespace viz_poset {

    
    // draws all 2D + 3D meshes for the given poset nodes
    void register_poset(const df::InputData& D, const pst::Poset1& P1);

    void register_poset1_node_meshes(const df::InputData& D, const pst::Poset1& P1);
    void register_poset1_edges(const pst::Poset1& P1);


    // draws the ImGui controls for the poset
    void poset_ui();

    // registers / update a second edge network that visualizes <=2 cover edges
    // cover_out convention: u -> v means u <=2 v (bottom->top)
    // we reverse for drawing (top->down)
    void register_poset2_cover_edges(const pst2::Poset2& P2);

    void rebuild_downflip_network_filtered();

    void set_poset_enabled(bool enabled);
    void apply_poset_visibility_from_masks();
    bool node_visible(int idx);



    const pst2::Poset2& get_poset2();
    const pst::Poset1& get_poset1();
    struct PosetLayout {
        std::vector<int> level;
        std::vector<std::vector<int>> byLevel;
        std::vector<glm::vec3> centers;
        int max_level = 0;
    };

    enum class ActiveLayout { P1, P2 };

    extern PosetLayout g_layout_p1;
    extern PosetLayout g_layout_p2;
    extern ActiveLayout g_active_layout;





    

} // namespace viz_poset

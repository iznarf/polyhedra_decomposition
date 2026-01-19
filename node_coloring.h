#pragma once
#include <glm/glm.hpp>

namespace viz_poset {

// base coloring helpers
void reset_node_coloring();
void color_node(int idx, glm::vec3 c2, glm::vec3 c3);

// Möbius coloring modes (new API)
void color_nodes_by_mobius_anchor_min_to_all_poset1();
void color_nodes_by_mobius_anchor_min_to_all_poset2();
void color_nodes_by_mobius_all_to_anchor_max_poset1();
void color_nodes_by_mobius_all_to_anchor_max_poset2();

} // namespace viz_poset


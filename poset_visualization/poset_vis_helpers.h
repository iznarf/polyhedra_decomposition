   
#pragma once

#include "poset_utils.h"

#include <vector>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <queue>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <array>
#include <glm/glm.hpp>

namespace pst_vis{

    std::vector<int> compute_levels_longest_from_roots_cover_up(const std::vector<std::vector<int>>& cover_up);

    std::vector<std::vector<int>>group_nodes_by_level(const std::vector<int>& levels);

    std::vector<glm::vec3> grid_for_poset_nodes(const std::vector<std::vector<int>>& byLevel,int number_of_nodes,glm::vec3 center,float xSpacing,float ySpacing);
      
}
    
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




#include "render_helpers.h"   
#include "input.h"            

namespace pst_vis{

    
    // levels = longest distance from ANY root in cover_up
    // cover_up[u] contains immediate parents of u (edge u -> parent)
    // roots are nodes with indegree 0 in cover_up
    // throws if graph has a directed cycle (via topo_sort_kahn)
    std::vector<int> compute_levels_longest_from_roots_cover_up(const std::vector<std::vector<int>>& cover_up){
        const int N = (int)cover_up.size();
        std::vector<int> level(N, 0);
        if (N == 0) return level;

        const std::vector<int> topo = pst::topo_sort_kahn(cover_up);

        // indegree in cover_up
        std::vector<int> indeg(N, 0);
        for (int u = 0; u < N; ++u) {
            for (int p : cover_up[u]) {
                if (p < 0 || p >= N) throw std::runtime_error("levels_cover_up: endpoint out of range");
                indeg[p]++;
            }
        }

        const int NEG_INF = std::numeric_limits<int>::min() / 4;
        std::vector<int> dist(N, NEG_INF);

        // initialize all roots (indeg==0) with distance 0
        for (int i = 0; i < N; ++i) {
            if (indeg[i] == 0) dist[i] = 0;
        }

        // DP for longest paths along cover_up
        for (int u : topo) {
            if (dist[u] == NEG_INF) continue;
            for (int p : cover_up[u]) {
                dist[p] = std::max(dist[p], dist[u] + 1);
            }
        }

        for (int i = 0; i < N; ++i) {
            level[i] = (dist[i] == NEG_INF) ? 0 : dist[i];
        }
        return level;
    }


    // groups nodes by their level:
    // result[l] = list of node indices i with levels[i] == l
    std::vector<std::vector<int>>group_nodes_by_level(const std::vector<int>& levels){
        int maxL = -1;
        for (int lv : levels) maxL = std::max(maxL, lv);

        std::vector<std::vector<int>> buckets((size_t)maxL + 1);
        for (int i = 0; i < (int)levels.size(); ++i) {
            int lv = levels[i];
            if (lv >= 0) buckets[(size_t)lv].push_back(i); // ignore -1 if you use that
        }
        return buckets;
    }



    // byLevel[l] = node indices at level l
    std::vector<glm::vec3> grid_for_poset_nodes(const std::vector<std::vector<int>>& byLevel,
                        int number_of_nodes,
                        glm::vec3 center,
                        float xSpacing,
                        float ySpacing)
    {
        std::vector<glm::vec3> pos((size_t)number_of_nodes, center);

        const int L = (int)byLevel.size();
        if (L == 0) return pos;

        const float totalH = (L - 1) * ySpacing;
        const float zTop   = center.z + 0.5f * totalH;

        for (int l = 0; l < L; ++l) {
            const auto& nodes = byLevel[(size_t)l];
            const int k = (int)nodes.size();

            const float totalW = (k > 0 ? (k - 1) * xSpacing : 0.f);
            const float xLeft  = -0.5f * totalW;

            const float z = zTop - l * ySpacing;

            for (int j = 0; j < k; ++j) {
                int n = nodes[(size_t)j];
                float x = xLeft + j * xSpacing;

                pos[(size_t)n] = glm::vec3(
                    center.x + x,
                    center.y,
                    center.z + z
                );
            }
        }
        return pos;
    }

}
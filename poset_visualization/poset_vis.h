#pragma once

#include "input.h"
#include "poset.h"
#include "poset2.h"
#include "poset_tri_mesh.h"
#include "dedekind_poset.h"

#include <glm/glm.hpp>
#include <string>
#include <vector>

namespace pst_vis {

// for global triangulation mesh visulization
extern TriVisResult g_triVis_P1; // for P1
extern TriVisResult g_triVis_P2;   // for P2
extern TriVisResult g_triVis_P1_flips; // for P1 just flips visualization

// for global grid positions, we need this for interval overlay registration
extern std::vector<glm::vec3> g_gridPos_P1;
extern std::vector<glm::vec3> g_gridPos_P2;
extern std::vector<glm::vec3> g_gridPos_P1_flips;



extern float g_gridXSpacing_P1;
extern float g_gridZSpacing_P1;
extern float g_gridXSpacing_P2;
extern float g_gridZSpacing_P2;
extern float g_gridXSpacing_P1_flips;
extern float g_gridZSpacing_P1_flips;


// function to show/hide all trivis meshes in a TriVisResult
// example: show all P1 2D meshes 
void set_trivis_enabled(const pst_vis::TriVisResult& R, bool show2d, bool show3d);



void register_poset_as_grid(const std::string& name,
                            int number_of_nodes,
                            const std::vector<int>& levels,
                            glm::vec3 center,
                            float xSpacing,
                            float zSpacing,
                            glm::vec3 color);
                            

// register poset nodes as a colored grid (per-node colors), returns computed positions
std::vector<glm::vec3> register_poset_as_grid_colored(
    const std::string& name,
    int number_of_nodes,
    const std::vector<int>& levels,
    const glm::vec3& center,
    float xSpacing,
    float zSpacing,
    const std::vector<glm::vec3>& colors,
    float pointRadius = 0.01f
);


void register_cover_up_as_edges(const std::string& name, const std::vector<glm::vec3>& pos, const std::vector<std::vector<int>>& cover_up, glm::vec3 color);

// register P1 poset: triangulation, grid, edges
void register_poset1(const df::InputData& D, const pst::Poset1& P1,const glm::vec3& center, float xSpacing, float zSpacing, glm::vec3 color);

void register_poset1_flips(const df::InputData& D,
                             const pst::Poset_just_flips& P,
                             const glm::vec3& center,
                             float xSpacing, float zSpacing,
                             glm::vec3 color);

// register P2 poset: triangulation, grid, edges
void register_poset2(const df::InputData& D, const pst::Poset1& P1, const pst2::Poset2& P2,const glm::vec3& center, float xSpacing, float zSpacing, glm::vec3 color);




} // namespace pst_vis


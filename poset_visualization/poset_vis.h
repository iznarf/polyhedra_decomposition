#pragma once

#include "input.h"
#include "poset.h"
#include "poset2.h"
#include "poset_tri_mesh.h"
#include <glm/glm.hpp>
#include <string>
#include <vector>

namespace pst_vis {

// for global triangulation mesh visulization
extern TriVisResult g_triVis_P1; // for P1
extern TriVisResult g_triVis_P2;   // for P2

// for global grid positions, we need this for interval overlay registration
extern std::vector<glm::vec3> g_gridPos_P1;
extern std::vector<glm::vec3> g_gridPos_P2;

// function to show/hide all trivis meshes in a TriVisResult
// example: show all P1 2D meshes 
void set_trivis_enabled(const pst_vis::TriVisResult& R, bool show2d, bool show3d);

// register P1 poset: triangulation, grid, edges
void register_poset1(const df::InputData& D, const pst::Poset1& P1,const glm::vec3& center, float xSpacing, float zSpacing, glm::vec3 color);

// register P2 poset: triangulation, grid, edges
void register_poset2(const df::InputData& D, const pst::Poset1& P1, const pst2::Poset2& P2,const glm::vec3& center, float xSpacing, float zSpacing, glm::vec3 color);

} // namespace pst_vis


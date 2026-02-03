#pragma once

#include "input.h"   
#include "poset.h"   

#include <polyscope/surface_mesh.h>
#include <glm/vec3.hpp>
#include <string>
#include <vector>

namespace pst_vis {

struct TriVisResult {
    std::vector<polyscope::SurfaceMesh*> meshes2d;
    std::vector<polyscope::SurfaceMesh*> meshes3d;
    std::vector<glm::vec3> centers;

    std::vector<glm::vec3> pivot2;
    std::vector<glm::vec3> pivot3;

};

void remove_triangulation_vis(TriVisResult& R);

void register_poset_as_triangulation_P1(const std::string& name_prefix,
                                        const df::InputData& D,
                                        const pst::Poset1& P1,
                                        TriVisResult& out);
void apply_triangulation_centers(TriVisResult& R,
                                 const std::vector<glm::vec3>& centers,
                                 float uniformScale);



} // namespace pst_vis



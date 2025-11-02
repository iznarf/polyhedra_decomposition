#pragma once
#include <string>
#include <vector>
#include <array>
#include "input.h"

namespace viz {

// Initialize Polyscope (call once)
// Optional title is used to set polyscope's programName before initialization.
void init();

// Register (or re-register) the planar triangulation at z=0.
// It extracts vertices/faces from the given CGAL triangulation.
void register_planar_triangulation(const df::Tri2& tri,
                                   const std::string& name = "planar");

// Register (or re-register) the lifted triangulation at z=ω(x,y).
void register_lifted_triangulation(const df::Tri2& tri,
                                   const df::OmegaQuad& omega,
                                   const std::string& name = "lifted");

// Update the planar mesh faces/vertices from the current triangulation.
// (Safe to call even if the mesh wasn't registered yet; it will register.)
void update_planar_triangulation(const df::Tri2& tri,
                                 const std::string& name = "planar");

// Update the lifted mesh faces/vertices from the current triangulation + ω.
void update_lifted_triangulation(const df::Tri2& tri,
                                 const df::OmegaQuad& omega,
                                 const std::string& name = "lifted");

// Launch Polyscope UI (blocking).
void show();

} // namespace viz

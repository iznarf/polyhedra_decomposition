#pragma once



#include <array>
#include <vector>


namespace env_max {

// ---------------------------
// output mesh
// ---------------------------
struct Mesh {
  // vertices as double (x,y,z)
  std::vector<std::array<double, 3>> V;
  // triangles as indices into V
  std::vector<std::array<int, 3>> F;
  // per triangle: which input surface induced it
  // encoded as (mesh_id, face_id)
  std::vector<std::array<int, 2>> tri_tag;
};

// ---------------------------
// Input triangle with a tag
// ---------------------------
struct TriTag {
  int mesh_id = 0; // e.g. 1 for A, 2 for B
  int face_id = -1;
};

struct InputTriangle {
  // vertices in double; internally we convert to exact
  std::array<double, 3> p0;
  std::array<double, 3> p1;
  std::array<double, 3> p2;
  TriTag tag;
};

// Compute the upper envelope mesh of the given triangles.
// Returns true on success.
bool compute_lower_envelope(const std::vector<InputTriangle>& tris, Mesh& out);




} // namespace env_max

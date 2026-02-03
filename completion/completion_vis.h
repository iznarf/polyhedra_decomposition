#pragma once

#include <glm/glm.hpp>
#include <string>

namespace df { struct InputData; }
namespace pst { struct Poset1; }
namespace pst2 { struct Poset2; }

namespace dm_completion { struct DedekindPoset; }

namespace completion_vis {

// Build grid + edges + 2D overlay (triangulations for old nodes, diagrams for new nodes)
// which: 0 = P1 completion, 1 = P2 completion
void build(int which,
           const df::InputData& D,
           const pst::Poset1& P1,
           const pst2::Poset2& P2,
           const dm_completion::DedekindPoset& C,
           const glm::vec3& center,
           float xSpacing,
           float zSpacing);

// Show/hide groups for a completion visualization
void set_enabled(int which, bool showGrid, bool showEdges, bool show2D);

// Clear only this completion visualization
void clear(int which);

// Clear both completions
void clear_all();

} // namespace completion_vis

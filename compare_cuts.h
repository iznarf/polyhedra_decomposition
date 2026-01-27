#pragma once

#include "dedekind_poset.h"

namespace viz_poset {

// UI: overlay meshes of min(F) and max(I') for two Dedekind cuts
// Uses base poset triangulations (poset1 histories) like compare_nodes.
void compare_cuts_ui(const pst2::DedekindPoset& D);

} // namespace viz_poset

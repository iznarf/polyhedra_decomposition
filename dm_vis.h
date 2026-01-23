#pragma once

#include "poset2.h"
#include "dedekind_cut.h"
#include "poset2.h"

#include <vector>

namespace viz_dm {

// Build DM lattice graph (nodes = cuts, edges = cover relations) and register with Polyscope.
// Nodes are laid out by rank = |I'| (y-axis), spread in x within each rank.
void build_and_register_dm_graph(const pst2::Poset2& P2, const std::vector<pst2::DedekindCut>& cuts);

// show/hide DM graph
void set_dm_enabled(bool enabled);
bool dm_enabled();

// selection (driven from ImGui slider, not Polyscope picking)
void set_selected_dm_node(int idx);
int  selected_dm_node();

// read selected cut (for UI printing)
const pst2::DedekindCut* selected_cut();

} // namespace viz_dm


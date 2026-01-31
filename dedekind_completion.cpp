#include "dedekind_completion.h"

#include "poset_view.h"
#include "dedekind_cut.h"
#include "dedekind_poset.h"

#include <imgui.h>

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <glm/glm.hpp>

extern pst::Poset1  g_P1; // global poset1
extern pst2::Poset2 g_P2; // global poset2

// ------------------------------------------------------------
// Store the latest completion result
// ------------------------------------------------------------
static dm_completion::DedekindPoset g_completion;
static bool g_hasCompletion = false;
static std::string g_status;

// UI state: which input poset to complete (0=P1, 1=P2)
static int  g_selectedPoset = 0;
static bool g_showCompletion = false;

// ------------------------------------------------------------
// Build a PosetView for the completion
// ------------------------------------------------------------
static PosetView view_of_completion(const dm_completion::DedekindPoset& D) {
    PosetView V;
    V.n = (int)D.cuts.size();

    V.cover_up          = &D.cover_up;
    V.cover_down        = &D.cover_down;
    V.topo_up           = &D.topo_up;
    V.topo_down         = &D.topo_down;
    V.reachability_up   = &D.reachability_up;
    V.reachability_down = &D.reachability_down;
    V.levels            = &D.levels;

    return V;
}

// ------------------------------------------------------------
// Choose input poset
// ------------------------------------------------------------
static PosetView current_input_view() {
    if (g_selectedPoset == 0) return view_of(g_P1);
    return view_of(g_P2);
}

static std::string completion_nodes_name() {
    return (g_selectedPoset == 0) ? "P1 completion nodes" : "P2 completion nodes";
}

static std::string completion_edges_name() {
    return (g_selectedPoset == 0) ? "P1 completion edges" : "P2 completion edges";
}




// ------------------------------------------------------------
// UI
// ------------------------------------------------------------
void dedekind_completion_ui() {

    ImGui::SeparatorText("Dedekind completion");

    ImGui::TextUnformatted("Input poset:");
    ImGui::RadioButton("P1", &g_selectedPoset, 0);
    ImGui::SameLine();
    ImGui::RadioButton("P2", &g_selectedPoset, 1);

    ImGui::Spacing();

    if (ImGui::Button("Compute completion")) {

        g_status.clear();

        // 1) input poset as PosetView
        PosetView P = current_input_view();

        // 2) compute all cuts
        std::vector<dm_completion::DedekindCut> cuts = dm_completion::compute_all_cuts(P);

        // 3) build completion
        g_completion = dm_completion::build_dedekind_poset(P, std::move(cuts));
        g_hasCompletion = true;

        // 4) lattice check
        PosetView PC = view_of_completion(g_completion);
        bool ok = dm_completion::check_complete_lattice(PC, true);

        // 5) UI status
        g_status += "Completion nodes: " + std::to_string((int)g_completion.cuts.size()) + "\n";
        g_status += ok ? "OK: completion is a lattice (finite => complete lattice)\n"
                       : "FAIL: completion is NOT a lattice\n";
        g_status += "Max level: " + std::to_string(g_completion.maxLevel) + "\n";

    }
}


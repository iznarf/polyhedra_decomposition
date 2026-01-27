#include "input.h"
#include "replay.h"
#include "vis_poset.h"
#include "visualization.h"
#include "moebius.h"
#include "poset.h"
#include "compare_nodes.h"
#include "interval_ui.h"
#include "meet_join_ui.h"
#include "dedekind_cut_ui.h"
#include "dm_vis.h"

#include <imgui.h>

void combined_ui_callback() {

    // ------------------------------------------------------------
    // Flip algorithm (grouped)
    // ------------------------------------------------------------
    if (ImGui::CollapsingHeader("FLIP ALGORITHM")) {

        ImGui::PushID("ReplayUI");
        df::replay_ui();
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("DebugTetUI");
        viz::debug_tet_ui();
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("DecompositionUI");
        viz::flip_decomposition_ui();
        ImGui::PopID();
    }

    ImGui::Separator();

    // ------------------------------------------------------------
    // Poset + friends
    // ------------------------------------------------------------
     if (ImGui::CollapsingHeader("POSET")) {
        ImGui::PushID("Poset");
        viz_poset::poset_ui();
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("IntervalUI");
        viz_poset::interval_ui();
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("MeetJoinUI");
        viz_poset::meet_join_ui();
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("CompareNodes");
        viz_poset::compare_nodes_ui();
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("Moebius");
        mob::draw_mobius_compare_ui(viz_poset::get_poset1(), viz_poset::get_poset2());
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("Completion");
        viz_poset::dedekind_cut_ui();
        ImGui::PopID();

        ImGui::Separator();

        ImGui::PushID("dmVis");
        viz_dm::dm_ui();
        ImGui::PopID();
    }
}




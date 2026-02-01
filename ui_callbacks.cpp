#include "input.h"
#include "replay.h"
#include "poset_vis.h"
#include "visualization.h"
#include "poset.h"
#include "poset_vis_ui.h"
#include "interval_ui.h"
#include "meet_join_ui.h"
#include "moebius_ui.h"
#include "compare_nodes.h"
#include "dedekind_cut_ui.h"
#include "dedekind_completion.h"


#include <imgui.h>

void combined_ui_callback() {

    // ------------------------------------------------------------
    // Flip algorithm 
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
        pst_vis_ui::poset_ui();
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

        ImGui::PushID("Moebius");
        viz_poset::moebius_ui();
        ImGui::PopID();
        

        ImGui::Separator();

        ImGui::PushID("CompareNodes");
        viz_poset::compare_nodes_ui();
        ImGui::PopID();

        
        ImGui::Separator();

        ImGui::PushID("DedekindCut");
        viz_poset::dedekind_cut_ui();
        ImGui::PopID();


        ImGui::Separator();

        ImGui::PushID("DedekindCompletion");
        dedekind_completion_ui();
        ImGui::PopID();
        
    }
}




#include "input.h"
#include "replay.h"
#include "vis_poset.h"
#include "visualization.h"
#include "moebius.h"
#include "poset.h"
#include <imgui.h>

void combined_ui_callback() {

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

    ImGui::Separator();

    ImGui::PushID("Poset");
    viz_poset::poset_ui();
    ImGui::PopID();

    ImGui::Separator();

    
    ImGui::PushID("Moebius");
    mob::draw_mobius_compare_ui(
    viz_poset::get_poset1_nodes(),
    viz_poset::get_poset2()
    );

    ImGui::PopID();
    
}



#include "input.h"
#include "replay.h"
#include "vis_poset.h"
#include "visualization.h"
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


    ImGui::PushID("FlipPoset");
    viz_poset::poset_ui();
    ImGui::PopID();


}



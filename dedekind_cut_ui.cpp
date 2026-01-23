#include "dedekind_cut_ui.h"
#include "dedekind_cut.h"
#include "dm_vis.h"

#include "node_coloring.h"
#include "dedekind_cut.h"
#include "vis_poset.h"

#include <imgui.h>
#include <glm/glm.hpp>
#include <cstdlib>
#include <vector>

namespace viz_poset {

static int parse_int(const char* s) { return (s && *s) ? std::atoi(s) : 0; }

void dedekind_cut_ui() {

    static bool have_all = false;
    static std::vector<pst2::DedekindCut> allCuts;
    static int cut_idx = 0;


    if (!ImGui::CollapsingHeader("dedekind cut"))
        return;

    // Legend
    ImGui::BulletText("I  ideal of {x,y} : green");
    ImGui::BulletText("F  filter : black");
    ImGui::BulletText("I' : pink");
    ImGui::BulletText("{x,y} generators : yellow");
    ImGui::Separator();

        // --- persistent UI state ---
    static char bufA[32] = "0";
    static char bufB[32] = "0";
    static int show_mode = 2;          // 0=I, 1=F, 2=(I',F)
    static bool have_cut = false;

    static int lastA = 0, lastB = 0;
    static std::vector<int> cachedI, cachedF, cachedIp;

    // helper: recompute
    auto recompute = [&]() {
        int nodeA = parse_int(bufA);
        int nodeB = parse_int(bufB);

        lastA = nodeA; lastB = nodeB;

        const auto& P2 = get_poset2();
        std::vector<int> gens = { nodeA, nodeB };

        cachedI  = pst2::ideal_generated_by(P2, gens);
        cachedF  = pst2::filter_of_generators(P2, gens);
        cachedIp = pst2::Iprime_from_filter(P2, cachedF);

        have_cut = true;
    };

    // --- UI ---
    ImGui::InputText("x", bufA, sizeof(bufA));
    ImGui::InputText("y", bufB, sizeof(bufB));

    ImGui::RadioButton("show I", &show_mode, 0); ImGui::SameLine();
    ImGui::RadioButton("show F", &show_mode, 1); ImGui::SameLine();
    ImGui::RadioButton("show (I',F)", &show_mode, 2); 


    if (ImGui::Button("compute")) {
        recompute();
    }

    ImGui::SameLine();
    if (ImGui::Button("clear")) {
        have_cut = false;
        reset_node_coloring();
    }

    // --- APPLY COLORING EVERY FRAME  ---
    if (have_cut) {
        reset_node_coloring();

        // colors
        glm::vec3 green2 (0.30f, 1.00f, 0.40f);
        glm::vec3 green3 (0.10f, 0.80f, 0.20f);

        glm::vec3 pink2  (1.00f, 0.40f, 0.80f);
        glm::vec3 pink3  (1.00f, 0.20f, 0.70f);

        glm::vec3 black2 (0.05f, 0.05f, 0.05f);
        glm::vec3 black3 (0.00f, 0.00f, 0.00f);

        glm::vec3 yellow2(1.00f, 0.90f, 0.10f);
        glm::vec3 yellow3(1.00f, 0.70f, 0.05f);

        auto color_list = [&](const std::vector<int>& nodes, glm::vec3 c2, glm::vec3 c3) {
            for (int v : nodes) color_node(v, c2, c3);
        };

        if (show_mode == 0) {
            color_list(cachedI, green2, green3);
        } else if (show_mode == 1) {
            color_list(cachedF, black2, black3);
        } else if (show_mode == 2) {
            color_list(cachedIp, pink2, pink3);
            color_list(cachedF,  black2, black3);
        }

        // highlight generators always
        color_node(lastA, yellow2, yellow3);
        color_node(lastB, yellow2, yellow3);

        ImGui::Text("Sizes: |I|=%d  |F|=%d  |I'|=%d",
                    (int)cachedI.size(), (int)cachedF.size(), (int)cachedIp.size());
    }

    ImGui::Separator();

    if (ImGui::CollapsingHeader("DM completion", 0)) {

    static bool dm_graph_built = false;
    static int  dm_node        = 0;
    static bool color_base_from_dm = true;

    if (ImGui::Button("compute all cuts")) {
        const auto& P2 = get_poset2();
        allCuts = pst2::compute_all_cuts(P2);
        have_all = true;

        viz_dm::build_and_register_dm_graph(P2, allCuts);
        dm_graph_built = true;

        dm_node = 0;
        viz_dm::set_selected_dm_node(dm_node);
        viz_dm::set_dm_enabled(true);
    }

    if (have_all) {
        ImGui::Text("Total cuts: %d", (int)allCuts.size());

        // Toggle button (requested)
        if (dm_graph_built) {
            if (ImGui::Button(viz_dm::dm_enabled() ? "Hide DM graph" : "Show DM graph")) {
                viz_dm::set_dm_enabled(!viz_dm::dm_enabled());
            }
        }

        ImGui::Checkbox("color base poset from selected DM cut", &color_base_from_dm);

        ImGui::SliderInt("DM node", &dm_node, 0, (int)allCuts.size() - 1);
        if (ImGui::Button("Prev") && dm_node > 0) dm_node--;
        ImGui::SameLine();
        if (ImGui::Button("Next") && dm_node + 1 < (int)allCuts.size()) dm_node++;

        viz_dm::set_selected_dm_node(dm_node);

        const pst2::DedekindCut* Cp = viz_dm::selected_cut();
        if (Cp) {
            const auto& C = *Cp;

            ImGui::Text("Sizes: |F|=%d |I'|=%d  |minF|=%d |maxI'|=%d",
                (int)C.F.size(), (int)C.Iprime.size(),
                (int)C.minimal_elements_F.size(),
                (int)C.maximal_elements_Iprime.size());

            ImGui::TextUnformatted("minF:");
            for (int x : C.minimal_elements_F) { ImGui::SameLine(); ImGui::Text("%d", x); }

            ImGui::TextUnformatted("maxI':");
            for (int x : C.maximal_elements_Iprime) { ImGui::SameLine(); ImGui::Text("%d", x); }

            if (color_base_from_dm) {
                reset_node_coloring();

                glm::vec3 pink2  (1.00f, 0.40f, 0.80f);
                glm::vec3 pink3  (1.00f, 0.20f, 0.70f);
                glm::vec3 black2 (0.05f, 0.05f, 0.05f);
                glm::vec3 black3 (0.00f, 0.00f, 0.00f);

                for (int v : C.Iprime) color_node(v, pink2, pink3);
                for (int v : C.F)      color_node(v, black2, black3);
            }
        }
    }
}
}


} // namespace viz_poset

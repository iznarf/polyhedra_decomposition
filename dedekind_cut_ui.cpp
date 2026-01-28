#include "dedekind_cut_ui.h"
#include "dedekind_cut.h"
#include "dm_vis.h"
#include "dedekind_poset.h"
#include "poset2.h"

#include "node_coloring.h"
#include "vis_poset.h"
#include "compare_nodes.h"
#include "compare_cuts.h"

#include <imgui.h>
#include <glm/glm.hpp>
#include <cstdlib>
#include <vector>
#include <algorithm> 

namespace viz_poset {

static int parse_int(const char* s) { return (s && *s) ? std::atoi(s) : 0; }

void dedekind_cut_ui() {


    static std::vector<pst2::DedekindCut> allCuts;
    static pst2::DedekindPoset g_DM;  
    static bool have_all =false;



    if (ImGui::CollapsingHeader("dedekind single cut", 0)) {

        // ------------------------------------------------------------
        // Cut visualizer (x,y -> generated cut)
        // ------------------------------------------------------------
        ImGui::BulletText("I  ideal of {x,y} : green");
        ImGui::BulletText("F  filter : black");
        ImGui::BulletText("I' : pink");
        ImGui::BulletText("{x,y} generators : yellow");
        ImGui::Separator();

        static char bufA[32] = "0";
        static char bufB[32] = "0";
        static int  show_mode = 2;          // 0=I, 1=F, 2=(I',F)
        static bool have_cut  = false;

        static int lastA = 0, lastB = 0;
        static std::vector<int> cachedI, cachedF, cachedIp;

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

        if (have_cut) {
            reset_node_coloring();

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
            } else {
                color_list(cachedIp, pink2, pink3);
                color_list(cachedF,  black2, black3);
            }

            color_node(lastA, yellow2, yellow3);
            color_node(lastB, yellow2, yellow3);

            ImGui::Text("Sizes: |I|=%d  |F|=%d  |I'|=%d", (int)cachedI.size(), (int)cachedF.size(), (int)cachedIp.size());
    
            
        }
    }
    ImGui::Separator();

    // ------------------------------------------------------------
    // DM lattice (compute + show/hide)
    // ------------------------------------------------------------

    if (ImGui::CollapsingHeader("dedekind completion in poset2", 0)) {
        static bool dm_graph_built = false;
        static int  dm_node        = 0;
        static bool color_base_from_dm = true;

        

            if (ImGui::Button("compute completion")) {
                const auto& P2 = get_poset2();

                // compute once
                allCuts = pst2::compute_all_cuts(P2);

                // build DM poset (do NOT redeclare g_DM here!)
                g_DM = pst2::build_dedekind_poset(P2, allCuts);
                pst2::check_complete_lattice(g_DM, true);
                have_all = !g_DM.cuts.empty();

                viz_dm::build_and_register_dm_graph(g_DM);
                dm_graph_built = true;

                dm_node = 0;
                viz_dm::set_selected_dm_node(dm_node);
                viz_dm::set_dm_enabled(true);
            }


            if (!have_all) {
                ImGui::TextUnformatted("(no cuts computed yet)");
            } else {
                ImGui::Text("Total cuts: %d", (int)g_DM.cuts.size());
                ImGui::SliderInt("dm cut", &dm_node, 0, (int)g_DM.cuts.size() - 1);
                viz_dm::set_selected_dm_node(dm_node);


                if (dm_graph_built) {
                    if (ImGui::Button(viz_dm::dm_enabled() ? "hide dm lattice" : "show dm lattice")) {
                        viz_dm::set_dm_enabled(!viz_dm::dm_enabled());
                    }
                }

                // turn on/off coloring from dm selection
                if (ImGui::Checkbox("color poset from selected cut", &color_base_from_dm)) {
                    if (!color_base_from_dm) {
                        reset_node_coloring();
                        // optional, if you want to restore your default poset visibility rules:
                        // viz_poset::apply_poset_visibility_from_masks();
                    }
                }


                if (const pst2::DedekindCut* Cp = viz_dm::selected_cut()) {
                const auto& C = *Cp;

                ImGui::Text("sizes: |F|=%d |I'|=%d  |minF|=%d |maxI'|=%d",
                    (int)C.F.size(), (int)C.Iprime.size(),
                    (int)C.minimal_elements_F.size(),
                    (int)C.maximal_elements_Iprime.size());

                // --- show elements (like you asked) ---
                ImGui::TextUnformatted("min(F):");
                for (int x : C.minimal_elements_F) { ImGui::SameLine(); ImGui::Text("%d", x); }

                ImGui::TextUnformatted("max(I'):");
                for (int x : C.maximal_elements_Iprime) { ImGui::SameLine(); ImGui::Text("%d", x); }

                ImGui::TextUnformatted("max(I'): white, min(F): yellow, I': pink, F: black");

                // --- coloring (match the other section's style) ---
                if (color_base_from_dm) {
                    reset_node_coloring();

                    // same palette as your generator cut UI
                    glm::vec3 pink2  (1.00f, 0.40f, 0.80f);
                    glm::vec3 pink3  (1.00f, 0.20f, 0.70f);

                    glm::vec3 black2 (0.05f, 0.05f, 0.05f);
                    glm::vec3 black3 (0.00f, 0.00f, 0.00f);

                    glm::vec3 yellow2(1.00f, 0.90f, 0.10f);
                    glm::vec3 yellow3(1.00f, 0.70f, 0.05f);

                    glm::vec3 white2(1.00f, 1.00f, 1.00f);
                    glm::vec3 white3(1.00f, 1.00f, 1.00f);



                    // base sets
                    for (int v : C.Iprime) color_node(v, pink2, pink3);
                    for (int v : C.F)      color_node(v, black2, black3);

                    // highlight extrema (like "generators" highlight)
                    for (int v : C.minimal_elements_F)        color_node(v, yellow2, yellow3);
                    for (int v : C.maximal_elements_Iprime)  color_node(v, white2, white3);
                }
            }
        }
    }
    ImGui::Separator();
    viz_poset::compare_cuts_ui(g_DM);

}


   

} // namespace viz_poset



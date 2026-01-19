#include "interval_ui.h"

#include "poset2.h"
#include "vis_poset.h"

#include <imgui.h>
#include <cstdlib>
#include <iostream>
#include <vector>


#include "interval_ui.h"

#include "poset2.h"
#include "vis_poset.h"

#include <imgui.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <utility>   
#include <cstdio>    

extern int  g_interval_x;
extern int  g_interval_y;
extern bool g_interval_active;
extern std::vector<char> g_interval_mask;
extern char g_interval_x_buf[32];
extern char g_interval_y_buf[32];

extern bool g_focus_active;
extern std::vector<char> g_focus_mask;

namespace viz_poset {

// defined in vis_poset.cpp
void apply_poset_visibility_from_masks();

static void apply_node_list_as_interval_mask(const std::vector<int>& nodes) {
    const auto& P1 = viz_poset::get_poset1();
    const int n = (int)P1.nodes.size();

    g_interval_mask.assign(n, 0);
    for (int v : nodes) {
        if (0 <= v && v < n) g_interval_mask[v] = 1;
    }
    g_interval_active = true;

    viz_poset::apply_poset_visibility_from_masks();
}

void interval_ui() {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("poset2 interval visualization using cover_up", ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }

    const auto& P2 = viz_poset::get_poset2();
    const auto& P1 = viz_poset::get_poset1();

    const int N = (int)P2.cover_up.size();
    if (N == 0 || (int)P2.cover_down.size() != N) {
        ImGui::TextUnformatted("poset2 not loaded / invalid.");
        return;
    }

    // --- input ---
    ImGui::InputText("x", g_interval_x_buf, IM_ARRAYSIZE(g_interval_x_buf));
    ImGui::InputText("y", g_interval_y_buf, IM_ARRAYSIZE(g_interval_y_buf));

    g_interval_x = std::atoi(g_interval_x_buf);
    g_interval_y = std::atoi(g_interval_y_buf);

    auto valid_idx = [&](int v) { return 0 <= v && v < N; };
    bool x_ok = valid_idx(g_interval_x);
    bool y_ok = valid_idx(g_interval_y);

    if (!x_ok || !y_ok) {
        ImGui::Text("Valid indices: 0 .. %d", N - 1);
    }

    // --- options ---
    static bool clear_focus_on_apply = true;
    ImGui::Checkbox("clear focus when applying", &clear_focus_on_apply);

    // --- quick actions row ---
    if (ImGui::Button("swap x/y")) {
        std::swap(g_interval_x, g_interval_y);
        std::snprintf(g_interval_x_buf, IM_ARRAYSIZE(g_interval_x_buf), "%d", g_interval_x);
        std::snprintf(g_interval_y_buf, IM_ARRAYSIZE(g_interval_y_buf), "%d", g_interval_y);
        x_ok = valid_idx(g_interval_x);
        y_ok = valid_idx(g_interval_y);
    }
    ImGui::SameLine();

    if (ImGui::Button("show full poset")) {
        g_interval_active = false;
        g_interval_mask.clear();
        viz_poset::apply_poset_visibility_from_masks();
    }

    // Optional: jump to unique min/max if they exist
    int min0 = pst2::unique_min_or_minus1(P2);
    int max1 = pst2::unique_max_or_minus1(P2);

    if (min0 != -1) {
        ImGui::SameLine();
        if (ImGui::Button("set x=min")) {
            g_interval_x = min0;
            std::snprintf(g_interval_x_buf, IM_ARRAYSIZE(g_interval_x_buf), "%d", g_interval_x);
            x_ok = valid_idx(g_interval_x);
        }
    }
    if (max1 != -1) {
        ImGui::SameLine();
        if (ImGui::Button("set y=max")) {
            g_interval_y = max1;
            std::snprintf(g_interval_y_buf, IM_ARRAYSIZE(g_interval_y_buf), "%d", g_interval_y);
            y_ok = valid_idx(g_interval_y);
        }
    }

    ImGui::Separator();

    // --- main actions ---
    if (ImGui::Button("visualize interval [x,y]")) {
        if (!x_ok || !y_ok) {
            std::cout << "[vis_poset] invalid x/y\n";
        } else {
            if (clear_focus_on_apply) {
                g_focus_active = false;
                g_focus_mask.clear();
            }

            auto nodes = pst2::interval_xy(P2, g_interval_x, g_interval_y);

            // If empty, tell user why in a UI-friendly way.
            // We can’t distinguish “x not <= y” vs “empty for other reason”
            // without exposing upper_set; but in a poset that’s the reason.
            if (nodes.empty() && g_interval_x != g_interval_y) {
                std::cout << "[vis_poset] interval empty: likely x not <= y\n";
            }

            apply_node_list_as_interval_mask(nodes);

            std::cout << "[vis_poset] interval [" << g_interval_x << "," << g_interval_y
                      << "] size=" << nodes.size() << "\n";
        }
    }

    ImGui::SameLine();

    if (ImGui::Button("visualize down(x)")) {
        if (!x_ok) {
            std::cout << "[vis_poset] invalid x\n";
        } else {
            if (clear_focus_on_apply) { g_focus_active = false; g_focus_mask.clear(); }

            // If you decide to remove interval_min_x later, replace with
            // a new pst2::down_cone_nodes(P2,x) wrapper.
            auto nodes = pst2::interval_min_x(P2, g_interval_x);

            apply_node_list_as_interval_mask(nodes);

            std::cout << "[vis_poset] down(" << g_interval_x << ") size=" << nodes.size() << "\n";
        }
    }

    ImGui::SameLine();

    if (ImGui::Button("visualize up(x)")) {
        if (!x_ok) {
            std::cout << "[vis_poset] invalid x\n";
        } else {
            if (clear_focus_on_apply) { g_focus_active = false; g_focus_mask.clear(); }

            auto nodes = pst2::interval_x_max(P2, g_interval_x);

            apply_node_list_as_interval_mask(nodes);

            std::cout << "[vis_poset] up(" << g_interval_x << ") size=" << nodes.size() << "\n";
        }
    }

    // small status line
    ImGui::Separator();
    ImGui::Text("Interval mode: %s", g_interval_active ? "ON" : "OFF");
    ImGui::Text("poset2 nodes: %d   poset1 nodes: %d", N, (int)P1.nodes.size());
}

} // namespace viz_poset

#include "meet_join_ui.h"
#include "vis_poset.h"
#include "poset2.h"
#include "node_coloring.h"
#include "poset_utils.h"

#include <imgui.h>
#include <cstdlib>
#include <iostream>
#include <vector>

// globals from vis_poset.cpp
extern pst::Poset1 g_poset1;
extern pst2::Poset2 g_poset2;

extern char g_pair_a_buf[32];
extern char g_pair_b_buf[32];

extern bool g_focus_active;
extern std::vector<char> g_focus_mask;

extern bool g_interval_active;
extern std::vector<char> g_interval_mask;

namespace viz_poset {

void apply_poset_visibility_from_masks();

void meet_join_ui() {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("poset2 meet/join visualization", ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }

    const auto& P1 = viz_poset::get_poset1();
    const auto& P2 = viz_poset::get_poset2();

    const int n = (int)P1.nodes.size();
    if (n == 0 || (int)P2.cover_up.size() != n) {
        ImGui::TextUnformatted("poset not loaded / size mismatch.");
        return;
    }

    ImGui::InputText("x", g_pair_a_buf, IM_ARRAYSIZE(g_pair_a_buf));
    ImGui::InputText("y", g_pair_b_buf, IM_ARRAYSIZE(g_pair_b_buf));

    int a = std::atoi(g_pair_a_buf);
    int b = std::atoi(g_pair_b_buf);

    auto valid = [&](int v) { return 0 <= v && v < n; };

    if (!valid(a) || !valid(b)) {
        ImGui::Text("Valid indices: 0 .. %d", n - 1);
    }

    // small UX helpers
    if (ImGui::Button("swap a/b")) {
        std::swap(a, b);
        std::snprintf(g_pair_a_buf, IM_ARRAYSIZE(g_pair_a_buf), "%d", a);
        std::snprintf(g_pair_b_buf, IM_ARRAYSIZE(g_pair_b_buf), "%d", b);
    }

    ImGui::SameLine();
    static bool clear_interval_on_apply = true;
    ImGui::Checkbox("clear interval when applying", &clear_interval_on_apply);

    auto mark_path = [&](const std::vector<int>& path) {
        for (int v : path) if (0 <= v && v < n) g_focus_mask[v] = 1;
    };

    auto start_focus_view = [&]() {
        g_focus_active = true;
        viz_poset::apply_poset_visibility_from_masks();
    };

    if (ImGui::Button("visualize meet")) {
        if (!valid(a) || !valid(b)) {
            std::cout << "[vis_poset] meet: invalid a/b\n";
            return;
        }

        if (clear_interval_on_apply) {
            g_interval_active = false;
            g_interval_mask.clear();
        }

        viz_poset::reset_node_coloring();
        auto meets = pst2::meet_candidates_xy(P2, a, b);

        g_focus_mask.assign(n, 0);
        g_focus_mask[a] = 1;
        g_focus_mask[b] = 1;
        for (int m : meets) if (valid(m)) g_focus_mask[m] = 1;

        // meet m is below a and b => show paths m -> a and m -> b along UP edges
        for (int m : meets) {
            if (!valid(m)) continue;
            mark_path(pst::find_path_up_cover_up(P2.cover_up, m, a));
            mark_path(pst::find_path_up_cover_up(P2.cover_up, m, b));
        }

        start_focus_view();

        // coloring
        viz_poset::color_node(a, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));
        viz_poset::color_node(b, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));
        for (int m : meets) {
            if (!valid(m)) continue;
            viz_poset::color_node(m, glm::vec3(1.0f, 0.2f, 0.6f), glm::vec3(1.0f, 0.2f, 0.6f));
        }

        std::cout << "[vis_poset] meet candidates for (" << a << "," << b << "): " << meets.size() << "\n";
    }

    ImGui::SameLine();

    if (ImGui::Button("visualize join")) {
        if (!valid(a) || !valid(b)) {
            std::cout << "[vis_poset] join: invalid a/b\n";
            return;
        }

        if (clear_interval_on_apply) {
            g_interval_active = false;
            g_interval_mask.clear();
        }

        viz_poset::reset_node_coloring();
        auto joins = pst2::join_candidates_xy(P2, a, b);

        g_focus_mask.assign(n, 0);
        g_focus_mask[a] = 1;
        g_focus_mask[b] = 1;
        for (int j : joins) if (valid(j)) g_focus_mask[j] = 1;

        // show a -> j and b -> j along UP edges
        for (int j : joins) {
            if (!valid(j)) continue;
            mark_path(pst::find_path_up_cover_up(P2.cover_up, a, j));
            mark_path(pst::find_path_up_cover_up(P2.cover_up, b, j));
        }

        start_focus_view();

        // coloring
        viz_poset::color_node(a, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));
        viz_poset::color_node(b, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));
        for (int j : joins) {
            if (!valid(j)) continue;
            viz_poset::color_node(j, glm::vec3(0.2f, 1.0f, 0.4f), glm::vec3(0.2f, 1.0f, 0.4f));
        }

        std::cout << "[vis_poset] join candidates for (" << a << "," << b << "): " << joins.size() << "\n";
    }

    ImGui::SameLine();

    if (ImGui::Button("reset meet/join view")) {
        g_focus_active = false;
        g_focus_mask.clear();
        viz_poset::apply_poset_visibility_from_masks();
        viz_poset::reset_node_coloring();
    }
}

} // namespace viz_poset

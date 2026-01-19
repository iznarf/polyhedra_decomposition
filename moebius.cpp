#include "moebius.h"
#include "poset_utils.h"

#include <imgui.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cstdint>

namespace mob {
CoverUpView build_view(const std::vector<std::vector<int>>& cover_up) {
    CoverUpView V;
    V.cover_up = &cover_up;
    V.topo = pst::topo_sort_kahn(cover_up);
    V.R    = pst::compute_reachability(cover_up, V.topo);
    return V;
}

bool leq(const CoverUpView& V, int a, int b) {
    if (a == b) return true;
    if (!V.valid()) return false;
    const int N = (int)V.R.size();
    if (a < 0 || a >= N || b < 0 || b >= N) return false;
    return V.R[a].test(b);
}

}

namespace viz_poset {
    void reset_node_coloring();
    void color_nodes_by_mobius_anchor_min_to_all_poset1();
    void color_nodes_by_mobius_anchor_min_to_all_poset2();
    void color_nodes_by_mobius_all_to_anchor_max_poset1();
    void color_nodes_by_mobius_all_to_anchor_max_poset2();
}

namespace mob {

int unique_min_cover(const std::vector<std::vector<int>>& cover_down) {
    int N = (int)cover_down.size();
    int min = -1;
    for (int v = 0; v < N; ++v) {
        if (cover_down[v].empty()) {
            if (min != -1) return -1;
            min = v;
        }
    }
    return min;
}

int unique_max_cover(const std::vector<std::vector<int>>& cover_up) {
    int N = (int)cover_up.size();
    int max = -1;
    for (int v = 0; v < N; ++v) {
        if (cover_up[v].empty()) {
            if (max != -1) return -1;
            max = v;
        }
    }
    return max;
}


// ---------- helpers ----------
static std::vector<int> interval_xy_from_view(const mob::CoverUpView& V, int x, int y) {
    const int N = (int)V.R.size();
    if (x < 0 || x >= N || y < 0 || y >= N) return {};
    if (!mob::leq(V, x, y)) return {};

    std::vector<int> interval;
    interval.reserve(N);
    for (int z = 0; z < N; ++z) {
        if (mob::leq(V, x, z) && mob::leq(V, z, y)) interval.push_back(z);
    }
    return interval;
}

static std::int64_t mobius_xy_from_view(const mob::CoverUpView& V, int x, int y) {
    const int N = (int)V.R.size();
    if (N == 0) return 0;
    if (x < 0 || x >= N || y < 0 || y >= N) return 0;
    if (!mob::leq(V, x, y)) return 0;
    if (x == y) return 1;

    std::vector<int> I = interval_xy_from_view(V, x, y);
    if (I.empty()) return 0;

    std::vector<char> inI(N, 0);
    for (int z : I) inI[z] = 1;

    std::vector<int> topoI;
    topoI.reserve(I.size());
    for (int v : V.topo) if (inI[v]) topoI.push_back(v);

    std::vector<std::int64_t> mu(N, 0);
    mu[x] = 1;

    for (int v : topoI) {
        if (v == x) continue;

        std::int64_t sum = 0;
        for (int z : topoI) {
            if (z == v) break;            // earlier in topoI only
            if (!mob::leq(V, z, v)) continue;
            sum += mu[z];
        }
        mu[v] = -sum;
    }

    return mu[y];
}


// ---------- UI ----------
void draw_mobius_compare_ui(
    const pst::Poset1& P1, const pst2::Poset2& P2,
    const mob::CoverUpView& V1, const mob::CoverUpView& V2)
{
    static char x_buf[32] = "0";
    static char y_buf[32] = "0";

    enum class WhichPoset { P1, P2, Both };
    static WhichPoset which = WhichPoset::Both;

    static bool computed = false;
    static std::int64_t mu1 = 0;
    static std::int64_t mu2 = 0;

    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Möbius μ(x,y) using cover_up", ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }

    const int N1 = (int)P1.cover_up.size();
    const int N2 = (int)P2.cover_up.size();

    // inputs
    ImGui::InputText("x (lower)", x_buf, IM_ARRAYSIZE(x_buf));
    ImGui::InputText("y (upper)", y_buf, IM_ARRAYSIZE(y_buf));
    int x = std::atoi(x_buf);
    int y = std::atoi(y_buf);

    // quick: swap
    if (ImGui::Button("swap x/y")) {
        std::swap(x, y);
        std::snprintf(x_buf, IM_ARRAYSIZE(x_buf), "%d", x);
        std::snprintf(y_buf, IM_ARRAYSIZE(y_buf), "%d", y);
    }

    ImGui::Separator();

    // choose poset(s)
    ImGui::TextUnformatted("Compute in:");
    if (ImGui::RadioButton("P1", which == WhichPoset::P1)) which = WhichPoset::P1;
    ImGui::SameLine();
    if (ImGui::RadioButton("P2", which == WhichPoset::P2)) which = WhichPoset::P2;
    ImGui::SameLine();
    if (ImGui::RadioButton("Both", which == WhichPoset::Both)) which = WhichPoset::Both;

    // diagnostics
    auto valid_idx = [&](int v, int N) { return 0 <= v && v < N; };

    bool x_ok_1 = valid_idx(x, N1), y_ok_1 = valid_idx(y, N1);
    bool x_ok_2 = valid_idx(x, N2), y_ok_2 = valid_idx(y, N2);

    if (which != WhichPoset::P2) {
        ImGui::Text("P1 size: %d  | valid indices: 0..%d", N1, std::max(0, N1 - 1));
        if (x_ok_1 && y_ok_1) {
            ImGui::Text("P1: x <= y ? %s", mob::leq(V1, x, y) ? "YES" : "NO");
        }
    }
    if (which != WhichPoset::P1) {
        ImGui::Text("P2 size: %d  | valid indices: 0..%d", N2, std::max(0, N2 - 1));
        if (x_ok_2 && y_ok_2) {
            ImGui::Text("P2: x <= y ? %s", mob::leq(V2, x, y) ? "YES" : "NO");
        }
    }

    ImGui::Separator();

    if (ImGui::Button("compute μ(x,y)")) {
        computed = true;
        mu1 = 0; mu2 = 0;

        if (which == WhichPoset::P1 || which == WhichPoset::Both) {
            if (x_ok_1 && y_ok_1) {
                mu1 = mobius_xy_from_view(V1, x, y);
                std::cout << "[moebius] P1 mu(" << x << "," << y << ")=" << mu1 << "\n";
            } else {
                std::cout << "[moebius] P1 invalid x/y\n";
            }
        }
        if (which == WhichPoset::P2 || which == WhichPoset::Both) {
            if (x_ok_2 && y_ok_2) {
                mu2 = mobius_xy_from_view(V2, x, y);
                std::cout << "[moebius] P2 mu(" << x << "," << y << ")=" << mu2 << "\n";
            } else {
                std::cout << "[moebius] P2 invalid x/y\n";
            }
        }
    }

    ImGui::SameLine();

    if (ImGui::Button("print interval nodes")) {
        if (which == WhichPoset::P1 || which == WhichPoset::Both) {
            if (x_ok_1 && y_ok_1) {
                auto I = interval_xy_from_view(V1, x, y);
                std::cout << "[moebius] P1 interval [" << x << "," << y << "] size=" << I.size() << "\n";
                for (int z : I) std::cout << z << " ";
                std::cout << "\n";
            }
        }
        if (which == WhichPoset::P2 || which == WhichPoset::Both) {
            if (x_ok_2 && y_ok_2) {
                auto I = interval_xy_from_view(V2, x, y);
                std::cout << "[moebius] P2 interval [" << x << "," << y << "] size=" << I.size() << "\n";
                for (int z : I) std::cout << z << " ";
                std::cout << "\n";
            }
        }
    }

    if (computed) {
        if (which == WhichPoset::P1 || which == WhichPoset::Both)
            ImGui::Text("P1: μ(%d,%d) = %lld", x, y, (long long)mu1);
        if (which == WhichPoset::P2 || which == WhichPoset::Both)
            ImGui::Text("P2: μ(%d,%d) = %lld", x, y, (long long)mu2);
    }

    ImGui::Separator();

    // Coloring actions
    ImGui::TextUnformatted("Color nodes by Möbius (global):");
    if (ImGui::Button("P2: color μ(min, v)")) viz_poset::color_nodes_by_mobius_anchor_min_to_all_poset2();
    ImGui::SameLine();
    if (ImGui::Button("P2: color μ(v, max)")) viz_poset::color_nodes_by_mobius_all_to_anchor_max_poset2();

    if (ImGui::Button("P1: color μ(min, v)")) viz_poset::color_nodes_by_mobius_anchor_min_to_all_poset1();
    ImGui::SameLine();
    if (ImGui::Button("P1: color μ(v, max)")) viz_poset::color_nodes_by_mobius_all_to_anchor_max_poset1();

    if (ImGui::Button("reset node coloring")) viz_poset::reset_node_coloring();

    ImGui::TextUnformatted("Colors: black=0, red=+1, blue=-1 (fallback: sign for |μ|!=1)");
}


void mob::draw_mobius_compare_ui(const pst::Poset1& P1, const pst2::Poset2& P2) {
    static mob::CoverUpView V1;
    static mob::CoverUpView V2;
    static bool built = false;

    // rebuild only when sizes change (cheap and safe)
    if (!built || (int)V1.R.size() != (int)P1.cover_up.size() || (int)V2.R.size() != (int)P2.cover_up.size()) {
        V1 = mob::build_view(P1.cover_up);
        V2 = mob::build_view(P2.cover_up);
        built = true;
    }

    // call your 4-arg UI
    mob::draw_mobius_compare_ui(P1, P2, V1, V2);
}


} // namespace mob


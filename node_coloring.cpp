#include "node_coloring.h"
#include "vis_poset.h"
#include "moebius.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <glm/glm.hpp>
#include <iostream>
#include <cstdint>
#include <cstdio>
#include <vector>

// extern globals (as before)
extern std::vector<polyscope::SurfaceMesh*> g_poset_meshes_2d;
extern std::vector<polyscope::SurfaceMesh*> g_poset_meshes_3d;
extern bool g_show_poset_2d;
extern bool g_show_poset_3d;

extern pst::Poset1 g_poset1;
extern pst2::Poset2 g_poset2;

extern bool g_interval_active;
extern std::vector<char> g_interval_mask;
extern bool g_focus_active;
extern std::vector<char> g_focus_mask;

namespace viz_poset {
void rebuild_downflip_network_filtered();
void register_poset2_cover_edges(const pst2::Poset2& P2);

static void reset_visibility_to_full() {
    const int n = (int)g_poset1.nodes.size();
    g_interval_active = false; g_interval_mask.clear();
    g_focus_active = false;    g_focus_mask.clear();

    for (int i = 0; i < n; ++i) {
        if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setEnabled(g_show_poset_2d);
        if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setEnabled(g_show_poset_3d);
    }
    rebuild_downflip_network_filtered();
    register_poset2_cover_edges(g_poset2);
}

static void apply_mu_color(int idx, std::int64_t mu, int& cnt0, int& cntP, int& cntN, int& cntOther) {
    glm::vec3 c;
    if (mu == 0) { c = glm::vec3(0.f, 0.f, 0.f); cnt0++; }
    else if (mu == 1) { c = glm::vec3(1.f, 0.f, 0.f); cntP++; }
    else if (mu == -1) { c = glm::vec3(0.f, 0.f, 1.f); cntN++; }
    else { c = (mu > 0) ? glm::vec3(1.f, 0.f, 0.f) : glm::vec3(0.f, 0.f, 1.f); cntOther++; }

    if (g_poset_meshes_2d[idx]) g_poset_meshes_2d[idx]->setSurfaceColor(c);
    if (g_poset_meshes_3d[idx]) g_poset_meshes_3d[idx]->setSurfaceColor(c);
}

static int unique_min_from_cover_down(const std::vector<std::vector<int>>& cover_down) {
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

static int unique_max_from_cover_up(const std::vector<std::vector<int>>& cover_up) {
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

// compute μ(x, v) for all v using one DP (in up-closure of x)
static std::vector<std::int64_t> mobius_row_from_view(const mob::CoverUpView& V, int x) {
    const int N = (int)V.R.size();
    std::vector<std::int64_t> mu(N, 0);
    if (x < 0 || x >= N) return mu;

    std::vector<char> inUp(N, 0);
    for (int v = 0; v < N; ++v) {
        if (mob::leq(V, x, v)) inUp[v] = 1;
    }

    std::vector<int> topoU;
    topoU.reserve(N);
    for (int v : V.topo) if (inUp[v]) topoU.push_back(v);

    mu[x] = 1;

    for (int v : topoU) {
        if (v == x) continue;
        std::int64_t sum = 0;
        for (int z : topoU) {
            if (z == v) break;
            if (!mob::leq(V, z, v)) continue;
            sum += mu[z];
        }
        mu[v] = -sum;
    }
    return mu;
}

// ---------- PUBLIC COLORING API ----------
// P2: color μ(min, v)
void color_nodes_by_mobius_anchor_min_to_all_poset2() {
    const int n = (int)g_poset1.nodes.size();
    if (n == 0) return;

    reset_visibility_to_full();

    mob::CoverUpView V = mob::build_view(g_poset2.cover_up);
    int min0 = unique_min_from_cover_down(g_poset2.cover_down);
    if (min0 == -1) {
        std::cout << "[moebius] P2: no unique minimum, cannot color μ(min,v)\n";
        return;
    }

    auto mu = mobius_row_from_view(V, min0);

    int cnt0=0, cntP=0, cntN=0, cntOther=0;
    for (int v = 0; v < n; ++v) {
        apply_mu_color(v, mu[v], cnt0, cntP, cntN, cntOther);
    }

    std::cout << "[moebius] P2 colored nodes by μ(" << min0 << ",v): "
              << "0=" << cnt0 << " +1=" << cntP << " -1=" << cntN;
    if (cntOther) std::cout << " (|μ|!=1: " << cntOther << ")";
    std::cout << "\n";
}

// P1: color μ(min, v)
void color_nodes_by_mobius_anchor_min_to_all_poset1() {
    const int n = (int)g_poset1.nodes.size();
    if (n == 0) return;

    reset_visibility_to_full();

    mob::CoverUpView V = mob::build_view(g_poset1.cover_up);
    int min0 = unique_min_from_cover_down(g_poset1.cover_down);
    if (min0 == -1) {
        std::cout << "[moebius] P1: no unique minimum, cannot color μ(min,v)\n";
        return;
    }

    auto mu = mobius_row_from_view(V, min0);

    int cnt0=0, cntP=0, cntN=0, cntOther=0;
    for (int v = 0; v < n; ++v) {
        apply_mu_color(v, mu[v], cnt0, cntP, cntN, cntOther);
    }

    std::cout << "[moebius] P1 colored nodes by μ(" << min0 << ",v): "
              << "0=" << cnt0 << " +1=" << cntP << " -1=" << cntN;
    if (cntOther) std::cout << " (|μ|!=1: " << cntOther << ")";
    std::cout << "\n";
}

// μ(v,max): easiest via dual view = build_view(cover_down) (treating down edges as up edges in the dual)
static std::vector<std::int64_t> mobius_to_max_via_dual(
    const std::vector<std::vector<int>>& cover_down_as_up,
    int max_in_original)
{
    // In the dual, "max(original)" becomes "min(dual)".
    // μ_original(v, max) = μ_dual(max, v).
    mob::CoverUpView Vdual = mob::build_view(cover_down_as_up);
    return mobius_row_from_view(Vdual, max_in_original);
}

void color_nodes_by_mobius_all_to_anchor_max_poset2() {
    const int n = (int)g_poset1.nodes.size();
    if (n == 0) return;

    reset_visibility_to_full();

    int max1 = unique_max_from_cover_up(g_poset2.cover_up);
    if (max1 == -1) {
        std::cout << "[moebius] P2: no unique maximum, cannot color μ(v,max)\n";
        return;
    }

    auto mu_dual_row = mobius_to_max_via_dual(g_poset2.cover_down, max1);
    // mu_dual_row[v] = μ_dual(max, v) = μ_original(v, max)

    int cnt0=0, cntP=0, cntN=0, cntOther=0;
    for (int v = 0; v < n; ++v) {
        apply_mu_color(v, mu_dual_row[v], cnt0, cntP, cntN, cntOther);
    }

    std::cout << "[moebius] P2 colored nodes by μ(v," << max1 << "): "
              << "0=" << cnt0 << " +1=" << cntP << " -1=" << cntN;
    if (cntOther) std::cout << " (|μ|!=1: " << cntOther << ")";
    std::cout << "\n";
}

void color_nodes_by_mobius_all_to_anchor_max_poset1() {
    const int n = (int)g_poset1.nodes.size();
    if (n == 0) return;

    reset_visibility_to_full();

    int max1 = unique_max_from_cover_up(g_poset1.cover_up);
    if (max1 == -1) {
        std::cout << "[moebius] P1: no unique maximum, cannot color μ(v,max)\n";
        return;
    }

    auto mu_dual_row = mobius_to_max_via_dual(g_poset1.cover_down, max1);

    int cnt0=0, cntP=0, cntN=0, cntOther=0;
    for (int v = 0; v < n; ++v) {
        apply_mu_color(v, mu_dual_row[v], cnt0, cntP, cntN, cntOther);
    }

    std::cout << "[moebius] P1 colored nodes by μ(v," << max1 << "): "
              << "0=" << cnt0 << " +1=" << cntP << " -1=" << cntN;
    if (cntOther) std::cout << " (|μ|!=1: " << cntOther << ")";
    std::cout << "\n";
}


void color_node(int idx, glm::vec3 c2, glm::vec3 c3) {
    if (idx < 0 || idx >= (int)g_poset_meshes_2d.size()) return;
    if (g_poset_meshes_2d[idx]) g_poset_meshes_2d[idx]->setSurfaceColor(c2);
    if (g_poset_meshes_3d[idx]) g_poset_meshes_3d[idx]->setSurfaceColor(c3);
}


void reset_node_coloring() {
    for (int i = 0; i < (int)g_poset_meshes_2d.size(); ++i) {
        if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setSurfaceColor(glm::vec3(0.6f, 0.8f, 1.0f));
        if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setSurfaceColor(glm::vec3(0.2f, 0.4f, 0.8f));
    }
}


} // namespace viz_poset


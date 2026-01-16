#include "moebius.h"

#include <boost/dynamic_bitset.hpp>
#include <queue>
#include <iostream>
#include <algorithm>
#include <cstdlib>

#include <imgui.h>

namespace viz_poset {
    void color_nodes_by_mobius_0T();
    void reset_node_coloring();
    void color_nodes_by_mobius_poset1_0T();
}


namespace mob {


using Bitset = boost::dynamic_bitset<>;

static void sort_unique(std::vector<std::vector<int>>& out) {
    for (auto& nbrs : out) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
}

static std::vector<int> topo_sort_kahn(const std::vector<std::vector<int>>& out) {
    const int N = (int)out.size();
    std::vector<int> indeg(N, 0);

    for (int u = 0; u < N; ++u) {
        for (int v : out[u]) {
            if (v < 0 || v >= N) continue;
            indeg[v]++;
        }
    }

    std::queue<int> q;
    for (int i = 0; i < N; ++i) if (indeg[i] == 0) q.push(i);

    std::vector<int> topo;
    topo.reserve(N);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        topo.push_back(u);
        for (int v : out[u]) {
            if (v < 0 || v >= N) continue;
            if (--indeg[v] == 0) q.push(v);
        }
    }

    // If not all nodes appear, graph had a cycle or invalid edges; still return what we have.
    return topo;
}

static std::vector<Bitset> compute_reachability(
    const std::vector<std::vector<int>>& out,
    const std::vector<int>& topo)
{
    const int N = (int)out.size();
    std::vector<Bitset> R;
    R.reserve(N);
    for (int i = 0; i < N; ++i) R.emplace_back(N);

    for (int ti = (int)topo.size() - 1; ti >= 0; --ti) {
        int u = topo[ti];
        Bitset ru(N);
        for (int v : out[u]) {
            if (v < 0 || v >= N) continue;
            ru |= R[v];
            ru.set(v);
        }
        R[u] = std::move(ru);
    }
    return R;
}

// Build interval [x,y] as the set of z with x <= z <= y (using reachability)
// Assumes cover_out is upward
static std::vector<int> interval_xy_from_cover(
    const std::vector<std::vector<int>>& cover_out,
    int x, int y)
{
    const int N = (int)cover_out.size();
    if (x < 0 || x >= N || y < 0 || y >= N) return {};

    // reachability on the whole poset
    auto topo = topo_sort_kahn(cover_out);
    if ((int)topo.size() != N) {
        // if topo sort failed, still try but results may be meaningless
        std::cout << "[moebius] warning: topo size != N (cycle or invalid edges?)\n";
    }
    auto R = compute_reachability(cover_out, topo);

    // x <= y iff y reachable from x OR x==y
    const bool x_leq_y = (x == y) || R[x].test(y);
    if (!x_leq_y) return {};

    std::vector<int> interval;
    interval.reserve(N);

    for (int z = 0; z < N; ++z) {
        bool x_leq_z = (x == z) || R[x].test(z);
        bool z_leq_y = (z == y) || R[z].test(y);
        if (x_leq_z && z_leq_y) interval.push_back(z);
    }
    return interval;
}

// computation of möbius function from cover graph
std::int64_t mobius_xy_from_cover(const std::vector<std::vector<int>>& cover_out,int x, int y){
    // N is number of nodes in the poset
    // check if empty poset or node inputs are invalid
    const int N = (int)cover_out.size();
    if (N == 0) return 0;
    if (x < 0 || x >= N || y < 0 || y >= N) return 0;

    // build interval nodes 
    // get all nodes in the interval [x,y]
    std::vector<int> interval = interval_xy_from_cover(cover_out, x, y);
    if (interval.empty()) return 0; // incomparable => mu=0 

    const int M = (int)interval.size();

    // global -> local
    std::vector<int> id(N, -1);
    for (int i = 0; i < M; ++i) id[interval[i]] = i;

    const int lx = id[x]; // local index of x in interval
    const int ly = id[y]; // local index of y in interval
    if (lx < 0 || ly < 0) return 0;


    // build the subgraph [x,y] of the poset using local indices
    std::vector<std::vector<int>> adj(M);
    std::vector<int> indeg(M, 0);

    // collect edges in the interval
    for (int u_glob : interval) {
        int u = id[u_glob];
        for (int v_glob : cover_out[u_glob]) {
            if (v_glob < 0 || v_glob >= N) continue;
            int v = id[v_glob];
            if (v != -1) {
                adj[u].push_back(v);
                indeg[v]++;
            }
        }
    }

    // topo sort induced interval graph
    std::queue<int> q;
    for (int i = 0; i < M; ++i) if (indeg[i] == 0) q.push(i);

    // check if interval graph is a DAG
    std::vector<int> topo;
    topo.reserve(M);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        topo.push_back(u);
        for (int v : adj[u]) {
            if (--indeg[v] == 0) q.push(v);
        }
    }
    if ((int)topo.size() != M) {
        std::cout << "[moebius] interval induced graph not a DAG (unexpected)\n";
        return 0;
    }

    // reachability in interval graph using a bitset
    std::vector<Bitset> R(M, Bitset(M));
    for (int ti = M - 1; ti >= 0; --ti) {
        int u = topo[ti];
        R[u].set(u);
        for (int v : adj[u]) R[u] |= R[v];
    }

    // compute möbius function using the formula:
    // mu(x,v) = 0 if 
    // mu(x,x)= 1
    // mu(x,v) = - sum_{x<=z<v} mu(x,z)

    // stores function values mu(x,v) for v in interval
    std::vector<std::int64_t> mu(M, 0);
    mu[lx] = 1; // mu(x,x)= 1

    for (int v : topo) {
        if (v == lx) continue;

        std::int64_t sum = 0;
        for (int z = 0; z < M; ++z) {
            if (z == v) continue;
            if (R[z].test(v)) { // z <= v
                sum += mu[z];
            }
        }
        mu[v] = -sum;
    }

    return mu[ly];  
}

std::vector<std::vector<int>> build_cover_up_from_poset1_nodes(const std::vector<pst::Node>& nodes){
    const int N = (int)nodes.size();
    std::vector<std::vector<int>> cover_up(N);

    // nodes[u].children are down edges u(top)->v(bottom)
    // upward cover edge is v -> u
    for (int u = 0; u < N; ++u) {
        for (int v : nodes[u].children) {
            if (v < 0 || v >= N) continue;
            cover_up[v].push_back(u);
        }
    }

    sort_unique(cover_up);
    return cover_up;
}

// möbius comparison UI
void draw_mobius_compare_ui(const std::vector<pst::Node>& poset1_nodes,const pst2::Poset2& poset2) {
    static char x_buf[32] = "0";
    static char y_buf[32] = "0";

    static bool computed = false;
    static std::int64_t mu1 = 0;
    static std::int64_t mu2 = 0;

    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (ImGui::CollapsingHeader("möb(x,y) compare: poset1 vs poset2", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::InputText("x (lower)", x_buf, IM_ARRAYSIZE(x_buf));
        ImGui::InputText("y (upper)", y_buf, IM_ARRAYSIZE(y_buf));

        const int x = std::atoi(x_buf);
        const int y = std::atoi(y_buf);

        if (ImGui::Button("compute möb(x,y) in p1 & p2")) {
            auto cover1_up = build_cover_up_from_poset1_nodes(poset1_nodes);

            mu1 = mobius_xy_from_cover(cover1_up, x, y);
            mu2 = mobius_xy_from_cover(poset2.cover_out, x, y);
            computed = true;

            std::cout << "[moebius] poset1 moeb(" << x << "," << y << ") = " << mu1 << "\n";
            std::cout << "[moebius] poset2 moeb(" << x << "," << y << ") = " << mu2 << "\n";
        }

        ImGui::SameLine();
        if (ImGui::Button("print interval nodes")) {
            auto cover1_up = build_cover_up_from_poset1_nodes(poset1_nodes);

            auto I1 = interval_xy_from_cover(cover1_up, x, y);
            auto I2 = interval_xy_from_cover(poset2.cover_out, x, y);

            std::cout << "[moebius] interval poset1 [" << x << "," << y << "], size=" << I1.size() << "\n";
            // print z in one row
            for (int z : I1) std::cout << z << " ";
            std::cout << "\n";
           


            std::cout << "[moebius] interval poset2 [" << x << "," << y << "], size=" << I2.size() << "\n";
            // print z in one row
            for (int z : I2) std::cout << z << " ";
            std::cout << "\n";
        
        }

        ImGui::Separator();

        if (computed) {
            ImGui::Text("poset1: moeb(%d,%d) = %lld", x, y, (long long)mu1);
            ImGui::Text("poset2: moeb(%d,%d) = %lld", x, y, (long long)mu2);
            //ImGui::Text("difference: möb2 - möb1 = %lld", (long long)(mu2 - mu1));
        } else {
            //ImGui::TextUnformatted("Enter x,y then compute. Compare mu-values between poset1 and poset2.");
        }


        if (ImGui::Button("compute moeb(T,0) for all T in poset2")) {
            viz_poset::color_nodes_by_mobius_0T();
        }


        if (ImGui::Button("compute moeb(T,0) for all T in poset1")) {
            viz_poset::color_nodes_by_mobius_poset1_0T();
        }


        if (ImGui::Button("reset node coloring")) {
            viz_poset::reset_node_coloring();
        }

        ImGui::Text("black = 0, red = +1, blue = -1");


    }

}

} // namespace mob

#pragma once
#include <vector>
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <queue>
#include <stdexcept>
#include <limits>
#include <iostream>

namespace pst {

    
    using Bitset = boost::dynamic_bitset<>; 

    // sorts adjacency lists and removes duplicates
    // isntended for cover graphs (poset1 / poset2)
    inline void sort_unique_adjacency(std::vector<std::vector<int>>& adj) {
        for (auto& nbrs : adj) {
            std::sort(nbrs.begin(), nbrs.end());
            nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        }
    }
    // kahn topo sort, throws if a cycle is detected
    inline std::vector<int> topo_sort_kahn(const std::vector<std::vector<int>>& out) {
        const int N = static_cast<int>(out.size());
        std::vector<int> indeg(N, 0);

        for (int u = 0; u < N; ++u) {
            for (int v : out[u]) {
                if (v < 0 || v >= N) {
                    throw std::runtime_error("topo_sort_kahn: edge endpoint out of range");
                }
                indeg[v]++;
            }
        }

        std::queue<int> q;
        for (int i = 0; i < N; ++i) {
            if (indeg[i] == 0) q.push(i);
        }

        std::vector<int> topo;
        topo.reserve(N);

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            topo.push_back(u);

            for (int v : out[u]) {
                if (--indeg[v] == 0) q.push(v);
            }
        }

        if ((int)topo.size() != N) {
            throw std::runtime_error("topo_sort_kahn: graph has a directed cycle");
        }
        return topo;
    }

    // reachability bitsets R[u][v] = 1 iff v reachable from u
    inline std::vector<Bitset> compute_reachability(
        const std::vector<std::vector<int>>& out,
        const std::vector<int>& topo)
    {
        const int N = static_cast<int>(out.size());
        std::vector<Bitset> R;
        R.reserve(N);
        for (int i = 0; i < N; ++i) R.emplace_back(N);

        for (int ti = (int)topo.size() - 1; ti >= 0; --ti) {
            int u = topo[ti];
            Bitset ru(N);
            ru.set(u); // reflexive reachability (<=)

            for (int v : out[u]) {
                if (v < 0 || v >= N) continue;
                ru |= R[v];
                ru.set(v);
            }
            R[u] = std::move(ru);
        }
        return R;
    }

    // find a path from start to goal by going up along cover edges
    // cover_up[u] lists immediate successors above u (u < v)
    // used by meet/ join visualization
    inline std::vector<int> find_path_up_cover_up(
        const std::vector<std::vector<int>>& cover_up,
        int start,
        int goal,
        const std::vector<char>* allowed_mask = nullptr) {
        const int n = (int)cover_up.size();
        if (start < 0 || start >= n || goal < 0 || goal >= n) return {};
        if (start == goal) return { start };

        auto ok = [&](int v) {
            return !allowed_mask || (v >= 0 && v < (int)allowed_mask->size() && (*allowed_mask)[v]);
        };
        if (!ok(start) || !ok(goal)) return {};

        std::vector<int> parent(n, -1);
        std::queue<int> q;

        parent[start] = start;
        q.push(start);

        bool found = false;

        while (!q.empty() && !found) {
            int u = q.front(); q.pop();

            for (int v : cover_up[u]) { // go UP: u -> v
                if (v < 0 || v >= n) continue;
                if (!ok(v)) continue;
                if (parent[v] != -1) continue;

                parent[v] = u;

                if (v == goal) {
                    found = true;
                    break;
                }
                q.push(v);
            }
        }

        if (parent[goal] == -1) return {}; // no path found

        // reconstruct
        std::vector<int> path;
        for (int cur = goal; cur != start; cur = parent[cur]) {
            path.push_back(cur);
        }
        path.push_back(start);
        std::reverse(path.begin(), path.end());
        return path;
    }



inline std::vector<int>compute_levels_longest_from_root_cover_down(
    const std::vector<std::vector<int>>& cover_down,
    int root = 0){
    const int n = (int)cover_down.size();
    std::vector<int> level(n, -1);
    if (n == 0) return level;
    if (root < 0 || root >= n) return level;

    std::vector<int> topo;
    try {
        topo = pst::topo_sort_kahn(cover_down);
    } catch (const std::exception& e) {
        // fallback: BFS levels from root (works even with cycles)
        std::vector<int> dist(n, -1);
        std::queue<int> q;
        dist[root] = 0;
        q.push(root);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : cover_down[u]) {
                if (v < 0 || v >= n) continue;
                if (dist[v] != -1) continue;
                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }

        int max_lv = 0;
        for (int i = 0; i < n; ++i) {
            if (dist[i] != -1) {
                level[i] = dist[i];
                max_lv = std::max(max_lv, level[i]);
            }
        }

        // place all unreachable after reachable block, stable order
        int next = max_lv + 1;
        for (int i = 0; i < n; ++i) {
            if (level[i] == -1) level[i] = next++;
        }

        std::cerr << "[vis_poset] WARNING: cover_down has a directed cycle; "
                     "using BFS fallback for levels. (" << e.what() << ")\n";
        return level;
    }

    // Longest path from root in DAG (standard DP over topo)
    const int NEG_INF = std::numeric_limits<int>::min();
    std::vector<int> dist(n, NEG_INF);
    dist[root] = 0;

    for (int u : topo) {
        if (dist[u] == NEG_INF) continue;
        for (int v : cover_down[u]) {
            if (v < 0 || v >= n) continue;
            dist[v] = std::max(dist[v], dist[u] + 1);
        }
    }

    int max_reachable = 0;
    for (int i = 0; i < n; ++i) {
        if (dist[i] != NEG_INF) {
            level[i] = dist[i];
            max_reachable = std::max(max_reachable, level[i]);
        }
    }

    // Unreachable nodes: keep non-negative and deterministic.
    // We’ll compute a “component depth” within the subgraph induced by unreachable nodes.
    // Simple approach: topo DP but only on unreachable nodes, starting from their sources.
    std::vector<int> indeg_un(n, 0);
    for (int u = 0; u < n; ++u) {
        if (level[u] != -1) continue; // only unreachable
        for (int v : cover_down[u]) {
            if (v < 0 || v >= n) continue;
            if (level[v] != -1) continue;
            indeg_un[v]++;
        }
    }

    std::queue<int> q;
    for (int i = 0; i < n; ++i) {
        if (level[i] == -1 && indeg_un[i] == 0) q.push(i);
    }

    // if there are unreachable cycles (shouldn’t if whole graph is DAG, but safe):
    // we’ll still assign them sequentially at the end.
    std::vector<int> dist_un(n, -1);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        if (dist_un[u] == -1) dist_un[u] = 0;
        for (int v : cover_down[u]) {
            if (v < 0 || v >= n) continue;
            if (level[v] != -1) continue; // only unreachable subgraph
            dist_un[v] = std::max(dist_un[v], dist_un[u] + 1);
            if (--indeg_un[v] == 0) q.push(v);
        }
    }

    int base = max_reachable + 1;
    int next = base;

    for (int i = 0; i < n; ++i) {
        if (level[i] != -1) continue;
        if (dist_un[i] != -1) {
            level[i] = base + dist_un[i];
        } else {
            // unreachable cycle / weirdness: place after everything
            level[i] = next++;
        }
    }

    return level;
}



// transitive reduction for a DAG using reachability bitsets
// for each edge u -> v, look at all neighbors w of u and check if w -> v 
static std::vector<std::vector<int>> transitive_reduction(
    const std::vector<std::vector<int>>& out,
    const std::vector<Bitset>& R)
{
    const int N = static_cast<int>(out.size());
    std::vector<std::vector<int>> cover_out(N);

    for (int u = 0; u < N; ++u) {
        // for each edge u->v, check if there exists w!=v in out[u] with w -> v
        for (int v : out[u]) {
            bool redundant = false;

            for (int w : out[u]) {
                if (w == v) continue;
                if (R[w].test(v)) {
                    redundant = true;
                    break;
                }
            }

            if (!redundant) {
                cover_out[u].push_back(v);
            }
        }
    }

    pst::sort_unique_adjacency(cover_out);
    return cover_out;
}





} 

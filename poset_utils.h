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
    inline std::vector<Bitset> compute_reachability(const std::vector<std::vector<int>>& out,const std::vector<int>& topo){
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

    
    // transitive reduction for a DAG using reachability bitsets
    // for each edge u -> v, look at all neighbors w of u and check if w -> v 
    inline std::vector<std::vector<int>> transitive_reduction(const std::vector<std::vector<int>>& out,const std::vector<Bitset>& R){
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

} 

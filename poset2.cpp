
#include "poset2.h"
#include "compare.h"
#include "input.h"

#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <queue>
#include <stdexcept>
#include <algorithm>
#include <functional>


namespace pst2 {

using Bitset = boost::dynamic_bitset<>;


// helper function to sort adjacency lists
// sorts neighbors and removes duplicates
// same edge could be added more than once 
static void sort_unique_adjacency(std::vector<std::vector<int>>& out) {
    for (auto& nbrs : out) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
}

// kahn topo sort, throws if a cycle is detected
static std::vector<int> topo_sort_kahn(const std::vector<std::vector<int>>& out) {
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
            indeg[v]--;
            if (indeg[v] == 0) q.push(v);
        }
    }

    if ((int)topo.size() != N) {
        throw std::runtime_error("topo_sort_kahn: graph has a directed cycle");
    }
    return topo;
}

// reachability bitsets R[u][v] = 1 iff v reachable from u 
static std::vector<Bitset> compute_reachability(
    const std::vector<std::vector<int>>& out,
    const std::vector<int>& topo)
{
    const int N = static_cast<int>(out.size());
    std::vector<Bitset> R;
    R.reserve(N);
    for (int i = 0; i < N; ++i) R.emplace_back(N);

    // reverse topo
    for (int ti = N - 1; ti >= 0; --ti) {
        int u = topo[ti];
        Bitset ru(N);

        for (int v : out[u]) {
            ru |= R[v];
            ru.set(v);
        }
        R[u] = std::move(ru);
    }
    return R;
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

    sort_unique_adjacency(cover_out);
    return cover_out;
}

// build out1 adjacency from pst::Node::children
static std::vector<std::vector<int>> build_out1_from_nodes(const std::vector<pst::Node>& nodes)
{
    const int N = static_cast<int>(nodes.size());
    std::vector<std::vector<int>> out1(N);
    for (int u = 0; u < N; ++u) {
        out1[u] = nodes[u].children; // <=_1 down edges
    }
    sort_unique_adjacency(out1);
    return out1;
}


int print_cover_relations(const pst2::Poset2& P2) {
    int counter = 0;
    std::cout << "\n=== <=2 cover relations top -> down ===\n";
    for (int u = 0; u < (int)P2.cover_out.size(); ++u) {
        for (int v : P2.cover_out[u]) {
            // stored as u <=2 v, but print as v -> u (top -> down)
            std::cout << v << " -> " << u << "\n";
            counter++;
        }
    }
    std::cout << "=== end of <=2 cover relations ===\n\n";
    return counter;
}


// fuunction to build the <=2 poset from triangulation comparison
Poset2 build_poset2(const df::InputData& D, const std::vector<pst::Node>& nodes, bool keep_full_out2){
    const int N = static_cast<int>(nodes.size());
    if (N == 0) return Poset2{};

    // 1) <=1 adjacency and reachability 
    // out1[u] = list of v with u -> v by down edges
    // to check which pairs are already ordered by <=1
    std::vector<std::vector<int>> out1 = build_out1_from_nodes(nodes);

    // topo sort <=1 graph so that all edges go forward 
    // confirm that graph is a DAG
    // topo1 vector of node indices in topo order
    std::vector<int> topo1 = topo_sort_kahn(out1);

    // R1[u] is a bitset of N bits with R1[u][v] = 1 iff v reachable from u by <=1 down edges
    std::vector<Bitset> R1 = compute_reachability(out1, topo1);

    // 2) build <=2 edge set out2
    std::vector<std::vector<int>> out2(N);

    // add all <=_1 implied relations for free 
    // adding all transitive pairs as edges would be huge
    // so we only add the existing <=_1 edges here (out1)
    // and rely on the comparator for pairs not related by <=_1 reachability
        
    // out2 starts empty; fill it from comparisons + optionally out1 edges
    // initialize out2 with <=1 edges (REVERSED: child <=2 parent)
    for (int u = 0; u < N; ++u) {
        for (int v : out1[u]) {
            out2[v].push_back(u);
        }
    }


    // (b) resolve <=_1-incomparable pairs with compare function
    // loop unordered pairs u<v
    for (int u = 0; u < N; ++u) {
        for (int v = u + 1; v < N; ++v) {

            // skip if <=1 already orders them, meaning R1[u][v] = 1 or R1[v][u] = 1
            if (R1[u].test(v) || R1[v].test(u)) {
                // R1[u].test(v) means: u -> v by down edges  =>  v <=1 u  =>  v <=2 u already implied
                continue;
            }

            bool result1 = pst2::compare(u, v, D, nodes);
            bool result2 = pst2::compare(v, u, D, nodes);

            // result1 and result2 are false: incomparable
            // if result1 is true and result2 is false: t1 <_2 t2
            // if result1 is false and result2 is true: t2 <_2 t1

            
            // result == false: 
            if (result1 == false && result2 == false) {
                // incomparable: no edge
            } else if (result1 == true && result2 == false) {
                out2[u].push_back(v);
            } else if (result1 == false && result2 == true) {
                out2[v].push_back(u);
            } else if (result1 == true && result2 == true) {
                // print warning: should not happen unless triangulations are equal 
                std::cerr << "Warning: build_poset2: triangulations " << u << " and " << v << " compare equal in both directions\n";
            }

        }
    }
    // deleting double edges and sorting (sorting not necessary but cleaner)
    sort_unique_adjacency(out2);

    // 3) topo sort <=2 graph (must be DAG)
    std::vector<int> topo2 = topo_sort_kahn(out2);

    // 4) reachability for <=_2 graph 
    // R2[u][v] = 1 iff v reachable from u by <=2 edges
    std::vector<Bitset> R2 = compute_reachability(out2, topo2);

    // 5) transitive reduction => covering edges 
    // cover_out[u] = list of v with u -> v by covering edges
    // cover_out edges go upward: u <=2 v stored as u -> v
    // stored at u as outgoing edge to v
    std::vector<std::vector<int>> cover_out = transitive_reduction(out2, R2);

    Poset2 P2;
    P2.cover_out = std::move(cover_out);

    
    if (keep_full_out2) {
        P2.out2 = std::move(out2);
    }
    
    return P2;
}

// this section is for computing intervals in the poset
// build cover_in adjacency from cover_out
// cover_in[v] = list of u with u -> v by covering edges
// cover_in edges go downward: u <=2 v stored as u -> v, 
// stored at v as incoming edge from u
static std::vector<std::vector<int>> build_cover_in(const Poset2& P) {
    const int N = (int)P.cover_out.size();
    std::vector<std::vector<int>> cover_in(N);
    for (int u = 0; u < N; ++u) {
        for (int v : P.cover_out[u]) {
            if (v < 0 || v >= N) continue;
            cover_in[v].push_back(u);
        }
    }
    for (auto& nbrs : cover_in) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
    return cover_in;
}

// compute upper set of x by going up along <=2 cover edges
// uses cover_out (up edges), start at x and do BFS upwards 
static std::vector<char> upper_set_x1(const Poset2& P, int x){
    const int N = (int)P.cover_out.size();
    // visited array
    std::vector<char> vis(N, 0);
    // BFS queue
    std::queue<int> q;

    // start at x
    vis[x] = 1;
    // x visited 
    q.push(x);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : P.cover_out[u]) { // u <=2 v (go up)
            if (!vis[v]) { vis[v] = 1; q.push(v); }
        }
    }
    // return list of visited nodes 
    return vis;
}

// compute lower set of y by going down along <=2 cover edges
// uses cover_in (down edges), start at y and do BFS downwards
static std::vector<char> lower_set_0y(const std::vector<std::vector<int>>& cover_in, int y){
    const int N = (int)cover_in.size();
    std::vector<char> vis(N, 0);
    std::queue<int> q;

    vis[y] = 1;
    q.push(y);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : cover_in[u]) { // v <=2 u (go down)
            if (!vis[v]) { vis[v] = 1; q.push(v); }
        }
    }
    return vis;
}


// build interval [x,y] = { z | x <=2 z <=2 y }
// this is a subgraph of the poset2 defined by cover_out
std::vector<int> interval_xy(const Poset2& P, int x, int y) {
    int N = (int)P.cover_out.size();
    if (x < 0 || x >= N || y < 0 || y >= N) {
        std::cout << "invalid node indices\n";
        return {};
    }

    // compute reachability
    auto up_x = upper_set_x1(P, x);
    auto up_y = upper_set_x1(P, y);

    
    bool x_le_y = up_x[y];
    bool y_le_x = up_y[x];

   // check if x <=2 y or y <=2 x
    std::cout << "relation: ";
    if (x == y) std::cout << "x == y\n";        
    else if (x_le_y) std::cout << "x <=2 y\n";
    else if (y_le_x) std::cout << "y <=2 x\n";
    else std::cout << "incomparable\n";
    

    // interval [x,y] only exists if x <=2 y (or x==y)
    if (!(x == y || x_le_y)) {
        std::cout << "interval [" << x << "," << y << "] is empty\n";
        return {};
    }

    // [x,1] intersection [0,y]
    auto cover_in = build_cover_in(P);
    auto low_y = lower_set_0y(cover_in, y);

    std::vector<int> res;
    for (int i = 0; i < N; ++i)
        if (up_x[i] && low_y[i])
            res.push_back(i);

    std::cout << "interval [" << x << "," << y << "] has "              
    << res.size() << " nodes\n";

    return res;
}



// compute meet candidates of x and y: maximal elements of [0,x] intersection [0,y]
std::vector<int> meet_candidates_xy(const Poset2& P, int x, int y) {
    int N = (int)P.cover_out.size();
    if (x < 0 || x >= N || y < 0 || y >= N) {
        std::cout << "meet_candidates_xy: invalid indices\n";
        return {};
    }

    // L = [0,x] intersection [0,y]
    auto cover_in = build_cover_in(P);
    auto low_x = lower_set_0y(cover_in, x);
    std::cout << "[0," << x << "]: ";
    for (int i = 0; i < N; ++i) {
        if (low_x[i]) std::cout << i << " ";
    }
    std::cout << "\n";

    auto low_y = lower_set_0y(cover_in, y);
    std::cout << "[0," << y << "]: ";
    for (int i = 0; i < N; ++i) {
        if (low_y[i]) std::cout << i << " ";
    }
    std::cout << "\n";

    std::vector<char> inL(N, 0);
    for (int i = 0; i < N; ++i) {
        if (low_x[i] && low_y[i]) inL[i] = 1;
    }

    // print L
    std::cout << "L = [0," << x << "] intersection [0," << y << "]: ";
    for (int i = 0; i < N; ++i) {
        if (inL[i]) std::cout << i << " ";
    }
    std::cout << "\n";

    // maximal elements of L: no outgoing cover edge to another element in L
    std::vector<int> maxima;
    for (int z = 0; z < N; ++z) {
        if (!inL[z]) continue;

        bool has_bigger_in_L = false;
        for (int w : P.cover_out[z]) {      // z <2 w
            if (w >= 0 && w < N && inL[w]) { // bigger and still in L
                has_bigger_in_L = true;
                break;
            }
        }

        if (!has_bigger_in_L) maxima.push_back(z);
    }

    return maxima;
}

// compute join candidates of x and y: minimal elements of [x,1] intersection [y,1]
std::vector<int> join_candidates_xy(const Poset2& P, int x, int y) {
    // check indices of x and y
    int N = (int)P.cover_out.size();
    if (x < 0 || x >= N || y < 0 || y >= N) {
        std::cout << "join_candidates_xy: invalid indices\n";
        return {};
    }

    // U = [x,1] intersection [y,1]
    auto up_x = upper_set_x1(P, x);
    std::cout << "[" << x << ",1]: ";
    for (int i = 0; i < N; ++i) {
        if (up_x[i]) std::cout << i << " ";
    }
    std::cout << "\n";

    auto up_y = upper_set_x1(P, y);
    std::cout << "[" << y << ",1]: ";
    for (int i = 0; i < N; ++i) {
        if (up_y[i]) std::cout << i << " ";
    }
    std::cout << "\n";

    // intersection
    std::vector<char> inU(N, 0);
    for (int i = 0; i < N; ++i) {
        if (up_x[i] && up_y[i]) inU[i] = 1;
    }
    // print U
    std::cout << "U = [" << x << ",1] intersection [" << y << ",1]: ";
    for (int i = 0; i < N; ++i) {
        if (inU[i]) std::cout << i << " ";
    }
    std::cout << "\n";

    // minima of U: no incoming cover edge from another element in U
    auto cover_in = build_cover_in(P);

    std::vector<int> minima;
    for (int z = 0; z < N; ++z) {
        if (!inU[z]) continue;

        bool has_smaller_in_U = false;
        for (int w : cover_in[z]) {          // w <2 z
            if (w >= 0 && w < N && inU[w]) { // smaller and still in U
                has_smaller_in_U = true;
                break;
            }
        }

        if (!has_smaller_in_U) minima.push_back(z);
    }

    return minima;
}



} // namespace pst2

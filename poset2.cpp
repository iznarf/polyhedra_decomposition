
#include "poset2.h"
#include "compare.h"
#include "input.h"
#include "poset.h"
#include "poset_utils.h"

#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <queue>
#include <stdexcept>
#include <algorithm>
#include <functional>


namespace pst2 {

using pst::Bitset;
 

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


// function to build the <=2 poset from triangulation comparison
//
// GOAL / CONVENTION (same as Poset1):
//   cover_down[u] = list of v such that u covers v  (v is immediately below u)
//   cover_up[u]   = list of v such that v covers u  (v is immediately above u)
//
// IMPORTANT:
//   We BUILD everything using DOWN edges (like Poset1.cover_down), because the rest of the code / intuition
//   and (likely) the comparator is easier to use in that direction.
//   In the end we derive cover_up as the reverse adjacency of cover_down.
Poset2 build_poset2(const df::InputData& D, const pst::Poset1& P1) {
    const int N = static_cast<int>(P1.nodes.size());
    if (N == 0) return Poset2{};

    const auto& nodes = P1.nodes;

    // 1) <=1 adjacency and reachability (DOWN direction)
    // out1_down[u] = list of v with u -> v by <=1 covering edges going DOWN
    // to check which pairs are already ordered by <=1
    std::vector<std::vector<int>> out1_down = P1.cover_down;
    pst::sort_unique_adjacency(out1_down);

    // topo sort <=1 graph so that all edges go forward
    // confirm that graph is a DAG
    // topo1 vector of node indices in topo order
    std::vector<int> topo1 = pst::topo_sort_kahn(out1_down);

    // R1[u] is a bitset of N bits with R1[u][v] = 1 iff v reachable from u by <=1 DOWN edges
    // i.e. u ->* v by down steps  <=>  v <=1 u
    std::vector<Bitset> R1 = pst::compute_reachability(out1_down, topo1);

    // 2) build <=2 edge set out2_down (DOWN direction)
    // out2_down[u] = list of v with u -> v meaning u is ABOVE v in <=2 (v <=2 u)
    std::vector<std::vector<int>> out2_down(N);

    // add all <=_1 implied relations for free
    // adding all transitive pairs as edges would be huge
    // so we only add the existing <=_1 cover edges here (out1_down == P1.cover_down)
    // and rely on the comparator for pairs not related by <=_1 reachability
    out2_down = out1_down;

    // (b) resolve <=_1-incomparable pairs with compare function
    // loop unordered pairs u<v
    for (int u = 0; u < N; ++u) {
        for (int v = u + 1; v < N; ++v) {

            // skip if <=1 already orders them, meaning R1[u][v] = 1 or R1[v][u] = 1
            // R1[u].test(v) means: u ->* v by down edges  =>  v <=1 u  =>  v <=2 u already implied
            if (R1[u].test(v) || R1[v].test(u)) {
                continue;
            }

            bool result1 = pst2::compare(u, v, D, nodes);
            bool result2 = pst2::compare(v, u, D, nodes);

            // result1 and result2 are false: incomparable
            // if result1 is true and result2 is false: u <=2 v
            // if result1 is false and result2 is true: v <=2 u
            //
            // NOTE (DOWN-edge storage):
            // if u <=2 v (u below v), then we store the DOWN edge v -> u.
            // if v <=2 u (v below u), then we store the DOWN edge u -> v.

            if (result1 == false && result2 == false) {
                // incomparable: no edge
            } else if (result1 == true && result2 == false) {
                // u <=2 v  =>  add down edge v -> u
                out2_down[v].push_back(u);
            } else if (result1 == false && result2 == true) {
                // v <=2 u  =>  add down edge u -> v
                out2_down[u].push_back(v);
            } else if (result1 == true && result2 == true) {
                // print warning: should not happen unless triangulations are equal
                std::cerr << "Warning: build_poset2: triangulations " << u << " and " << v
                          << " compare equal in both directions\n";
            }
        }
    }

    // deleting double edges and sorting (sorting not necessary but cleaner)
    pst::sort_unique_adjacency(out2_down);

    // 3) topo sort <=2 graph (must be DAG)
    // NOTE: this is still a DAG in DOWN direction (edges go from bigger to smaller)
    std::vector<int> topo2 = pst::topo_sort_kahn(out2_down);

    // 4) reachability for <=_2 graph
    // R2[u][v] = 1 iff v reachable from u by <=2 DOWN edges
    std::vector<Bitset> R2 = pst::compute_reachability(out2_down, topo2);

    // 5) transitive reduction => covering edges
    // cover_down[u] = list of v with u -> v by covering edges (DOWN direction)
    // i.e. u covers v in <=2 (v is immediately below u)
    std::vector<std::vector<int>> cover_down = transitive_reduction(out2_down, R2);

    // build cover_up as reverse adjacency of cover_down
    // cover_up[x] = list of y such that y covers x (y immediately above x)
    std::vector<std::vector<int>> cover_up(N);
    for (int u = 0; u < N; ++u) {
        for (int v : cover_down[u]) {
            cover_up[v].push_back(u);
        }
    }
    pst::sort_unique_adjacency(cover_up);

    Poset2 P2;
    P2.cover_down = std::move(cover_down);
    P2.cover_up   = std::move(cover_up);

    return P2;
}



void debug_compare_poset1(const df::InputData& D, const pst::Poset1& P1) {
    const auto& nodes = P1.nodes;
    const int n = static_cast<int>(nodes.size());

    std::cout << "\n=== comparator test on <=1 edges ===\n";

    int fail = 0;
    int total = 0;

    // ------------------------------------------------------------
    // 1) Check all <=1 (down) cover edges
    // ------------------------------------------------------------
    for (int u = 0; u < n; ++u) {
        for (int v : P1.cover_down[u]) {
            ++total;

            // v <=1 u   (down edge u -> v)
            if (!pst2::compare(v, u, D, nodes)) {
                ++fail;
                std::cout
                    << "FAIL: v<=1u but compare(v,u)=false: "
                    << v << " <=1 " << u
                    << "   (edge stored as " << u << " -> " << v << ")\n";
            }
        }
    }

    std::cout
        << "checked " << total
        << " <=1 edges, failures = " << fail << "\n";


    // ------------------------------------------------------------
    // 2) Full compare table
    // ------------------------------------------------------------
    std::cout << "\n=== full compare table (<=2 relation) ===\n";

    for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
            bool uv = pst2::compare(u, v, D, nodes);
            bool vu = pst2::compare(v, u, D, nodes);

            if (uv && vu) {
                std::cout << u << " == " << v << "\n";
            } else if (uv) {
                std::cout << u << " <2 " << v << "\n";
            } else if (vu) {
                std::cout << v << " <2 " << u << "\n";
            } else {
                std::cout << u << " || " << v << "\n";
            }
        }
    }
}

static int poset_size_checked(const Poset2& P) {
    int Nu = (int)P.cover_up.size();
    int Nd = (int)P.cover_down.size();
    if (Nu != Nd) {
        std::cout << "[poset2] WARNING: cover_up.size=" << Nu
                  << " != cover_down.size=" << Nd << "\n";
    }
    return std::max(Nu, Nd);
}


std::vector<int> minima(const Poset2& P) {
  int N = (int)P.cover_up.size();
  std::vector<int> mins;
  for (int v = 0; v < N; ++v)
    if (P.cover_down[v].empty()) mins.push_back(v);
  return mins;
}

std::vector<int> maxima(const Poset2& P) {
  int N = (int)P.cover_up.size();
  std::vector<int> maxs;
  for (int v = 0; v < N; ++v)
    if (P.cover_up[v].empty()) maxs.push_back(v);
  return maxs;
}

int unique_min_or_minus1(const Poset2& P) {
  auto mins = minima(P);
  return mins.size() == 1 ? mins[0] : -1;
}

int unique_max_or_minus1(const Poset2& P) {
  auto maxs = maxima(P);
  return maxs.size() == 1 ? maxs[0] : -1;
}

static std::vector<char> upper_set(const Poset2& P, int x) {
    const int N = poset_size_checked(P);
    std::vector<char> vis(N, 0);
    if (x < 0 || x >= N) return vis;

    std::queue<int> q;
    vis[x] = 1;
    q.push(x);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        if (u < 0 || u >= (int)P.cover_up.size()) continue;
        for (int v : P.cover_up[u]) { // go UP
            if (v < 0 || v >= N) continue;
            if (!vis[v]) { vis[v] = 1; q.push(v); }
        }
    }
    return vis;
}

static std::vector<char> lower_set(const Poset2& P, int y) {
    const int N = poset_size_checked(P);
    std::vector<char> vis(N, 0);
    if (y < 0 || y >= N) return vis;

    std::queue<int> q;
    vis[y] = 1;
    q.push(y);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        if (u < 0 || u >= (int)P.cover_down.size()) continue;
        for (int v : P.cover_down[u]) { // go DOWN
            if (v < 0 || v >= N) continue;
            if (!vis[v]) { vis[v] = 1; q.push(v); }
        }
    }
    return vis;
}


std::vector<int> interval_xy(const Poset2& P, int x, int y) {
  int N = (int)P.cover_up.size();
  if (x<0||x>=N||y<0||y>=N) return {};
  auto upx = upper_set(P, x);
  if (!upx[y] && x != y) return {};          // x not <= y
  auto lowy = lower_set(P, y);
  std::vector<int> res;
  for (int i=0;i<N;++i) if (upx[i] && lowy[i]) res.push_back(i);
  return res;
}

std::vector<int> interval_min_x(const Poset2& P, int x) {
  int N = (int)P.cover_up.size();
  if (x<0||x>=N) return {};
  auto lowx = lower_set(P, x);
  std::vector<int> res;
  for (int i=0;i<N;++i) if (lowx[i]) res.push_back(i);
  return res;
}

std::vector<int> interval_x_max(const Poset2& P, int x) {
  int N = (int)P.cover_up.size();
  if (x<0||x>=N) return {};
  auto upx = upper_set(P, x);
  std::vector<int> res;
  for (int i=0;i<N;++i) if (upx[i]) res.push_back(i);
  return res;
}



// meet candidates of x and y: maximal elements of lower_set(x) ∩ lower_set(y)
std::vector<int> meet_candidates_xy(const Poset2& P, int x, int y) {
    int N = (int)P.cover_up.size();
    if (x < 0 || x >= N || y < 0 || y >= N) {
        std::cout << "meet_candidates_xy: invalid indices\n";
        return {};
    }

    auto low_x = lower_set(P, x);  // down-closure of x
    auto low_y = lower_set(P, y);  // down-closure of y

    // intersection L
    std::vector<char> inL(N, 0);
    for (int i = 0; i < N; ++i) inL[i] = (low_x[i] && low_y[i]);

    // debug prints (optional)
    std::cout << "down(" << x << ") intersection down(" << y << "): ";
    for (int i = 0; i < N; ++i) if (inL[i]) std::cout << i << " ";
    std::cout << "\n";

    // maximal elements of L: no UP neighbor still in L
    std::vector<int> maxima;
    for (int z = 0; z < N; ++z) {
        if (!inL[z]) continue;

        bool has_bigger_in_L = false;
        for (int w : P.cover_up[z]) {         // z < w
            if (w >= 0 && w < N && inL[w]) {  // still in L
                has_bigger_in_L = true;
                break;
            }
        }
        if (!has_bigger_in_L) maxima.push_back(z);
    }

    return maxima;
}


// join candidates of x and y: minimal elements of upper_set(x) ∩ upper_set(y)
std::vector<int> join_candidates_xy(const Poset2& P, int x, int y) {
    int N = (int)P.cover_up.size();
    if (x < 0 || x >= N || y < 0 || y >= N) {
        std::cout << "join_candidates_xy: invalid indices\n";
        return {};
    }

    auto up_x = upper_set(P, x);  // up-closure of x
    auto up_y = upper_set(P, y);  // up-closure of y

    // intersection U
    std::vector<char> inU(N, 0);
    for (int i = 0; i < N; ++i) inU[i] = (up_x[i] && up_y[i]);

    // debug prints (optional)
    std::cout << "up(" << x << ")  up(" << y << "): ";
    for (int i = 0; i < N; ++i) if (inU[i]) std::cout << i << " ";
    std::cout << "\n";

    // minimal elements of U: no DOWN neighbor still in U
    std::vector<int> minima;
    for (int z = 0; z < N; ++z) {
        if (!inU[z]) continue;

        bool has_smaller_in_U = false;
        for (int w : P.cover_down[z]) {       // w < z
            if (w >= 0 && w < N && inU[w]) {
                has_smaller_in_U = true;
                break;
            }
        }
        if (!has_smaller_in_U) minima.push_back(z);
    }

    return minima;
}

}
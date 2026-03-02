#include "ear_star_tri_word.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace earstar {

// hash for sorted face triple
struct FaceKey {
    df::vertex_id a, b, c; // always sorted

    bool operator==(const FaceKey& o) const {
        return a == o.a && b == o.b && c == o.c;
    }
};

struct FaceKeyHash {
    std::size_t operator()(FaceKey const& k) const noexcept {
        std::size_t h = 0;
        auto mix = [&](std::size_t x) {
            h ^= x + 0x9e3779b97f4a7c15ull + (h<<6) + (h>>2);
        };
        mix((std::size_t)k.a);
        mix((std::size_t)k.b);
        mix((std::size_t)k.c);
        return h;
    }
};

static inline FaceKey make_face_key(df::vertex_id x, df::vertex_id y, df::vertex_id z) {
    df::vertex_id t[3] = {x,y,z};
    std::sort(t, t+3);
    return FaceKey{t[0], t[1], t[2]};
}

//precompute faces + star adjacency
struct Precomp {
    std::unordered_set<FaceKey, FaceKeyHash> faces;
    std::vector<char> star_adj; // size n_boundary, star_adj[x]=1 iff star adjacent to x
};

static Precomp precompute(const df::Tri2& tri, int n_boundary, df::vertex_id star_id) {
    Precomp P;
    P.faces.reserve((std::size_t)tri.number_of_faces() * 2);
    P.star_adj.assign(n_boundary, 0);

    for (auto f = tri.finite_faces_begin(); f != tri.finite_faces_end(); ++f) {
        df::vertex_id v0 = f->vertex(0)->info();
        df::vertex_id v1 = f->vertex(1)->info();
        df::vertex_id v2 = f->vertex(2)->info();

        P.faces.insert(make_face_key(v0, v1, v2));

        // star adjacency: if face contains star and boundary x, mark star_adj[x]
        if (v0 == star_id) {
            if ((int)v1 >= 0 && (int)v1 < n_boundary) P.star_adj[(int)v1] = 1;
            if ((int)v2 >= 0 && (int)v2 < n_boundary) P.star_adj[(int)v2] = 1;
        } else if (v1 == star_id) {
            if ((int)v0 >= 0 && (int)v0 < n_boundary) P.star_adj[(int)v0] = 1;
            if ((int)v2 >= 0 && (int)v2 < n_boundary) P.star_adj[(int)v2] = 1;
        } else if (v2 == star_id) {
            if ((int)v0 >= 0 && (int)v0 < n_boundary) P.star_adj[(int)v0] = 1;
            if ((int)v1 >= 0 && (int)v1 < n_boundary) P.star_adj[(int)v1] = 1;
        }
    }

    return P;
}

static inline bool is_face(const Precomp& P, df::vertex_id a, df::vertex_id b, df::vertex_id c) {
    return P.faces.find(make_face_key(a,b,c)) != P.faces.end();
}

//does triangulation contain the star vertex?
static bool has_star_vertex(const df::Tri2& tri, df::vertex_id star_id) {
    for (auto v = tri.finite_vertices_begin(); v != tri.finite_vertices_end(); ++v) {
        if (v->info() == star_id) return true;
    }
    return false;
}

//boundary cycle bookkeeping
struct BoundaryCycle {
    int n;
    std::vector<int> prev;
    std::vector<int> next;
    std::vector<char> alive;

    // initialize as 0--1--2--...--(n-1) in a cycle, all alive
    BoundaryCycle(int n_) : n(n_), prev(n_), next(n_), alive(n_, 1) {
        for (int i = 0; i < n; ++i) {
            // prev[i] = i-1, next[i] = i+1 in cyclic order 
            // modulo because we never erase endpoints 0 and n-1, so they will always be there to link the cycle together
            prev[i] = (i - 1 + n) % n;
            next[i] = (i + 1) % n;
        }
    }

    // erase vertex v from the cycle (mark as not alive and link its neighbors together)
    void erase(int v) {
        if (!alive[v]) return;
        int a = prev[v];
        int b = next[v];
        next[a] = b;
        prev[b] = a;
        alive[v] = 0;
    }

    // count how many vertices are still alive in the cycle
    int alive_count() const {
        int c = 0;
        for (int i = 0; i < n; ++i) if (alive[i]) ++c;
        return c;
    }
};

//Phase A step: remove min ordinary ear if exists
static bool remove_one_ear(const Precomp& P, BoundaryCycle& B, std::vector<std::string>& word) {
    int best = -1;

    // scan increasing => first found is min(E)
    // we can only remove an ear v if (prev[v], v, next[v]) is a face in the triangulation
    for (int v = 1; v <= B.n - 2; ++v) { // excludes 0 and n-1 by range
        if (!B.alive[v]) continue;
        int a = B.prev[v];
        int b = B.next[v];
        if (is_face(P, (df::vertex_id)a, (df::vertex_id)v, (df::vertex_id)b)) {
            best = v;
            break;
        }
    }

    // if no ear found, return false to indicate failure
    if (best == -1) return false;

    word.push_back(std::to_string(best));
    B.erase(best);
    return true;
}

//full algorithm 
std::vector<std::string> triangulation_to_word(const df::Tri2& tri, int n_boundary, df::vertex_id star_id) {
    Precomp P = precompute(tri, n_boundary, star_id);
    BoundaryCycle B(n_boundary);

    std::vector<std::string> word;

    // PHASE A
    while (remove_one_ear(P, B, word)) {
        // keep removing
    }

    // PHASE B (star collapse once) only if star is present
    if (has_star_vertex(tri, star_id)) {

        // output star once
        word.push_back("*");

        // S = { x alive : star_adj[x] }
        int u = -1, v = -1;
        for (int x = 0; x < n_boundary; ++x) {
            if (!B.alive[x]) continue;
            if (!P.star_adj[x]) continue;
            if (u == -1) u = x;
            v = x;
        }

        // collapse S \ {u,v}, keep min/max, never erase endpoints
        if (u != -1) {
            for (int x = 0; x < n_boundary; ++x) {
                if (!B.alive[x]) continue;
                if (!P.star_adj[x]) continue;
                if (x == u || x == v) continue;
                if (x == 0 || x == n_boundary - 1) continue;
                B.erase(x);
            }
        }
    }

    // PHASE C
    while (B.alive_count() > 2) {
        bool ok = remove_one_ear(P, B, word);
        if (!ok) break; // avoid infinite loop if inconsistent
    }

    return word;
}

std::string word_to_string(const std::vector<std::string>& w) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < w.size(); ++i) {
        if (i) oss << ' ';
        oss << w[i];
    }
    return oss.str();
}

void print_words_for_poset(
    const df::Tri2& tri_start,
    const pst_es::FlipPoset& Pposet,
    int n_boundary,
    df::vertex_id star_id,
    const df::P2& star_point
) {
    for (int idx = 0; idx < (int)Pposet.nodes.size(); ++idx) {
        df::Tri2 tri;
        bool ok = pst_es::reconstruct_triangulation(tri_start, Pposet, idx, star_id, star_point, tri);
        if (!ok) {
            std::cout << "T" << idx << ": <reconstruction failed>\n";
            continue;
        }

        auto w = triangulation_to_word(tri, n_boundary, star_id);
        std::cout << "T" << idx << ": " << word_to_string(w) << "\n";
    }
}

} // namespace earstar
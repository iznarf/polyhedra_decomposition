#pragma once

#include "input.h"

#include <vector>
#include <array>
#include <unordered_map>

namespace pst_es {

// ------------------------------------------------------------
// Signature = sorted list of normalized finite edges
// ------------------------------------------------------------

struct TriSignature {
    std::vector<std::array<df::vertex_id,2>> edges;

    bool operator==(TriSignature const& o) const {
        return edges == o.edges;
    }
};

struct TriSignatureHash {
    std::size_t operator()(TriSignature const& s) const noexcept;
};

// ------------------------------------------------------------
// Move types
// ------------------------------------------------------------

enum class MoveKind {
    Flip,
    InsertStar,
    DeleteStar
};

struct Move {
    MoveKind kind;
    df::vertex_id a = 0;
    df::vertex_id b = 0;
};

// ------------------------------------------------------------

struct Node {
    TriSignature sig;
    int parent = -1;
    Move move_from_parent;
    std::vector<int> children;
};

struct FlipPoset {
    std::vector<Node> nodes;
    std::unordered_map<TriSignature,int,TriSignatureHash> sig_to_idx;
};

// ------------------------------------------------------------

FlipPoset build_flip_poset(
    const df::Tri2& tri_start,
    df::vertex_id star_id,
    const df::P2& star_point
);

bool reconstruct_triangulation(
    const df::Tri2& tri_start,
    const FlipPoset& P,
    int idx,
    df::vertex_id star_id,
    const df::P2& star_point,
    df::Tri2& out_tri
);

} // namespace pst_es
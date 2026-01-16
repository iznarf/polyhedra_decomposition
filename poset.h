#pragma once 
#include "input.h"
#include <vector>
#include <array>
#include <cstddef> 

namespace pst {

    using FaceTriple = std::array<df::vertex_id, 3>;

    struct TriSignature {
        // list of face triples (sorted)
        std::vector<FaceTriple> faces;
    };

    // equality for use in std::unordered_map (signature map)
    inline bool operator==(const TriSignature& a, const TriSignature& b) {
        return a.faces == b.faces;
    }

    // poset node structure
    struct Node {
        std::vector<df::StepRecord> history;  // steps from upper to here
        TriSignature                signature; // triangulation signature at this node
        std::vector<int> parents;          // all poset parents
        std::vector<int> children;         // all poset children
        std::vector<df::StepRecord> child_steps; // steps to each child

    };

    TriSignature make_signature(const df::Tri2& T);

    void replay_history_poset(df::Tri2& tri, const std::vector<df::StepRecord>& history, const df::InputData& D);

    // builds the conforming flip poset from upper to lower triangulation
    void build_poset(const df::InputData& D, std::vector<Node>& nodes);

    std::vector<int> nodes_with_no_incoming_down_flips(const std::vector<Node>& nodes); 
 
    struct PosetTriIndices {
        int upper  = -1;
        int current = -1;
        int lower   = -1;
    };

    PosetTriIndices find_special_triangulations_in_poset(const df::InputData& D,const std::vector<pst::Node>& nodes);

    bool exists_path_via_children(const std::vector<pst::Node>& nodes,int source_idx, int target_idx);

    bool find_conforming_down_path_in_global_poset(
        const df::InputData& D,
        const std::vector<pst::Node>& nodes,
        int source_idx,
        int target_idx,
        std::vector<int>& out_node_path,             // global poset indices (mesh numbers)
        std::vector<df::StepRecord>& out_step_path   // steps along that path
    );
        
    void replay_step_poset(const df::StepRecord& step, df::Tri2& tri, const df::InputData& D);

    void replay_history_poset(df::Tri2& tri, const std::vector<df::StepRecord>& history, const df::InputData& D);


} // namespace pst


// hash function for signatures for use in unordered_map
namespace std {
    template <>
    struct hash<pst::TriSignature> {
        std::size_t operator()(const pst::TriSignature& sig) const noexcept {
            std::size_t h = 0;
            std::hash<df::vertex_id> hv;

            for (const auto& f : sig.faces) {
                std::size_t hf = hv(f[0]);
                hf ^= hv(f[1]) + 0x9e3779b9 + (hf << 6) + (hf >> 2);
                hf ^= hv(f[2]) + 0x9e3779b9 + (hf << 6) + (hf >> 2);

                h ^= hf + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };
}

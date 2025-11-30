#pragma once 
#include "input.h"
#include <vector>
#include <array>
#include <cstddef> // for std::size_t

namespace pst {

    using FaceTriple = std::array<df::vertex_id, 3>;

    struct TriSignature {
        // list of face triples (sorted)
        std::vector<FaceTriple> faces;
    };

    // equality for use in std::map / std::unordered_map
    inline bool operator==(const TriSignature& a, const TriSignature& b) {
        return a.faces == b.faces;
    }

    struct Node {
        std::vector<df::StepRecord> history;  // steps from upper to here
        TriSignature                signature; // triangulation signature at this node
        std::vector<int> parents;          // all poset parents
        std::vector<int> children;         // all poset children

    };

    void build_poset(const df::InputData& D);

} // namespace pst


// hash function for TriSignature for use in unordered_map
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

// we take input data struct, we begin with current triangulation = upper triangulation
// we begin with collecting all possible down-flips. so we try to insert all possible points and we try do apply all possible edge flips
// lets say this gives us a list of edges to flip with k elements
// then we create k lists, each list corresponds to one path in the poset to the lower triangulation
// we also have to create k triangulations for comparison
// best would be to create a edge and face list for every triangulation so we can compare them easily
// we have to sort the edge pairs by their global vertex ids and face triples also by their global vertex ids
// then we can compare them by comparing these lists
// so we do not create CGAL objects, we have one CGAL object which is the upper triangulation and then we record the flips in the lists
// so each list is a different path to the lower regular triangulation
// how to record each flip? 
// step record and flip record structs already exist, also the replay function already exists. we can use this function to create then each triangulation
// for every flip we need to check every triangulation if it is already contained in the poset


#include "input.h"
#include "check_edges.h"
#include "insertion.h"
#include "geometry_utils.h"

#include "poset.h"
#include <algorithm>
#include <unordered_map>

namespace pst {

    TriSignature make_signature(const df::Tri2& T) {
        TriSignature sig;
        sig.faces.clear();

        for (auto f = T.finite_faces_begin(); f != T.finite_faces_end(); ++f) {
            df::vertex_id ids[3] = {
                f->vertex(0)->info(),
                f->vertex(1)->info(),
                f->vertex(2)->info()
            };
            std::sort(std::begin(ids), std::end(ids)); // sort within face
            sig.faces.push_back({ ids[0], ids[1], ids[2] });
        }

        std::sort(sig.faces.begin(), sig.faces.end()); // sort faces globally
        return sig;
    }

    int apply_edge_flip_poset(df::vertex_id ia, df::vertex_id ib, int parent_idx, df::Tri2& tri, std::vector<Node>& nodes, std::unordered_map<TriSignature, int>& sig_to_node) {
        // 1) Find the edge (ia, ib) in the current triangulation by vertex ids
        df::Tri2::Face_handle fh;
        int ei = -1;
        bool found = false;

        for (auto e = tri.finite_edges_begin(); e != tri.finite_edges_end(); ++e) {
            auto f  = e->first;
            int  i  = e->second;

            auto va = f->vertex(tri.cw(i));
            auto vb = f->vertex(tri.ccw(i));

            df::vertex_id ja = va->info();
            df::vertex_id jb = vb->info();

            if ((ja == ia && jb == ib) || (ja == ib && jb == ia)) {
                fh    = f;
                ei    = i;
                found = true;
                break;
            }
        }

        if (!found) {
            std::cerr << "[flip] ERROR: edge (" << ia << "," << ib
                    << ") not found in triangulation\n";
            return -1;
        }

        // 2) boundary / infinite check
        if (tri.is_infinite(fh) || tri.is_infinite(fh->neighbor(ei))) {
            std::cerr << "[flip] WARNING: edge (" << ia << "," << ib
                    << ") is on boundary / infinite, not flipping\n";
            return -1;
        }

        // 3) determine the opposite vertices c,d *before* the flip 
        auto vc = fh->vertex(ei);           // opposite vertex in fh
        auto gh = fh->neighbor(ei);         // neighbor across the edge
        int  j  = tri.mirror_index(fh, ei); // mirrored index in neighbor
        auto vd = gh->vertex(j);            // opposite vertex in gh

        df::vertex_id ic = vc->info();
        df::vertex_id id = vd->info();

        // 4) actually perform the flip
        tri.flip(fh, ei);

        // 5) build the step record for this flip
        df::StepRecord step;
        step.kind = df::StepKind::EdgeFlip;
        step.a = ia;
        step.b = ib;
        step.c = ic;
        step.d = id;

        // 6) compute signature of the new triangulation
        TriSignature sig = make_signature(tri);

        // 7) see if this triangulation already exists as a node
        auto it = sig_to_node.find(sig);
        int child_idx;

        if (it == sig_to_node.end()) {
            // --- new triangulation: create child node ---
            child_idx = static_cast<int>(nodes.size());

            Node child;
            child.history = nodes[parent_idx].history;  // copy parent history
            child.history.push_back(step);              // add this flip
            child.signature = std::move(sig);

            child.parents.push_back(parent_idx);        // poset edge parent -> child

            nodes.push_back(std::move(child));
            sig_to_node[nodes[child_idx].signature] = child_idx;
        } else {
            // --- already seen triangulation ---
            child_idx = it->second;

            // create poset edge parent -> existing child
            nodes[child_idx].parents.push_back(parent_idx);
        }

        // in either case, parent gets this child
        nodes[parent_idx].children.push_back(child_idx);

        return child_idx;
    }



    int apply_vertex_insertion_poset(df::vertex_id v, int parent_idx, df::Tri2& tri, const df::InputData& D, std::vector<Node>& nodes, std::unordered_map<TriSignature, int>& sig_to_node) {
        // 1) 2D position of the vertex to be inserted
        const df::P2& p2 = D.points2d[v];

        // 2) locate in current triangulation
        df::Tri2::Locate_type lt;
        int li = -1;
        df::Tri2::Face_handle fh = tri.locate(p2, lt, li);

        if (lt != df::Tri2::FACE) {
            std::cerr << "[poset] WARNING: insertion point " << v
                    << " is not strictly inside a face, skipping\n";
            return -1;
        }

        // 3) face containing the point BEFORE insertion
        auto va = fh->vertex(0);
        auto vb = fh->vertex(1);
        auto vc = fh->vertex(2);

        df::vertex_id ia = va->info();
        df::vertex_id ib = vb->info();
        df::vertex_id ic = vc->info();

        // 4) actually insert the vertex in this face
        df::Tri2::Vertex_handle vh_new = tri.insert_in_face(p2, fh);
        vh_new->info() = v; // set global id

        // 5) build the step record
        df::StepRecord step;
        step.kind = df::StepKind::VertexInsertion;
        step.a = ia;
        step.b = ib;
        step.c = ic;
        step.d = v;   // the inserted vertex

        // 6) compute signature of the new triangulation
        TriSignature sig = make_signature(tri);

        // 7) see if this triangulation already exists as a node
        auto it = sig_to_node.find(sig);
        int child_idx;

        if (it == sig_to_node.end()) {
            // --- new triangulation: create child node ---
            child_idx = static_cast<int>(nodes.size());

            Node child;
            child.history = nodes[parent_idx].history;  // copy parent history
            child.history.push_back(step);              // add this insertion
            child.signature = std::move(sig);

            child.parents.push_back(parent_idx);

            nodes.push_back(std::move(child));
            sig_to_node[nodes[child_idx].signature] = child_idx;
        } else {
            // --- already seen triangulation ---
            child_idx = it->second;
            nodes[child_idx].parents.push_back(parent_idx);
        }

        nodes[parent_idx].children.push_back(child_idx);
        return child_idx;
    }

    static const char* step_kind_name(df::StepKind k) {
        switch (k) {
            case df::StepKind::EdgeFlip:         return "EdgeFlip";
            case df::StepKind::VertexInsertion:  return "VertexInsertion";
        
        }
    }

    void debug_print_poset(const std::vector<Node>& nodes) {
        std::cout << "[poset] total nodes: " << nodes.size() << "\n";

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            const Node& n = nodes[i];

            std::cout << "Node " << i
                    << " | parents: "  << n.parents.size()
                    << ", children: "   << n.children.size()
                    << ", history size: " << n.history.size()
                    << "\n";

            // print the step history for this node
            for (std::size_t s = 0; s < n.history.size(); ++s) {
                const df::StepRecord& step = n.history[s];

                std::cout << "  step " << s
                        << " : " << step_kind_name(step.kind)
                        << " (a=" << step.a
                        << ", b=" << step.b
                        << ", c=" << step.c
                        << ", d=" << step.d
                        << ")\n";
            }

        std::cout << "\n";
        }
    }








    void build_poset(const df::InputData& D) {
        // we strart with upper triangulation
        df::Tri2 tri_root = D.tri_poset;
        df::Tri2 tri_lower = D.tri_lower;

        // create root node
        std::vector<Node> nodes;
        nodes.reserve(300); // arbitrary

        Node root;
        root.history.clear();  // no steps from upper to upper
        root.signature = make_signature(tri_root);

        nodes.push_back(root);

        // map from signature to node index
        // does this have to be a global variable?
        std::unordered_map<TriSignature, int> sig_to_node;
        sig_to_node.emplace(root.signature, 0);


        // collect all possible down-flips
        // first edge flips
        std::vector<std::array<df::vertex_id,2>> down_flip_edges = df::reg::find_locally_non_regular_edges(tri_root);
        // then vertex insertions
        std::vector<df::vertex_id> missing_vertices = df::find_missing_vertices(tri_root, tri_lower);

        // for every down-flip we create a new triangulation which is a node in the poset
        for (const auto& edge : down_flip_edges) {
            tri_root = D.tri_poset; // reset to root triangulation
            df::vertex_id ia = edge[0];
            df::vertex_id ib = edge[1];
            int child_idx = apply_edge_flip_poset(ia, ib, /*parent_idx=*/0,
            tri_root, nodes, sig_to_node);
            // apply every edge flip in the list to the current triangulation
            // create a node for every new triangulation
            // create signature for each new triangulation
            // at this point we do not have to check if signature already exists since this are the first triangulations after the root
            // add the flip to the history of the node
            // update parents (root) 
            // update children of root
        }

        // we have to do the same for every missing vertex (insertion)
        for (df::vertex_id v : missing_vertices) {
            df::Tri2 tri = tri_root; // start from upper triangulation
            int child_idx = apply_vertex_insertion_poset(v, /*parent_idx=*/0, tri, D, nodes, sig_to_node);
        }

        // after the edge + insertion loops
        std::cout << "[poset] after first layer:\n";
        debug_print_poset(nodes);

    }
    
}
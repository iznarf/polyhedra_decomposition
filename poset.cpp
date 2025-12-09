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
#include "conforming.h"
#include "conforming_insertion.h"

#include "poset.h"
#include <algorithm>
#include <unordered_map>

namespace pst {
    // signature of a triangulation: sorted list of face triples (each triple sorted)
    // we have this so that we can compare triangulations 
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
            // print the sig.faces
        }

        std::sort(sig.faces.begin(), sig.faces.end()); // sort faces globally
        // print the full signature
       
        return sig;
    }

    // apply an edge flip (ia, ib) to triangulation 'tri' at node 'parent_idx' in the poset
    // checks if the resulting triangulation already exists in the poset
    // if not, creates a new node
    // if yes, just creates a new poset edge
    int apply_edge_flip_poset(df::vertex_id ia, df::vertex_id ib,
                          int current_idx,
                          df::Tri2& tri,
                          std::vector<Node>& nodes,
                          std::unordered_map<TriSignature, int>& sig_to_node){
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

    // 3) extract the full quad (a,b,c,d) around this edge BEFORE the flip
    auto va = fh->vertex(tri.cw(ei));      // endpoint a
    auto vb = fh->vertex(tri.ccw(ei));     // endpoint b
    auto vc = fh->vertex(ei);              // opposite c in f
    auto gh = fh->neighbor(ei);
    int  j  = tri.mirror_index(fh, ei);
    auto vd = gh->vertex(j);               // opposite d in neighbor

    df::vertex_id a_id = va->info();
    df::vertex_id b_id = vb->info();
    df::vertex_id c_id = vc->info();
    df::vertex_id d_id = vd->info();

    const df::P2& a2 = va->point();
    const df::P2& b2 = vb->point();
    const df::P2& c2 = vc->point();
    const df::P2& d2 = vd->point();

    df::P3 a3 = df::lift(a2);
    df::P3 b3 = df::lift(b2);
    df::P3 c3 = df::lift(c2);
    df::P3 d3 = df::lift(d2);

    // 4) decide local up/down orientation using the same convention as insertion
    CGAL::Orientation orient = df::oriented_height_sign(a2, b2, c2,a3, b3, c3, d3);

    // 5) perform the flip
    tri.flip(fh, ei);

    // 6) build the step record (a,b) edge and (c,d) opposite vertices
    df::StepRecord step;
    if (orient == CGAL::NEGATIVE) {
        // print down flip
        std::cout << "[flip] DOWN flip\n";
        // DOWN flip
        step.kind = df::StepKind::EdgeFlip_down;
    } else if (orient == CGAL::POSITIVE) {
        // UP flip
        std::cout << "[flip] UP flip\n";
        step.kind = df::StepKind::EdgeFlip_up;
    }
    
    step.a = a_id;
    step.b = b_id;
    step.c = c_id;
    step.d = d_id;

    // 7) compute signature of the new triangulation
    TriSignature sig = make_signature(tri);

    
    // 8) dedup / create node
    auto it = sig_to_node.find(sig);
    int idx;

    if (orient == CGAL::NEGATIVE) {
        // DOWN flip: new node is a child of current_idx
        // we have not seen this triangulation before -> create new node
        if (it == sig_to_node.end()) {

            // print that we create a new node
            std::cout << "[flip] Creating new node\n";
            idx = static_cast<int>(nodes.size());

            Node child;
            child.history = nodes[current_idx].history;
            child.history.push_back(step);
            child.signature = std::move(sig);
            child.parents.push_back(current_idx);  // parent -> child
            nodes.push_back(std::move(child));
            sig_to_node[nodes[idx].signature] = idx;
        // we already know the triangulation and found a new path to it -> update parents 
        // we have to check if the parent is already contained 
        } else {
            std::cout << "[flip] Reusing existing node\n";
            idx = it->second;
            nodes[idx].parents.push_back(current_idx); 
        }
        
        nodes[current_idx].children.push_back(idx);
        nodes[current_idx].child_steps.push_back(step);
       

    } else if (orient == CGAL::POSITIVE) {
        // UP flip: new node is a *parent* of current_idx
        if (it == sig_to_node.end()) {
            idx = static_cast<int>(nodes.size());

            Node parent;
            parent.history = nodes[current_idx].history;
            parent.history.push_back(step);
            parent.signature = std::move(sig);

            //parent.children.push_back(current_idx); // parent -> current_idx
             

            nodes.push_back(std::move(parent));
            sig_to_node[nodes[idx].signature] = idx;
        } else {
            std::cout << "[flip] Reusing existing node\n";
            idx = it->second;
            //nodes[idx].children.push_back(current_idx); // existing parent
           
        }

        //nodes[current_idx].parents.push_back(idx);
    } else {
        // COLLINEAR or degenerate (should not happen with generic input)
        std::cerr << "[flip] WARNING: oriented_height_sign == COLLINEAR, skipping\n";
        return -1;
    }

    return idx;
}


    // apply a vertex insertion of vertex v to triangulation tri at node parent_idx in the poset
    // checks if the resulting triangulation already exists in the poset
    // if not, creates a new node
    // if yes, just creates a new poset edge
    int apply_vertex_insertion_poset(df::vertex_id v, int current_idx, df::Tri2& tri, const df::InputData& D, std::vector<Node>& nodes, std::unordered_map<TriSignature, int>& sig_to_node) {
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


        // 5) check if insertion is up or downflip
        const df::P2& a2 = va->point();
        const df::P2& b2 = vb->point();
        const df::P2& c2 = vc->point();
        const df::P2& d2 = vh_new->point();

        df::P3 a3 = df::lift(a2);
        df::P3 b3 = df::lift(b2);
        df::P3 c3 = df::lift(c2);
        df::P3 d3 = df::lift(d2);

        CGAL::Orientation orient = df::oriented_height_sign(a2, b2, c2, a3, b3, c3, d3);

        // 6) build the step record
        df::StepRecord step;
        if (orient == CGAL::NEGATIVE) {
            std::cout << "[insertion] DOWN insertion\n";
            // DOWN insertion
            step.kind = df::StepKind::VertexInsertion_down;
        } else if (orient == CGAL::POSITIVE) {
            std::cout << "[insertion] UP insertion\n";
            // UP insertion
            step.kind = df::StepKind::VertexInsertion_up;
        }
        

        step.a = ia;
        step.b = ib;
        step.c = ic;
        step.d = v;   // the inserted vertex

        // 7) compute signature of the new triangulation
        TriSignature sig = make_signature(tri);

        // 8) dedup / create node
        auto it = sig_to_node.find(sig);
        int idx;

        if (orient == CGAL::NEGATIVE) {
            // DOWN insertion: new node is a child of current_idx
            // we have not seen this triangulation before -> create new node
            if (it == sig_to_node.end()) {

                // print that we create a new node
                std::cout << "[insertion] Creating new node\n";
                idx = static_cast<int>(nodes.size());

                Node child;
                child.history = nodes[current_idx].history;
                child.history.push_back(step);
                child.signature = std::move(sig);
                
                child.parents.push_back(current_idx);  // parent -> child
               

                nodes.push_back(std::move(child));
                sig_to_node[nodes[idx].signature] = idx;

            // we already know the triangulation and found a new path to it -> update parents 
            // we have to check if the parent is already contained 
            } else {
                std::cout << "[insertion] Reusing existing node\n";
                idx = it->second;
                nodes[idx].parents.push_back(current_idx); 
            }
        
            nodes[current_idx].children.push_back(idx);
            nodes[current_idx].child_steps.push_back(step);
           

            } else if (orient == CGAL::POSITIVE) {
                // UP insertion: new node is a *parent* of current_idx
                if (it == sig_to_node.end()) {
                    idx = static_cast<int>(nodes.size());

                    Node parent;
                    parent.history = nodes[current_idx].history;
                    parent.history.push_back(step);
                    parent.signature = std::move(sig);

                    //parent.children.push_back(current_idx); // parent -> current_idx
                  

                    nodes.push_back(std::move(parent));
                    sig_to_node[nodes[idx].signature] = idx;
                } else {
                    std::cout << "[insertion] Reusing existing node\n";
                    idx = it->second;
                    //nodes[idx].children.push_back(current_idx); // existing parent
                  
                }

                //nodes[current_idx].parents.push_back(idx);
        } else {
            // COLLINEAR or degenerate (should not happen with generic input)
            std::cerr << "[insertion] WARNING: oriented_height_sign == COLLINEAR, skipping\n";
            return -1;
        }
        return idx;
    }


    // apply a vertex deletion of vertex v to triangulation tri at node parent_idx in the poset
    // checks if the resulting triangulation already exists in the poset

      int apply_vertex_deletion_poset(df::vertex_id vid, int current_idx, df::Tri2& tri, const df::InputData& D, std::vector<Node>& nodes, std::unordered_map<TriSignature, int>& sig_to_node) {
            // 1) find vertex handle
            df::Tri2::Vertex_handle vh = nullptr;
            for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit) {
                if (vit->info() == vid) {
                    vh = vit;
                    break;
                }
            }
            if (vh == nullptr) {
                    std::cerr << "[deletion] ERROR: vertex " << vid << " not found in triangulation\n";
                    return -1;
            }

            // 2) collect 3 finite neighbors (we know degree == 3 and not boundary)
            std::array<df::Tri2::Vertex_handle, 3> neighbor_vertices;
            {
                auto vc = tri.incident_vertices(vh);
                auto start = vc;
                int k = 0;

                do {
                    if (!tri.is_infinite(vc)) {
                        if (k < 3) neighbor_vertices[k] = vc;
                        ++k;
                    }
                    ++vc;
                } while (vc != start);

                if (k != 3) {
                    std::cerr << "[deletion] vertex " << vid << " does not have 3 finite neighbors\n";
                    return -1;
                }
            }

            auto va = neighbor_vertices[0];
            auto vb = neighbor_vertices[1];
            auto vc = neighbor_vertices[2];

            df::vertex_id ia = va->info();
            df::vertex_id ib = vb->info();
            df::vertex_id ic = vc->info();

            // 3) collect points *before* deletion
            const df::P2& a2 = va->point();
            const df::P2& b2 = vb->point();
            const df::P2& c2 = vc->point();
            const df::P2& d2 = vh->point(); // to be deleted

            df::P3 a3 = df::lift(a2);
            df::P3 b3 = df::lift(b2);
            df::P3 c3 = df::lift(c2);
            df::P3 d3 = df::lift(d2);

            // the orientation here was not consistent -> fix this properly, not just randomly swap b and a!!!

            CGAL::Orientation orient = df::oriented_height_sign(a2, b2, c2, a3, b3, c3, d3);

            if (orient == CGAL::COLLINEAR) {
                std::cerr << "[deletion] WARNING: oriented_height_sign == COLLINEAR, skipping\n";
                return -1;
            }

            // 4) delete the vertex
            tri.remove_degree_3(vh); 

            // 6) build step record
            df::StepRecord step;
            if (orient == CGAL::NEGATIVE) {
                // print down deletion
                std::cout << "[deletion] UP deletion\n";
                // if d lies below the plane and we delete d, then we make an up flip
                step.kind = df::StepKind::VertexDeletion_up;
            } else if (orient == CGAL::POSITIVE) {
                // if d lies above the plane and we delete d, then we make a down flip
                std::cout << "[deletion] DOWN deletion\n";
                step.kind = df::StepKind::VertexDeletion_down;
            }
           
          
           
           
            step.a = ia;
            step.b = ib;
            step.c = ic;
            step.d = vid;   // deleted vertex id

            // 5) compute signature of new triangulation
            TriSignature sig = make_signature(tri);

            // 7) dedup / create node
            auto it = sig_to_node.find(sig);
            int idx;

            if (orient == CGAL::NEGATIVE) {
                // DOWN deletion: new node is child of current_idx
                if (it == sig_to_node.end()) {
                    std::cout << "[deletion] Creating new node\n";
                    idx = static_cast<int>(nodes.size());

                    Node child;
                    child.history = nodes[current_idx].history;
                    child.history.push_back(step);
                    child.signature = std::move(sig);
                    child.parents.push_back(current_idx);

                    nodes.push_back(std::move(child));
                    sig_to_node[nodes[idx].signature] = idx;
                } else {
                    std::cout << "[deletion] Reusing existing node\n";
                    idx = it->second;
                    nodes[idx].parents.push_back(current_idx);
                }
                nodes[current_idx].children.push_back(idx);
                nodes[current_idx].child_steps.push_back(step);

            } if (orient == CGAL::POSITIVE) { 
                // UP deletion: new node is *parent* of current_idx
                if (it == sig_to_node.end()) {
                    std::cout << "[deletion] Creating new node\n";
                    idx = static_cast<int>(nodes.size());

                    Node parent;
                    parent.history = nodes[current_idx].history;
                    parent.history.push_back(step);
                    parent.signature = std::move(sig);
                    //parent.children.push_back(current_idx);
                   

                    nodes.push_back(std::move(parent));
                    sig_to_node[nodes[idx].signature] = idx;
                } else {
                    std::cout << "[deletion] Reusing existing node\n";
                    idx = it->second;
                    //nodes[idx].children.push_back(current_idx);
                   
                }
                //nodes[current_idx].parents.push_back(idx);

            }

            return idx;
        }



    // replay a single step (edge flip or vertex insertion) on triangulation tri
    // we need this to reconstruct triangulations at poset nodes from step histories
    void replay_step_poset(const df::StepRecord& step, df::Tri2& tri, const df::InputData& D) {
        if (step.kind == df::StepKind::EdgeFlip_down || step.kind == df::StepKind::EdgeFlip_up) { // EdgeFlip
            df::vertex_id ia = step.a;
            df::vertex_id ib = step.b;

            df::Tri2::Face_handle fh;
            int ei = -1;
            bool found = false;

            for (auto e = tri.finite_edges_begin(); e != tri.finite_edges_end(); ++e) {
                auto f = e->first;
                int  i = e->second;

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
                std::cerr << "[replay] edge (" << ia << "," << ib
                        << ") not found during replay\n";
                return;
            }
            if (tri.is_infinite(fh) || tri.is_infinite(fh->neighbor(ei))) {
                std::cerr << "[replay] edge (" << ia << "," << ib
                        << ") is boundary / infinite during replay\n";
                return;
            }

            tri.flip(fh, ei);
        } 
        if (step.kind == df::StepKind::VertexInsertion_down || step.kind == df::StepKind::VertexInsertion_up) { // VertexInsertion
            df::vertex_id v  = step.d;
            const df::P2& p2 = D.points2d[v];

            df::Tri2::Locate_type lt;
            int li = -1;
            df::Tri2::Face_handle fh = tri.locate(p2, lt, li);

            
            if (lt != df::Tri2::FACE) {
                std::cerr << "[replay] WARNING: insertion point " << v
                        << " not in FACE during replay\n";
                return;
            }
            df::Tri2::Vertex_handle vh = tri.insert_in_face(p2, fh);
            vh->info() = v;
        }
        if (step.kind == df::StepKind::VertexDeletion_down || step.kind == df::StepKind::VertexDeletion_up) { // VertexDeletion
            df::vertex_id vid = step.d;

            // find vertex handle
            df::Tri2::Vertex_handle vh = nullptr;
            for (auto vit = tri.finite_vertices_begin(); vit != tri.finite_vertices_end(); ++vit) {
                if (vit->info() == vid) {
                    vh = vit;
                    break;
                }
            }
            if (vh == nullptr) {
                std::cerr << "[replay] ERROR: vertex " << vid << " not found during replay\n";
                return;
            }

            tri.remove_degree_3(vh); 
        }
    }

    // replay a sequence of steps on triangulation tri
    void replay_history_poset(df::Tri2& tri, const std::vector<df::StepRecord>& history, const df::InputData& D) {
        for (const auto& step : history) {
            replay_step_poset(step, tri, D);
        }
    }


    std::vector<std::array<df::vertex_id, 2>> find_flip_edges(const df::Tri2& tri) {
        std::vector<std::array<df::vertex_id, 2>> flippable_edges;

        for (auto edge = tri.finite_edges_begin(); edge != tri.finite_edges_end(); ++edge) {
            auto f = edge->first;  // incident face of the edge
            int  i = edge->second; // index of opposite vertex in face f

            // internal edge? both incident faces must be finite
            if (tri.is_infinite(f)) continue;
            auto g = f->neighbor(i);
            if (tri.is_infinite(g)) continue;

            // extract the quad around the edge
            auto va = f->vertex(tri.cw(i));
            auto vb = f->vertex(tri.ccw(i));
            auto vc = f->vertex(i);
            int  j  = tri.mirror_index(f, i);
            auto vd = g->vertex(j);

            const df::P2 &a2 = va->point(), &b2 = vb->point(), &c2 = vc->point(), &d2 = vd->point();

            if (!df::quad_strictly_convex(a2, b2, c2, d2)) continue;

            // check orientation of lifted points
            df::P3 a3 = df::lift(a2);
            df::P3 b3 = df::lift(b2);
            df::P3 c3 = df::lift(c2);
            df::P3 d3 = df::lift(d2);

            CGAL::Orientation orient = df::oriented_height_sign(a2, b2, c2,
                                                                a3, b3, c3, d3);

            if (orient == CGAL::COLLINEAR) continue; // degenerate, skip



            df::vertex_id ia = va->info();
            df::vertex_id ib = vb->info();
    
            flippable_edges.push_back({ ia, ib });
        }
        return flippable_edges;
    }

    bool is_boundary_vertex(const df::Tri2& tri, df::Tri2::Vertex_handle v) {
        if (tri.is_infinite(v))
            return false;

        auto fc = tri.incident_faces(v);
        if (fc == nullptr)
            return false;

        auto start = fc;

        while (true) {
            if (tri.is_infinite(fc))
                return true;

            ++fc;
            if (fc == start) break;
        }
        return false;
    }


    std::vector<df::vertex_id>
    find_deletion_vertices(const df::Tri2& current) {
        std::vector<df::vertex_id> deletion_vertices;
        for (auto v = current.finite_vertices_begin(); v != current.finite_vertices_end(); ++v) {
            if (is_boundary_vertex(current, v))
                continue;
            // check if v has degree three
            // if yes, then add it to the deletion vertices
            if (current.degree(v) == 3)
                deletion_vertices.push_back(v->info());
        }
        return deletion_vertices;
    }

    // builds the conforming flip poset from upper to lower triangulation
    // can easily be adapted to build full poset including non-conforming flips/insertions
    void build_poset(const df::InputData& D, std::vector<Node>& nodes) {
        df::Tri2 tri_root  = D.tri_poset;
        //df::Tri2 tri_root  = D.tri_current;
        df::Tri2 tri_lower = D.tri_lower;

        // clear any old content and reserve some space
        nodes.clear();
        nodes.reserve(600); // arbitrary initial guess

        Node root;
        root.history.clear();
        root.signature = make_signature(tri_root);
        nodes.push_back(root);

        std::unordered_map<TriSignature, int> sig_to_node;
        sig_to_node.emplace(root.signature, 0);

    
        std::size_t current_idx = 0; // start at root

        while (current_idx < nodes.size()) {
                /*
                if (current_idx > 50) {
                    break; // for testing: limit number of nodes
                }
                */

                Node& current_node = nodes[current_idx];

                // reconstruct triangulation for this node
                df::Tri2 tri = D.tri_poset;
                //df::Tri2 tri = D.tri_current;
                replay_history_poset(tri, current_node.history, D);


                // find all flippable edges -> up OR down to get full poset
                auto flip_edges = find_flip_edges(tri);
                auto missing_vertices = df::find_missing_vertices(tri, tri_lower);
                auto deletion_vertices = find_deletion_vertices(tri);

                // try all possible edge flips from current triangulation
                for (const auto& edge : flip_edges) {
                    df::Tri2 tri_child = tri;
                    df::vertex_id ia = edge[0];
                    df::vertex_id ib = edge[1];

                   
                    apply_edge_flip_poset(ia, ib, current_idx, tri_child,
                                        nodes, sig_to_node);
                }

                // try all possible vertex insertions from current triangulation
                for (df::vertex_id v : missing_vertices) {
                    df::Tri2 tri_child = tri;
                    apply_vertex_insertion_poset(v, current_idx, tri_child,
                                                D, nodes, sig_to_node);
                }

                // try all possible vertex deletions from current triangulation
                for (df::vertex_id v : deletion_vertices) {
                    df::Tri2 tri_child = tri;
                    apply_vertex_deletion_poset(v, current_idx, tri_child, D, nodes, sig_to_node);
                }

                ++current_idx;
               

        }

        std::cout << "[poset] finished building flip poset\n";
        /*
        //print signature of every node
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            std::cout << "Node " << i << " signature: ";
            for (const auto& face : nodes[i].signature.faces) {
                std::cout << "(" << face[0] << "," << face[1] << "," << face[2] << ") ";
            }
            std::cout << std::endl;
        }
        */

        //print step history of every node
        for (std::size_t i = 0; i < nodes.size(); ++i) {
            std::cout << "Node " << i << " history: ";
            for (const auto& step : nodes[i].history) {
                if (step.kind == df::StepKind::EdgeFlip_down || step.kind == df::StepKind::EdgeFlip_up) {
                    std::cout << "Flip(" << step.a << "," << step.b << ") ";
                } else if (step.kind == df::StepKind::VertexInsertion_down || step.kind == df::StepKind::VertexInsertion_up) {
                    std::cout << "Insert(" << step.d << ") ";
                } else if (step.kind == df::StepKind::VertexDeletion_down || step.kind == df::StepKind::VertexDeletion_up) {
                    std::cout << "Delete(" << step.d << ") ";
                }
            }
            std::cout << std::endl;
        }
        

    }


    // DFS for a conforming path from upper to lower triangulation
    // returns true if found, false if not found or aborted (max nodes/depth)
    // does not build the full poset, only explores reachable conforming nodes
    bool find_conforming_path_dfs(const df::InputData& D, std::vector<df::StepRecord>& out_history, std::size_t max_nodes, std::size_t max_depth) {
        // start from upper triangulation
        df::Tri2 tri_root  = D.tri_poset;
        df::Tri2 tri_lower = D.tri_lower;

        // signature of the lower triangulation
        TriSignature sig_lower = make_signature(tri_lower);

        // poset nodes: we will NOT build the full poset, enough to explore paths
        std::vector<Node> nodes;
        nodes.reserve(300);

        Node root;
        root.history.clear();              // no steps from upper to upper
        root.signature = make_signature(tri_root);
        nodes.push_back(root);

        // map signatures to node indices so we don't create duplicates
        std::unordered_map<TriSignature, int> sig_to_node;
        sig_to_node.emplace(root.signature, 0);

        // stack for DFS storing node indices
        std::vector<int> stack;
        stack.push_back(0); // start at root

        // track which indices we have "expanded"
        // but with this pattern each node is expanded at most once
        while (!stack.empty()) {

            int idx = stack.back();
            stack.pop_back();

            if (idx < 0 || idx >= (int)nodes.size())
                continue;

            Node& current_node = nodes[idx];

            // depth = number of flips/insertions from upper
            if (current_node.history.size() > max_depth) {
                continue; // don't expand beyond this depth
            }

            // if this node already is the lower triangulation we are done
            if (current_node.signature == sig_lower) {
                out_history = current_node.history;
                std::cout << "[dfs] found lower triangulation at node " << idx
                        << ", depth " << current_node.history.size() << "\n";
                return true;
            }

            // reconstruct triangulation for this node
            df::Tri2 tri = D.tri_poset;
            replay_history_poset(tri, current_node.history, D);

            // collect candidate moves from this triangulation
            auto down_flip_edges  = df::reg::find_locally_non_regular_edges(tri);
            auto missing_vertices = df::find_missing_vertices(tri, tri_lower);

            // conforming edge flip children
            for (const auto& edge : down_flip_edges) {
                if (nodes.size() >= max_nodes) {
                    std::cout << "[dfs] reached max_nodes = " << max_nodes << "\n";
                    return false; // aborted
                }

                df::vertex_id ia = edge[0];
                df::vertex_id ib = edge[1];

                // only explore conforming flips (your new check)
                if (!df::reg::is_flip_conforming(ia, ib, D, tri))
                    continue;

                df::Tri2 tri_child = tri;
                int child_idx = apply_edge_flip_poset(ia, ib, idx, tri_child, nodes, sig_to_node);

                if (child_idx < 0)
                    continue;

                // if this triangulation is new, push it on the stack
                if ((std::size_t)child_idx == nodes.size() - 1) {
                    stack.push_back(child_idx);
                }

                // early exit if child is lower
                if (nodes[child_idx].signature == sig_lower) {
                    out_history = nodes[child_idx].history;
                    std::cout << "[dfs] found lower via edge flip at node "
                            << child_idx << "\n";
                    return true;
                }
            }

            // conforming vertex insertion children
            for (df::vertex_id v : missing_vertices) {
                if (nodes.size() >= max_nodes) {
                    std::cout << "[dfs] reached max_nodes = " << max_nodes << "\n";
                    return false; // aborted
                }

                // 1) check down-flip condition
                if (!df::is_insertion_downflip(v, D, tri))
                    continue;

                // 2) check global conformance w.r.t. lower triangulation
                if (!df::reg::is_insertion_conforming(v, D, tri))
                    continue;

                df::Tri2 tri_child = tri;
                int child_idx = apply_vertex_insertion_poset(v, idx, tri_child, D, nodes, sig_to_node);
                
                if (child_idx < 0)
                    continue;

                if ((std::size_t)child_idx == nodes.size() - 1) {
                    stack.push_back(child_idx);
                }

                if (nodes[child_idx].signature == sig_lower) {
                    out_history = nodes[child_idx].history;
                    std::cout << "[dfs] found lower via insertion at node "
                            << child_idx << "\n";
                    return true;
                }
            }
        }

    // DFS exhausted all reachable conforming nodes from upper, never hit lower
    std::cout << "[dfs] no conforming path from upper to lower found.\n";
    return false;
}

  
}
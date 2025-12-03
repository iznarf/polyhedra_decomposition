// check the present global vertex ids of the current triangulation and compare them with the target triangulation
// if the target triangulation contains more vertices, then we have to do a vertex insertion
#include "insertion.h"
#include "input.h"
#include "geometry_utils.h"
#include <iostream>

#include <unordered_set>

namespace df {


    // helper function: build a set of global ids present in a triangulation
    static std::unordered_set<df::vertex_id>
    collect_ids(const df::Tri2& tri)
    {
        std::unordered_set<df::vertex_id> ids;
        ids.reserve(tri.number_of_vertices());
        for (auto v = tri.finite_vertices_begin();
             v != tri.finite_vertices_end(); ++v)
        {
            ids.insert(v->info()); // info() == global vertex_id
        }
        return ids;
    }

    std::vector<df::vertex_id>
    find_missing_vertices(const df::Tri2& current,
                          const df::Tri2& target)
    {
        // all ids in current triangulation
        auto current_ids = collect_ids(current);

        // all ids that exist in target but NOT in current
        std::vector<df::vertex_id> missing;
        missing.reserve(target.number_of_vertices());

        for (auto v = target.finite_vertices_begin();
             v != target.finite_vertices_end(); ++v)
        {
            df::vertex_id id = v->info();
            if (current_ids.find(id) == current_ids.end()) {
                missing.push_back(id);
            }
        }

        return missing;
    }

    // sort missing vertices by descending lifted height 
    std::vector<df::vertex_id>
    sorted_insertion_vertices(const std::vector<df::vertex_id>& missing,
                                const df::InputData& D)
    {
        std::vector<std::pair<double, df::vertex_id>> temp;
        temp.reserve(missing.size());

        // compute lifted height and store (height, id)
        for (auto id : missing) {
            const auto& p = D.points2d[id];
            double z = p.x()*p.x() + p.y()*p.y();  // lift 
            temp.emplace_back(z, id);
        }

        // sort by descending height
        std::sort(temp.begin(), temp.end(),
                [](const auto& a, const auto& b) {
                    return a.first > b.first; // largest height first
                });

        // extract ids only
        std::vector<df::vertex_id> sorted_ids;
        sorted_ids.reserve(temp.size());
        for (auto& pair : temp)
            sorted_ids.push_back(pair.second);

        return sorted_ids;
    }

   

    void apply_vertex_insertion(df::vertex_id id, df::InputData& D) {

        const auto& p = D.points2d[id];
        df::Tri2& T = D.tri_current;

        df::Tri2::Locate_type lt;
        int li; // index used if the point lies on an edge
        df::Tri2::Face_handle fh = T.locate(p, lt, li);

        if (lt == df::Tri2::FACE) {
            // p is strictly inside this triangle
            std::cout << "[locate] point is inside a face\n";

          
            auto va = fh->vertex(0);
            auto vb = fh->vertex(1);
            auto vc = fh->vertex(2);

            df::vertex_id ia = va->info();
            df::vertex_id ib = vb->info();
            df::vertex_id ic = vc->info();
           

            // record flip
            df::StepRecord s;
            s.kind = df::StepKind::VertexInsertion;
            s.a = ia; s.b = ib; s.c = ic; s.d = id;  // face (a,b,c) and new vertex d
            D.step_history.push_back(s);

            df::Tri2::Vertex_handle vh = T.insert_in_face(p, fh);
            vh->info() = id; // set the global vertex id 

            std::cout << "Inserted vertex id " << id
                    << " at point (" << p.x() << "," << p.y() << ")\n";
        }
    }

    bool is_insertion_downflip(vertex_id id, const InputData& D) {
        const Tri2& tri = D.tri_current;
        const P2& d2    = D.points2d[id];

        Tri2::Locate_type lt;
        int li;
        Tri2::Face_handle fh = tri.locate(d2, lt, li);

        if (lt != Tri2::FACE) {
            return false;
        }

        auto va = fh->vertex(0);
        auto vb = fh->vertex(1);
        auto vc = fh->vertex(2);

        const P2& a2 = va->point();
        const P2& b2 = vb->point();
        const P2& c2 = vc->point();

        // lift four points 
        P3 a3 = df::lift(a2);
        P3 b3 = df::lift(b2);
        P3 c3 = df::lift(c2);
        P3 d3 = df::lift(d2);

        CGAL::Orientation s = oriented_height_sign(a2, b2, c2, a3, b3, c3, d3);

        // NEGATIVE = d' is below the face -> allowed
        if (s == CGAL::NEGATIVE) {
            return true;
        }
        else {
            return false;
        }
    }

} // namespace df

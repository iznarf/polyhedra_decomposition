#include "flip.h"
#include "input.h"
#include <iostream>

namespace df {

// we now have to perform the edge flip on the current triangulation with given indices
// after the flip the vertex handles remain the same, we do not have to update the mapping from global index to vertex handles 
// since the flip only changes the connectivity of the triangulation, not the vertices themselves



void apply_edge_flip(df::vertex_id ia,
                     df::vertex_id ib,
                     df::InputData& D)
{
    df::Tri2& tri = D.tri_current;

    // Find the edge (ia, ib) in the current triangulation by vertex ids
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
                  << ") not found in tri_current\n";
        return; // or assert(false);
    }

    // boundary / infinite check
    if (tri.is_infinite(fh) || tri.is_infinite(fh->neighbor(ei))) {
        std::cerr << "[flip] WARNING: edge (" << ia << "," << ib
                  << ") is on boundary / infinite, not flipping\n";
        return;
    }

    // determine the opposite vertices c,d *before* the flip 
    auto vc = fh->vertex(ei);                 // opposite vertex in fh
    auto gh = fh->neighbor(ei);               // neighbor across the edge
    int  j  = tri.mirror_index(fh, ei);       // mirrored index in neighbor
    auto vd = gh->vertex(j);                  // opposite vertex in gh

    df::vertex_id ic = vc->info();
    df::vertex_id id = vd->info();

   
    D.flip_history.emplace_back(ia, ib, ic, id);
    df::StepRecord s;
    s.kind = df::StepKind::EdgeFlip;
    s.a = ia; s.b = ib;
    s.c = ic; s.d = id;
    D.step_history.push_back(s);

    
    tri.flip(fh, ei);
}



void print_flip_history(const df::InputData& D)
{
    std::cout << "\n=== edge flip history ===\n";
    int i = 0;
    for (const auto& rec : D.flip_history) {
        std::cout << i++ << ": (" 
                  << rec.a << "," << rec.b 
                  << ") -> (" 
                  << rec.c << "," << rec.d << ")\n";
    }
}



} // namespace df
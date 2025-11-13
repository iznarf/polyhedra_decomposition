#include "flip.h"




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

    // Optional: check we’re not on the boundary
    if (tri.is_infinite(fh) || tri.is_infinite(fh->neighbor(ei))) {
        std::cerr << "[flip] WARNING: edge (" << ia << "," << ib
                  << ") is on boundary / infinite, not flipping\n";
        return;
    }

    // Perform the edge flip
    tri.flip(fh, ei);
}

} // namespace df
#include "flip.h"




namespace df {

// we now have to perform the edge flip on the current triangulation with given indices
// after the flip the vertex handles remain the same, we do not have to update the mapping from global index to vertex handles 
// since the flip only changes the connectivity of the triangulation, not the vertices themselves


    void apply_edge_flip(df::vertex_id ia,
                    df::vertex_id ib,
                    df::InputData& D)
    {
        // get vertex handles from global indices
        auto va = D.index_to_vertex_handle_current.at(ia);
        auto vb = D.index_to_vertex_handle_current.at(ib);

        // perform the edge flip in the current triangulation
        df::Tri2::Face_handle fh; int i = -1;
        // we need the face handle and index of the edge to flip
        D.tri_current.is_edge(va, vb, fh, i); 

        // flip the edge
        D.tri_current.flip(fh, i);

        // do we have to update the index_to_vertex_handle_current map?
        // no, since vertex handles remain the same after the flip
    }
} // namespace df
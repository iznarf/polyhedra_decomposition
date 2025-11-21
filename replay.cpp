#include "replay.h"
#include "visualization.h"
#include <iostream>



namespace df {

// === edge flip (a,b,c,d) ===
// a,b : edge to be flipped 
// c,d : opposite vertices, only used for optional sanity checks/debug
void edge_flip_replay(vertex_id ia,
                      vertex_id ib,
                      vertex_id ic,
                      vertex_id id,
                      InputData& D)
{
    Tri2& T = D.tri_replay;

    Tri2::Face_handle fh;
    int ei = -1;
    bool found = false;

    // find the edge (ia, ib) in the current triangulation by vertex ids
    for (auto e = T.finite_edges_begin(); e != T.finite_edges_end(); ++e) {
        auto f  = e->first;
        int  i  = e->second;

        auto va = f->vertex(T.cw(i));
        auto vb = f->vertex(T.ccw(i));

        vertex_id ja = va->info();
        vertex_id jb = vb->info();

        if ((ja == ia && jb == ib) || (ja == ib && jb == ia)) {
            fh    = f;
            ei    = i;
            found = true;
            break;
        }
    }

    if (!found) {
        std::cerr << "[replay] edge flip: edge (" << ia << "," << ib
                  << ") not found in tri_current\n";
        return;
    }

    if (T.is_infinite(fh) || T.is_infinite(fh->neighbor(ei))) {
        std::cerr << "[replay] edge flip: edge (" << ia << "," << ib
                  << ") is boundary / infinite, skipping flip\n";
        return;
    }

    // optional sanity: check that opposite vertices match (c,d)
    auto vc = fh->vertex(ei);
    auto gh = fh->neighbor(ei);
    int  j  = T.mirror_index(fh, ei);
    auto vd = gh->vertex(j);

    vertex_id jc = vc->info();
    vertex_id jd = vd->info();

    // just a warning if mismatch; not fatal
    if (!((jc == ic && jd == id) || (jc == id && jd == ic))) {
        std::cerr << "[replay] edge flip: recorded opposite vertices ("
                  << ic << "," << id << ") but current quad is ("
                  << jc << "," << jd << ")\n";
    }

    T.flip(fh, ei);
}

// === vertex insertion (face a,b,c, new vertex d) ===
void vertex_insertion_replay(vertex_id ia, vertex_id ib, vertex_id ic, vertex_id id, InputData& D) {
    Tri2& T   = D.tri_replay;
    const P2& insertion_point = D.points2d[id];

    // Option 1: use locate (robust, uses actual geometry)
    Tri2::Locate_type lt;
    int li;
    Tri2::Face_handle fh = T.locate(insertion_point, lt, li);


    // Optional: sanity check that this face is indeed (a,b,c)
    vertex_id fa = fh->vertex(0)->info();
    vertex_id fb = fh->vertex(1)->info();
    vertex_id fc = fh->vertex(2)->info();

    // sort and compare as sets
    std::array<vertex_id,3> rec  = { ia, ib, ic };
    std::array<vertex_id,3> face = { fa, fb, fc };
    std::sort(rec.begin(),  rec.end());
    std::sort(face.begin(), face.end());

    if (rec != face) {
        std::cerr << "[replay] vertex insertion: recorded face ("
                  << ia << "," << ib << "," << ic
                  << "), but locate() found face ("
                  << fa << "," << fb << "," << fc << ")\n";
        // not fatal; we still insert in the geometric face
    }

    Tri2::Vertex_handle vh = T.insert_in_face(insertion_point, fh);
    vh->info() = id;
}




void apply_step(const StepRecord& step, InputData& D) {
    if (step.kind == StepKind::EdgeFlip) {

        // apply the flip on tri_replay
        edge_flip_replay(step.a, step.b, step.c, step.d, D);

    } else if (step.kind == StepKind::VertexInsertion) {

        // insert vertex d into face (a,b,c) on tri_replay
        vertex_insertion_replay(step.a, step.b, step.c, step.d, D);
    }

    // after modifying tri_replay, just update the replay meshes
    viz::show_or_update_replay(D);
}

}



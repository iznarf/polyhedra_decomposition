#include "replay.h"
#include "visualization.h"
#include "input.h"
#include <imgui.h> 
#include <iostream>

namespace {
    // pointer to the data we are replaying
    df::InputData* g_replay_data = nullptr;

    // triangulation at step 0 (no flips applied, this is just tri_upper)
    df::Tri2 g_replay_base;

    // current step index: number of steps applied on tri_replay (0..N)
    int g_step_idx = 0;

    // reset tri_replay to base and apply first k steps
    void rebuild_to_step(int k) {
        if (!g_replay_data) return;

        int n = static_cast<int>(g_replay_data->step_history.size());
        if (k < 0) k = 0;
        if (k > n) k = n;

        g_step_idx = k;

        // reset to base triangulation
        g_replay_data->tri_replay = g_replay_base;

        // apply first k steps
        for (int i = 0; i < k; ++i) {
            df::apply_step(g_replay_data->step_history[i], *g_replay_data);
        }

        // if k == 0, no apply_step was called, so we still need to show the base mesh
        if (k == 0) {
            viz::show_or_update_replay(*g_replay_data);
        }
    }

} // namespace 






namespace df {

// edge flip (a,b,c,d)
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

// vertex insertion (face a,b,c, new vertex d) 
void vertex_insertion_replay(vertex_id ia, vertex_id ib, vertex_id ic, vertex_id id, InputData& D) {
    Tri2& T   = D.tri_replay;
    const P2& insertion_point = D.points2d[id];

    // locate the face (a,b,c) in current triangulation
    Tri2::Locate_type lt;
    int li;
    Tri2::Face_handle fh = T.locate(insertion_point, lt, li);

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

void init_replay(InputData& D) {
    g_replay_data = &D;

    g_step_idx = 0;

    // tri_replay currently contains the step-0 triangulation
    g_replay_base = D.tri_replay;

    g_step_idx = 0;
    viz::show_or_update_replay(D);
}



void replay_ui() {
    if (!g_replay_data) {
        ImGui::Text("replay not initialized.");
        return;
    }

    int n = (int) g_replay_data->step_history.size();

    if (n == 0) {
        ImGui::Text("no steps recorded");
    } else if (g_step_idx == 0) {
        ImGui::Text("initial state (no flips yet)");
    } else {
        // last applied step is index g_step_idx - 1
        const StepRecord& step = g_replay_data->step_history[g_step_idx - 1];

        const char* kindStr =
            (step.kind == StepKind::EdgeFlip) ? "2-2 flip" : "1-3 flip";

        ImGui::Text("replay (green) step: %d / %d", g_step_idx, n);
        ImGui::Text("kind: %s", kindStr);
        ImGui::Text("vertices: (%zu, %zu, %zu, %zu)",
                    step.a, step.b, step.c, step.d);
    }

    ImGui::Separator();

    if (ImGui::Button("prev") && g_step_idx > 0) {
        rebuild_to_step(g_step_idx - 1);
    }
    ImGui::SameLine();
    if (ImGui::Button("next") && g_step_idx < n) {   // allow going up to n
        rebuild_to_step(g_step_idx + 1);
    }
}

}



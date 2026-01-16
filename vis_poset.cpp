#include "vis_poset.h"
#include "poset.h"
#include "visualization.h"      
#include "geometry_utils.h"   
#include "poset.h"
#include "poset2.h"
#include "moebius.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>  

#include <glm/glm.hpp>
#include <imgui.h>

#include <vector>
#include <array>
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <polyscope/view.h>


using glm::vec3;
using df::vertex_id;

namespace viz_poset {
    void rebuild_downflip_network_filtered();
    void register_poset2_cover_edges(const std::vector<std::vector<int>>& cover_out);

}

namespace {

    // faces as local index triples for planar triangulation
    static std::vector<std::array<int,3>> faces_from_triangles(const df::Tri2& t, const std::unordered_map<vertex_id,int>& to_local) {
        std::vector<std::array<int,3>> F;
        F.reserve(t.number_of_faces());
        for (auto f = t.finite_faces_begin(); f != t.finite_faces_end(); ++f) {
            int a = to_local.at(f->vertex(0)->info());
            int b = to_local.at(f->vertex(1)->info());
            int c = to_local.at(f->vertex(2)->info());
            F.push_back({a,b,c});
        }
        return F;
    }

    // attach global id scalar quantity to polyscope mesh
    static void add_global_id_quantity(polyscope::SurfaceMesh* mesh, const std::vector<vertex_id>& ids) {
        std::vector<double> values;
        values.reserve(ids.size());
        for (auto id : ids) {
            values.push_back(static_cast<double>(id));
        }
        mesh->addVertexScalarQuantity("global id", values);
    }


    // function for computing levels for global poset visualization
    // levels:
    //  (A) reachable from root (0): longest DOWN-path length from root
    //  (B) unreachable: place using children levels if possible,
    //      otherwise start new "block" below everything
    static std::vector<int> compute_levels_longest_root_with_unreachable(const std::vector<pst::Node>& nodes) {
        const int n = (int)nodes.size();
        std::vector<int> level(n, -1);
        if (n == 0) return level;

        auto is_down = [](df::StepKind k) {
            return k == df::StepKind::EdgeFlip_down
                || k == df::StepKind::VertexInsertion_down
                || k == df::StepKind::VertexDeletion_down;
        };

        // build down flip adjacency + indegrees for topo sort
        std::vector<std::vector<int>> adj(n);
        std::vector<int> indeg(n, 0);

        for (int u = 0; u < n; ++u) {
            const auto& kids  = nodes[u].children;
            const auto& steps = nodes[u].child_steps;
            const std::size_t m = std::min(kids.size(), steps.size());

            for (std::size_t e = 0; e < m; ++e) {
                if (!is_down(steps[e].kind)) continue;
                int v = kids[e];
                if (v < 0 || v >= n) continue;
                adj[u].push_back(v);
                indeg[v] += 1;
            }
        }

        // kahn topo order on down flip subgraph
        std::queue<int> q;
        std::vector<int> indeg_work = indeg;
        for (int i = 0; i < n; ++i) if (indeg_work[i] == 0) q.push(i);

        std::vector<int> topo;
        topo.reserve(n);
        while (!q.empty()) {
            int u = q.front(); q.pop();
            topo.push_back(u);
            for (int v : adj[u]) {
                if (--indeg_work[v] == 0) q.push(v);
            }
        }

        if ((int)topo.size() < n) {
            std::cerr << "[vis_poset] WARNING: down-edge subgraph is not a DAG "
                    << "(cycle detected); level layout may be invalid.\n";
            // fallback: keep everything at 0
            std::fill(level.begin(), level.end(), 0);
            return level;
        }

        // (A) longest path from root over topo
        const int NEG_INF = std::numeric_limits<int>::min();
        std::vector<int> dist(n, NEG_INF);
        dist[0] = 0;

        for (int u : topo) {
            if (dist[u] == NEG_INF) continue;
            for (int v : adj[u]) {
                dist[v] = std::max(dist[v], dist[u] + 1);
            }
        }

        int max_level = 0;
        for (int i = 0; i < n; ++i) {
            if (dist[i] != NEG_INF) {
                level[i] = dist[i];
                max_level = std::max(max_level, level[i]);
            }
        }

        // (B) place unreachable nodes in reverse topo (children processed before parents)
        // rule: if any child has a level, put u at (min_child_level - 1)
        // if u has no children with levels, start new block below
        for (int ti = (int)topo.size() - 1; ti >= 0; --ti) {
            int u = topo[ti];
            if (level[u] != -1) continue; // already placed (reachable)

            int min_child = std::numeric_limits<int>::max();
            for (int v : adj[u]) {
                if (level[v] != -1) {
                    min_child = std::min(min_child, level[v]);
                }
            }

            if (min_child != std::numeric_limits<int>::max()) {
                level[u] = min_child - 1;
            } else {
                // new unreachable 
                level[u] = ++max_level;
            }

            max_level = std::max(max_level, level[u]);
        }

        // safety: clamp anything still -1 (should not happen!!) to 0
        for (int& lv : level) if (lv < 0) lv = 0;

        return level;
    }

    // make planar vertex positions for a given triangulation centered at (cx,cz)
    // in the X-Z plane, with a small uniform scale
    static std::vector<vec3> make_planar_poset_vertices(const std::vector<vertex_id>& ids, const std::vector<df::P2>& points2d, float cx, float cz, float scale){
        double minx =  std::numeric_limits<double>::infinity();
        double maxx = -std::numeric_limits<double>::infinity();
        double miny =  std::numeric_limits<double>::infinity();
        double maxy = -std::numeric_limits<double>::infinity();

        for (auto id : ids) {
            const auto& p = points2d[id];
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            minx = std::min(minx, x);
            maxx = std::max(maxx, x);
            miny = std::min(miny, y);
            maxy = std::max(maxy, y);
        }

        double cx_local = 0.5 * (minx + maxx);
        double cy_local = 0.5 * (miny + maxy);

        // avoid degenerate scales
        if (!(maxx > minx)) maxx = minx + 1.0;
        if (!(maxy > miny)) maxy = miny + 1.0;

        std::vector<vec3> V;
        V.reserve(ids.size());

        for (auto id : ids) {
            const auto& p = points2d[id];
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());

            double dx = (x - cx_local);
            double dy = (y - cy_local);

            float X = cx + static_cast<float>(dx * scale);
            float Z = cz + static_cast<float>(dy * scale);

            // y=0 plane in Polyscope, (x,z) in our 2D sense
            V.emplace_back(X, 0.0f, Z);
        }

        return V;
    }

    // make lifted vertices above same grid center (cx,cz)
    static std::vector<vec3> make_lifted_poset_vertices(const std::vector<vertex_id>& ids, const std::vector<df::P2>& points2d, float cx, float cz, float scale_xy, float scale_z) {
        double minx =  std::numeric_limits<double>::infinity();
        double maxx = -std::numeric_limits<double>::infinity();
        double miny =  std::numeric_limits<double>::infinity();
        double maxy = -std::numeric_limits<double>::infinity();

        for (auto id : ids) {
            const auto& p = points2d[id];
            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            minx = std::min(minx, x);
            maxx = std::max(maxx, x);
            miny = std::min(miny, y);
            maxy = std::max(maxy, y);
        }

        double cx_local = 0.5 * (minx + maxx);
        double cy_local = 0.5 * (miny + maxy);

        std::vector<vec3> V;
        V.reserve(ids.size());

        for (auto id : ids) {
            const auto& p = points2d[id];
            df::P3 lp = df::lift(p);

            double x = CGAL::to_double(p.x());
            double y = CGAL::to_double(p.y());
            double z = CGAL::to_double(lp.z());

            double dx = (x - cx_local);
            double dy = (y - cy_local);

            float X = cx + static_cast<float>(dx * scale_xy);
            float Z = cz + static_cast<float>(dy * scale_xy);
            float Y = static_cast<float>(z * scale_z); // height

            // Polyscope convention (x,y,z)
            V.emplace_back(X, Y, Z);
        }

        return V;
    }

    // ------------------------------ 
    // global poset visualization data
    //-------------------------------

    // poset nodes meshes

    // global storage of registered poset meshes for UI toggling
    std::vector<polyscope::SurfaceMesh*> g_poset_meshes_2d;
    std::vector<polyscope::SurfaceMesh*> g_poset_meshes_3d;

    // down_flip edges visualization <=_1 relation
    std::vector<glm::vec3> g_node_centers; // center position for each node
    polyscope::CurveNetwork* g_downflip_network = nullptr;

    bool g_show_poset_2d = true;
    bool g_show_poset_3d = true;
    bool g_show_down_flips = true;

    // poset2 cover relations visualization
    polyscope::CurveNetwork* g_poset2_cover_network = nullptr;
    bool g_show_poset2_covers = true;

    // local poset UI data
    // keep input + nodes so the UI can rebuild local neighborhoods
    const df::InputData*       g_input_for_poset = nullptr;
    std::vector<pst::Node>     g_nodes_for_ui;
    int                        g_center_node_idx = 0;
    int                        g_max_local_depth = 5;

    // poset2 interval data 
    // store poset2 so UI can query intervals
    static pst2::Poset2 g_poset2;

    // interval UI state
    static int  g_interval_x = 0;
    static int  g_interval_y = 0;
    static bool g_interval_active = false;
    static std::vector<char> g_interval_mask; // size = #nodes, 1 => visible
    static char g_interval_x_buf[32] = "0";
    static char g_interval_y_buf[32] = "0";

    // poset2 meet/join 
    static char g_pair_a_buf[32] = "0";
    static char g_pair_b_buf[32] = "0";

    // focus mask for meet/join 
    static bool g_focus_active = false;
    static std::vector<char> g_focus_mask;

    // mesh comparison global variables
    // comparison view
    static char g_cmp_a_buf[32] = "0";
    static char g_cmp_b_buf[32] = "0";
    
    static bool g_cmp_show_A_2d = false;
    static bool g_cmp_show_A_3d = true;
    static bool g_cmp_show_B_2d = false;
    static bool g_cmp_show_B_3d = true;

    static bool g_compare_mode = false;


    static polyscope::SurfaceMesh* g_cmpA_2d = nullptr;
    static polyscope::SurfaceMesh* g_cmpB_2d = nullptr;
    static polyscope::SurfaceMesh* g_cmpA_3d = nullptr;
    static polyscope::SurfaceMesh* g_cmpB_3d = nullptr;

    static void clear_compare_meshes() {
        if (g_cmpA_2d) { g_cmpA_2d->remove(); g_cmpA_2d = nullptr; }
        if (g_cmpB_2d) { g_cmpB_2d->remove(); g_cmpB_2d = nullptr; }
        if (g_cmpA_3d) { g_cmpA_3d->remove(); g_cmpA_3d = nullptr; }
        if (g_cmpB_3d) { g_cmpB_3d->remove(); g_cmpB_3d = nullptr; }
    }

    static void set_compare_enabled() {
        if (g_cmpA_2d) g_cmpA_2d->setEnabled(g_cmp_show_A_2d);
        if (g_cmpA_3d) g_cmpA_3d->setEnabled(g_cmp_show_A_3d);
        if (g_cmpB_2d) g_cmpB_2d->setEnabled(g_cmp_show_B_2d);
        if (g_cmpB_3d) g_cmpB_3d->setEnabled(g_cmp_show_B_3d);
    }

    static void set_poset_enabled(bool enabled) {
        // node meshes
        for (auto* m : g_poset_meshes_2d) if (m) m->setEnabled(enabled && g_show_poset_2d);
        for (auto* m : g_poset_meshes_3d) if (m) m->setEnabled(enabled && g_show_poset_3d);

        // edge networks
        if (g_downflip_network) g_downflip_network->setEnabled(enabled && g_show_down_flips);
        if (g_poset2_cover_network) g_poset2_cover_network->setEnabled(enabled && g_show_poset2_covers);
    }

    // forward declarations (needed because helpers call these before their definitions)
    void rebuild_downflip_network_filtered();
    static bool node_visible(int idx);


    // restore visibility respecting interval/focus masks (same logic you already use elsewhere)
    static void apply_poset_visibility_from_masks() {
        int n = (int)g_nodes_for_ui.size();
        for (int i = 0; i < n; ++i) {
            bool on = node_visible(i);
            if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setEnabled(on && g_show_poset_2d);
            if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setEnabled(on && g_show_poset_3d);
        }

        viz_poset::rebuild_downflip_network_filtered();
        viz_poset::register_poset2_cover_edges(g_poset2.cover_out);
    }


    // helpers for meet/join UI coloring 
    static void reset_all_node_colors() {
        for (int i = 0; i < (int)g_poset_meshes_2d.size(); ++i) {
            if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setSurfaceColor(glm::vec3(0.6f, 0.8f, 1.0f));
            if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setSurfaceColor(glm::vec3(0.2f, 0.4f, 0.8f));
        }
    }

    static void color_node(int idx, glm::vec3 c2, glm::vec3 c3) {
        if (idx < 0 || idx >= (int)g_poset_meshes_2d.size()) return;
        if (g_poset_meshes_2d[idx]) g_poset_meshes_2d[idx]->setSurfaceColor(c2);
        if (g_poset_meshes_3d[idx]) g_poset_meshes_3d[idx]->setSurfaceColor(c3);
    }

    static std::vector<int> find_path_up_cover_out(const std::vector<std::vector<int>>& cover_out,int start, int goal, const std::vector<char>* allowed_mask = nullptr) {
        const int n = (int)cover_out.size();
        if (start < 0 || start >= n || goal < 0 || goal >= n) return {};
        if (start == goal) return {start};

        std::vector<int> parent(n, -1);
        std::queue<int> q;

        auto ok = [&](int v) {
            return !allowed_mask || (v >= 0 && v < (int)allowed_mask->size() && (*allowed_mask)[v]);
        };

        if (!ok(start) || !ok(goal)) return {};

        parent[start] = start;
        q.push(start);

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : cover_out[u]) { // go UP
                if (v < 0 || v >= n) continue;
                if (!ok(v)) continue;
                if (parent[v] != -1) continue;
                parent[v] = u;
                if (v == goal) break;
                q.push(v);
            }
        }

        if (parent[goal] == -1) return {}; // no path found

        // reconstruct
        std::vector<int> path;
        for (int cur = goal; cur != start; cur = parent[cur]) {
            path.push_back(cur);
        }
        path.push_back(start);
        std::reverse(path.begin(), path.end());
        return path;
    }

    // node visibility based on focus / interval masks 
    static bool node_visible(int idx) {
        if (g_focus_active) {
            return idx >= 0 && idx < (int)g_focus_mask.size() && g_focus_mask[idx];
        }
        if (g_interval_active) {
            return idx >= 0 && idx < (int)g_interval_mask.size() && g_interval_mask[idx];
        }
        return true;
    }


} // anonymous namespace


namespace viz_poset {

// this function constructs and registers all poset1 node meshes for global visualization
// returns number of down edges in relation <=1 registered

int register_poset(const df::InputData& D, const std::vector<pst::Node>& nodes) {
    int down_edge_count = 0;
    // clear previous poset meshes if any
    for (auto* m : g_poset_meshes_2d) {
        if (m) m->remove();
    }
    for (auto* m : g_poset_meshes_3d) {
        if (m) m->remove();
    }
    g_poset_meshes_2d.clear();
    g_poset_meshes_3d.clear();

    // clear previous down-flip network if any
    if (g_downflip_network) {
        g_downflip_network->remove();
        g_downflip_network = nullptr;
    }

    if (g_poset2_cover_network) {
        g_poset2_cover_network->remove();
        g_poset2_cover_network = nullptr;
    }


    // remember input + nodes for the local-poset UI
    g_input_for_poset = &D;
    g_nodes_for_ui    = nodes;
    g_center_node_idx = 0;


    if (nodes.empty()) {
        std::cout << "[vis_poset] no nodes to visualize.\n";
        return down_edge_count;
    }

    const int n = static_cast<int>(nodes.size());

    // positions for each node center (for curve network)
    g_node_centers.assign(n, glm::vec3(0.0f));

    // compute level for each node
    auto level = compute_levels_longest_root_with_unreachable(nodes);


    // keep this for the whole poset version -> visualize nodes that are not reachable from root
    int max_level = 0;
    for (int lv : level) {
        max_level = std::max(max_level, lv);
    }

    // build list of nodes per level
    std::vector<std::vector<int>> byLevel(max_level + 1);
    for (int i = 0; i < n; ++i) {
        int lv = level[i];
        if (lv < 0) lv = 0;
        byLevel[lv].push_back(i);
    }


    // grid layout parameters
    const float LEVEL_SPACING   = 4.0f;  // vertical spacing between levels (in Z)
    //const float LEVEL_SPACING   = 0.0f;  // vertical spacing between levels (in Z)
    const float NODE_SPACING    = 4.0f;  // horizontal spacing between nodes (in X)
    //const float NODE_SPACING    = 0.0f;  // horizontal spacing between nodes (in X)
    const float TRI_SCALE_PLAN  = 0.7f;  // scale for 2D shape
    const float TRI_SCALE_LIFT  = 0.7f;  // horizontal scale for 3D
    const float LIFT_HEIGHT_SCL = 0.5f;  // vertical scale for lift

    g_poset_meshes_2d.assign(n, nullptr);
    g_poset_meshes_3d.assign(n, nullptr);


    // for each node, reconstruct triangulation and register meshes
    for (int lv = 0; lv <= max_level; ++lv) {
        auto& nodesAtLevel = byLevel[lv];
        int k = static_cast<int>(nodesAtLevel.size());
        if (k == 0) continue;

        for (int idx_in_level = 0; idx_in_level < k; ++idx_in_level) {
            int node_idx = nodesAtLevel[idx_in_level];
            const auto& node = nodes[node_idx];

            float row_z  = static_cast<float>(lv) * LEVEL_SPACING;
            float base_x = (idx_in_level - 0.5f * (k - 1)) * NODE_SPACING;

            // store a simple center position for this node
            g_node_centers[node_idx] = glm::vec3(base_x, 0.0f, row_z);

            // reconstruct triangulation at this node
            df::Tri2 tri = D.tri_poset;                    // root is upper triangulation
            pst::replay_history_poset(tri, node.history, D);  // apply steps

            // collect ids and build faces
            auto ids      = viz::present_ids(tri);
            auto to_local = viz::make_local_index(ids);
            auto faces    = faces_from_triangles(tri, to_local);

            // planar vertices on grid
            auto V2 = make_planar_poset_vertices(ids, D.points2d,
                                                 base_x, row_z,
                                                 TRI_SCALE_PLAN);
            // lifted vertices above same grid cell
            auto V3 = make_lifted_poset_vertices(ids, D.points2d,
                                                 base_x, row_z,
                                                 TRI_SCALE_LIFT,
                                                 LIFT_HEIGHT_SCL);

            // unique names per node
            std::string name2d = "poset node " + std::to_string(node_idx) + " 2D";
            std::string name3d = "poset node " + std::to_string(node_idx) + " lifted";

            auto* m2 = polyscope::registerSurfaceMesh(name2d, V2, faces);
            add_global_id_quantity(m2, ids);
            m2->setEnabled(g_show_poset_2d);
            glm::vec3 col2(0.6f, 0.8f, 1.0f);
            m2->setSurfaceColor(col2);
            m2->setEdgeWidth(1.0f);                // non-zero => edges visible
            m2->setEdgeColor(glm::vec3(0, 0, 0));  // black edges

            auto* m3 = polyscope::registerSurfaceMesh(name3d, V3, faces);
            add_global_id_quantity(m3, ids);
            m3->setEnabled(g_show_poset_3d);
            glm::vec3 col3(0.2f, 0.4f, 0.8f);
            m3->setSurfaceColor(col3);
            m3->setTransparency(0.6f);
            m3->setEdgeWidth(1.0f);                // non-zero => edges visible
            m3->setEdgeColor(glm::vec3(0, 0, 0));  // black edges

            g_poset_meshes_2d[node_idx] = m2;
            g_poset_meshes_3d[node_idx] = m3;

        }
    }

       // build down-flip edges from parent/children and child_steps
        std::vector<glm::uvec2> down_edges;
        down_edges.reserve(nodes.size());

        for (int parent = 0; parent < n; ++parent) {
            const auto& kids  = nodes[parent].children;
            const auto& steps = nodes[parent].child_steps;

            // child_steps is parallel to children
            for (std::size_t e = 0; e < kids.size(); ++e) {
                int child = kids[e];
                const df::StepRecord& step = steps[e];

                bool is_down_step =
                    step.kind == df::StepKind::EdgeFlip_down ||
                    step.kind == df::StepKind::VertexInsertion_down ||
                    step.kind == df::StepKind::VertexDeletion_down;

                if (!is_down_step)
                    continue;

                down_edges.emplace_back(
                    static_cast<unsigned int>(parent),
                    static_cast<unsigned int>(child));
            }
        }

    if (!down_edges.empty()) {
        g_downflip_network = polyscope::registerCurveNetwork(
            "poset down flips", g_node_centers, down_edges);

        g_downflip_network->setEnabled(g_show_down_flips);
        g_downflip_network->setRadius(0.00037f, true);              // constant screen size

        g_downflip_network->setColor(glm::vec3(1.0f, 0.0f, 0.0f)); // 
    }

    std::cout << "[vis_poset] registered " << g_poset_meshes_2d.size()
              << " poset nodes (2D+3D meshes).\n";
    // print size of down_edges
    std::cout << "[vis_poset] registered " << down_edges.size()
              << " down-flip edges.\n";
    down_edge_count = static_cast<int>(down_edges.size());
    return down_edge_count;
}


static void rebuild_downflip_network_filtered() {
    const int n = (int)g_nodes_for_ui.size();
    if (n == 0) return;

    if (g_downflip_network) {
        g_downflip_network->remove();
        g_downflip_network = nullptr;
    }

    std::vector<glm::uvec2> down_edges;
    down_edges.reserve(n);

    for (int parent = 0; parent < n; ++parent) {

        if (!node_visible(parent)) continue;


        const auto& kids  = g_nodes_for_ui[parent].children;
        const auto& steps = g_nodes_for_ui[parent].child_steps;

        for (std::size_t e = 0; e < kids.size(); ++e) {
            int child = kids[e];
            const df::StepRecord& step = steps[e];

            bool is_down_step =
                step.kind == df::StepKind::EdgeFlip_down ||
                step.kind == df::StepKind::VertexInsertion_down ||
                step.kind == df::StepKind::VertexDeletion_down;

            if (!is_down_step) continue;
            if (child < 0 || child >= n) continue;

            if (!node_visible(child)) continue;


            down_edges.emplace_back((unsigned)parent, (unsigned)child);
        }
    }

    if (!down_edges.empty()) {
        g_downflip_network = polyscope::registerCurveNetwork("poset down flips", g_node_centers, down_edges);
        g_downflip_network->setEnabled(g_show_down_flips);
        g_downflip_network->setRadius(0.00037f, true);
        g_downflip_network->setColor(glm::vec3(1.0f, 0.0f, 0.0f));
    }
}


// ImGui UI for global poset visualization toggles
void poset_ui() {
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (ImGui::CollapsingHeader("poset1 visualization", ImGuiTreeNodeFlags_DefaultOpen)){
    
        if (g_poset_meshes_2d.empty() && g_poset_meshes_3d.empty()) {
            ImGui::Text("poset: no meshes registered");
            return;
        }

        if (ImGui::Checkbox("show all poset 2D meshes", &g_show_poset_2d)) {
            for (auto* m : g_poset_meshes_2d) {
                if (m) m->setEnabled(g_show_poset_2d);
            }
        }

        if (ImGui::Checkbox("show all poset 3D meshes", &g_show_poset_3d)) {
            for (auto* m : g_poset_meshes_3d) {
                if (m) m->setEnabled(g_show_poset_3d);
            }
        }
        if (ImGui::Checkbox("show down flips", &g_show_down_flips)) {
            if (g_downflip_network) {
                g_downflip_network->setEnabled(g_show_down_flips);
            }
        }
    }

    ImGui::Separator();
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (ImGui::CollapsingHeader("poset2 interval visualization [x,y]", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::InputText("x (lower)", g_interval_x_buf, IM_ARRAYSIZE(g_interval_x_buf));
        ImGui::InputText("y (upper)", g_interval_y_buf, IM_ARRAYSIZE(g_interval_y_buf));

        // parse input node numbers
        g_interval_x = std::atoi(g_interval_x_buf);
        g_interval_y = std::atoi(g_interval_y_buf);


        if (ImGui::Button("visualize [x,y]")) {

            // compute interval nodes
            auto nodes = pst2::interval_xy(g_poset2, g_interval_x, g_interval_y);

            //print nodes in interval
            std::cout << "[vis_poset] interval [" << g_interval_x << "," << g_interval_y << "] includes "
                    << nodes.size() << " nodes:\n";   
            for (int v : nodes) {
                std::cout << v << " ";
            }
            std::cout << "\n";

            // build mask
            int n = (int)g_nodes_for_ui.size();
            g_interval_mask.assign(n, 0);
            for (int v : nodes) {
                if (0 <= v && v < n) g_interval_mask[v] = 1;
            }
            g_interval_active = true;

            // enable only meshes in interval
            for (int i = 0; i < (int)g_poset_meshes_2d.size(); ++i) {
                if (g_poset_meshes_2d[i]) {
                    bool on = (i < n && g_interval_mask[i]) && g_show_poset_2d;
                    g_poset_meshes_2d[i]->setEnabled(on);
                }
            }
            for (int i = 0; i < (int)g_poset_meshes_3d.size(); ++i) {
                if (g_poset_meshes_3d[i]) {
                    bool on = (i < n && g_interval_mask[i]) && g_show_poset_3d;
                    g_poset_meshes_3d[i]->setEnabled(on);
                }
            }

            // rebuild edge networks restricted to interval
            rebuild_downflip_network_filtered();
            register_poset2_cover_edges(g_poset2.cover_out);
        }

        ImGui::SameLine();

        if (ImGui::Button("show full poset")) {
            g_interval_active = false;
            g_interval_mask.clear();

            // restore meshes
            for (auto* m : g_poset_meshes_2d) if (m) m->setEnabled(g_show_poset_2d);
            for (auto* m : g_poset_meshes_3d) if (m) m->setEnabled(g_show_poset_3d);

            // restore edge networks
            rebuild_downflip_network_filtered();             // will rebuild full when interval_active=false
            register_poset2_cover_edges(g_poset2.cover_out); // rebuild full
        }
    }

    ImGui::Separator();
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (ImGui::CollapsingHeader("poset2 meet/join visualization", ImGuiTreeNodeFlags_DefaultOpen)) {
  

        ImGui::InputText("a", g_pair_a_buf, IM_ARRAYSIZE(g_pair_a_buf));
        ImGui::InputText("b", g_pair_b_buf, IM_ARRAYSIZE(g_pair_b_buf));

        int a = std::atoi(g_pair_a_buf);
        int b = std::atoi(g_pair_b_buf);

        if (ImGui::Button("visualize meet")) {
            g_interval_active = false;
            g_interval_mask.clear();

            reset_all_node_colors();

            int a = std::atoi(g_pair_a_buf);
            int b = std::atoi(g_pair_b_buf);

            auto meets = pst2::meet_candidates_xy(g_poset2, a, b);

            int n = (int)g_nodes_for_ui.size();
            g_focus_mask.assign(n, 0);

            auto mark_path = [&](const std::vector<int>& path) {
                for (int v : path) if (0 <= v && v < n) g_focus_mask[v] = 1;
            };

            // always include inputs + candidates
            if (0 <= a && a < n) g_focus_mask[a] = 1;
            if (0 <= b && b < n) g_focus_mask[b] = 1;
            for (int m : meets) if (0 <= m && m < n) g_focus_mask[m] = 1;

            // paths: meet m is below a and b, so show m -> a and m -> b (UP edges)
            for (int m : meets) {
                mark_path(find_path_up_cover_out(g_poset2.cover_out, m, a));
                mark_path(find_path_up_cover_out(g_poset2.cover_out, m, b));
            }

            g_focus_active = true;

            // enable only focused meshes
            for (int i = 0; i < n; ++i) {
                bool on = g_focus_mask[i];
                if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setEnabled(on && g_show_poset_2d);
                if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setEnabled(on && g_show_poset_3d);
            }

            // rebuild edge networks restricted to focus
            rebuild_downflip_network_filtered();
            register_poset2_cover_edges(g_poset2.cover_out);

            // inputs yellow (same color)
            color_node(a, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));
            color_node(b, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));

            // meet candidates pink
            for (int m : meets) {
                color_node(m, glm::vec3(1.0f, 0.2f, 0.6f), glm::vec3(1.0f, 0.2f, 0.6f));
            }

            std::cout << "[vis_poset] number of meet candidates for (" << a << "," << b << "): "
                    << meets.size() << "\n";
            for (int m : meets) std::cout << m << " ";
            std::cout << "\n";
        }


        ImGui::SameLine();

        if (ImGui::Button("visualize join")) {
            g_interval_active = false;
            g_interval_mask.clear();

            reset_all_node_colors();

            int a = std::atoi(g_pair_a_buf);
            int b = std::atoi(g_pair_b_buf);

            auto joins = pst2::join_candidates_xy(g_poset2, a, b);

            // build focus mask
            int n = (int)g_nodes_for_ui.size();
            g_focus_mask.assign(n, 0);

            auto mark_path = [&](const std::vector<int>& path) {
                for (int v : path) if (0 <= v && v < n) g_focus_mask[v] = 1;
            };

            // always include inputs + candidates
            if (0 <= a && a < n) g_focus_mask[a] = 1;
            if (0 <= b && b < n) g_focus_mask[b] = 1;
            for (int j : joins) if (0 <= j && j < n) g_focus_mask[j] = 1;

            // add paths a->j and b->j
            for (int j : joins) {
                mark_path(find_path_up_cover_out(g_poset2.cover_out, a, j));
                mark_path(find_path_up_cover_out(g_poset2.cover_out, b, j));
            }

            g_focus_active = true;

            // enable only focused meshes
            for (int i = 0; i < n; ++i) {
                bool on = g_focus_mask[i];
                if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setEnabled(on && g_show_poset_2d);
                if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setEnabled(on && g_show_poset_3d);
            }

            // rebuild networks filtered by focus mask
            rebuild_downflip_network_filtered();
            register_poset2_cover_edges(g_poset2.cover_out);

            // color inputs yellow
            color_node(a, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));
            color_node(b, glm::vec3(1.0f, 1.0f, 0.1f), glm::vec3(1.0f, 1.0f, 0.1f));

            // color join candidates green
            for (int j : joins) {
                color_node(j, glm::vec3(0.2f, 1.0f, 0.4f), glm::vec3(0.2f, 1.0f, 0.4f));
            }

            std::cout << "[vis_poset] number of join candidates for (" << a << "," << b << "): "
                << joins.size() << "\n";
            for (int j : joins) std::cout << j << " ";
            std::cout << "\n";
        }


            ImGui::SameLine();

        
            if (ImGui::Button("reset meet/join view")) {
                g_focus_active = false;
                g_focus_mask.clear();

                // restore meshes
                for (auto* m : g_poset_meshes_2d) if (m) m->setEnabled(g_show_poset_2d);
                for (auto* m : g_poset_meshes_3d) if (m) m->setEnabled(g_show_poset_3d);

                rebuild_downflip_network_filtered();
                register_poset2_cover_edges(g_poset2.cover_out);
                reset_all_node_colors();
            }
    }
    ImGui::SetNextItemOpen(false, ImGuiCond_Once);
    if (ImGui::CollapsingHeader("compare nodes (overlay)", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::InputText("A", g_cmp_a_buf, IM_ARRAYSIZE(g_cmp_a_buf));
        ImGui::InputText("B", g_cmp_b_buf, IM_ARRAYSIZE(g_cmp_b_buf));

        bool changed = false;

        changed |= ImGui::Checkbox("A 2D", &g_cmp_show_A_2d);
        ImGui::SameLine();
        changed |= ImGui::Checkbox("A 3D", &g_cmp_show_A_3d);

        changed |= ImGui::Checkbox("B 2D", &g_cmp_show_B_2d);
        ImGui::SameLine();
        changed |= ImGui::Checkbox("B 3D", &g_cmp_show_B_3d);

        if (changed) set_compare_enabled();


        if (ImGui::Button("build compare overlay")) {

            clear_compare_meshes();

            if (!g_input_for_poset) {
                std::cerr << "[vis_poset] compare overlay: g_input_for_poset is null\n";
            } else {

            const int a = std::atoi(g_cmp_a_buf);
            const int b = std::atoi(g_cmp_b_buf);

            auto build_one = [&](int node_idx,polyscope::SurfaceMesh*& out2d, polyscope::SurfaceMesh*& out3d, const std::string& tag){
                if (node_idx < 0 || node_idx >= (int)g_nodes_for_ui.size()) {
                    std::cerr << "[vis_poset] compare overlay: invalid node " << node_idx << "\n";
                    return;
                }

                // reconstruct triangulation at node_idx
                df::Tri2 tri = g_input_for_poset->tri_poset;
                pst::replay_history_poset(tri, g_nodes_for_ui[node_idx].history, *g_input_for_poset);

                auto ids      = viz::present_ids(tri);
                auto to_local = viz::make_local_index(ids);
                auto faces    = faces_from_triangles(tri, to_local);

                float S = 20.0f; // try 10, 20, 50

                // (these helpers already re-center locally by bbox center)
                auto V2 = make_planar_poset_vertices(ids, g_input_for_poset->points2d,
                                                     0.0f, 20.0f, S);

                auto V3 = make_lifted_poset_vertices(ids, g_input_for_poset->points2d,
                                                     0.0f, 20.0f, S, 0.5f * S);

                std::string name2d = "compare " + tag + " node " + std::to_string(node_idx) + " 2D";
                std::string name3d = "compare " + tag + " node " + std::to_string(node_idx) + " lifted";

                out2d = polyscope::registerSurfaceMesh(name2d, V2, faces);
                out3d = polyscope::registerSurfaceMesh(name3d, V3, faces);

                add_global_id_quantity(out2d, ids);
                add_global_id_quantity(out3d, ids);

                // 2D: visible edges
                out2d->setEdgeWidth(1.0f);
                out2d->setEdgeColor(glm::vec3(0,0,0));

                // 3D: NOT transparent for compare
                out3d->setEdgeWidth(1.0f);
                out3d->setEdgeColor(glm::vec3(0,0,0));
            };

            build_one(a, g_cmpA_2d, g_cmpA_3d, "A");
            build_one(b, g_cmpB_2d, g_cmpB_3d, "B");

            // enforce initial visibility based on toggles
            set_compare_enabled();

            // hide the huge poset so camera/zoom acts on compare meshes only
            g_compare_mode = true;
            set_poset_enabled(false);


        }
    }

    ImGui::SameLine();
    if (ImGui::Button("clear compare overlay")) {
        clear_compare_meshes();

        // restore poset view
        g_compare_mode = false;
        apply_poset_visibility_from_masks(); // respects focus/interval if active
    }
}
}



// register / update a second edge network that visualizes <=2 cover edges
void register_poset2_cover_edges(const std::vector<std::vector<int>>& cover_out)
{
    const int n = (int)cover_out.size();
    if (n == 0) return;

    // store in global poset2 (so UI can compute intervals)
    g_poset2.cover_out = cover_out;

    // node centers
    if ((int)g_node_centers.size() != n) {
        std::cerr << "[vis_poset] register_poset2_cover_edges: node centers size mismatch.\n";
        return;
    }

    // remove old
    if (g_poset2_cover_network) {
        g_poset2_cover_network->remove();
        g_poset2_cover_network = nullptr;
    }

    // build edges
    std::vector<glm::uvec2> edges;
    edges.reserve(n);

    for (int u = 0; u < n; ++u) {

        if (!node_visible(u)) continue;

        for (int v : cover_out[u]) {
            if (v < 0 || v >= n) continue;

            if (!node_visible(v)) continue;

            // cover_out: u -> v means u <=2 v (bottom->top)
            // visualize top->down: reverse edge
            edges.emplace_back((unsigned)v, (unsigned)u);
        }
    }

    if (edges.empty()) return;

    g_poset2_cover_network = polyscope::registerCurveNetwork("poset2 covers", g_node_centers, edges);
    g_poset2_cover_network->setEnabled(g_show_poset2_covers);
    g_poset2_cover_network->setRadius(0.00037f, true);
    g_poset2_cover_network->setColor(glm::vec3(0.0f, 0.8f, 0.0f));
}

const std::vector<pst::Node>& get_poset1_nodes() {
    return g_nodes_for_ui;
}

const pst2::Poset2& get_poset2() {
    return g_poset2;
}


// color nodes by mobius function mu(0,T) in poset2
void color_nodes_by_mobius_0T() {
    const int n = (int)g_nodes_for_ui.size();
    if (n == 0) {
        std::cout << "[vis_poset] moebius coloring: no nodes.\n";
        return;
    }

    // make sure everything is visible (full poset), so you see all colors
    g_interval_active = false;
    g_interval_mask.clear();
    g_focus_active = false;
    g_focus_mask.clear();

    // enable node meshes based on current toggles
    for (int i = 0; i < n; ++i) {
        if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setEnabled(g_show_poset_2d);
        if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setEnabled(g_show_poset_3d);
    }

    // rebuild edge networks for full poset
    rebuild_downflip_network_filtered();
    register_poset2_cover_edges(g_poset2.cover_out);

    int cnt0 = 0, cntP = 0, cntN = 0, cntOther = 0;

    for (int t = 0; t < n; ++t) {
        std::int64_t mu = mob::mobius_xy_from_cover(g_poset2.cover_out, t, 0);


        glm::vec3 c;
        if (mu == 0) { c = glm::vec3(0.f, 0.f, 0.f); cnt0++; }          // black
        else if (mu == 1) { c = glm::vec3(1.f, 0.f, 0.f); cntP++; }     // red
        else if (mu == -1) { c = glm::vec3(0.f, 0.f, 1.f); cntN++; }    // blue
        else {
            // fallback: keep same sign coloring, but count it
            c = (mu > 0) ? glm::vec3(1.f, 0.f, 0.f) : glm::vec3(0.f, 0.f, 1.f);
            cntOther++;
        }

        // color both 2D and 3D mesh of the node
        if (t >= 0 && t < (int)g_poset_meshes_2d.size()) {
            if (g_poset_meshes_2d[t]) g_poset_meshes_2d[t]->setSurfaceColor(c);
            if (g_poset_meshes_3d[t]) g_poset_meshes_3d[t]->setSurfaceColor(c);
        }
    }

    std::cout << "[moebius] poset2 colored nodes by mu(0,T): "
              << "mu=0 black=" << cnt0 << ", mu=+1 red=" << cntP
              << ", mu=-1 blue=" << cntN;
    if (cntOther) std::cout << " (|mu|!=1: " << cntOther << ")";
    std::cout << "\n";
}

// color nodes by mobius function mu(0,T) in poset1
void color_nodes_by_mobius_poset1_0T() {
    const int n = (int)g_nodes_for_ui.size();
    if (n == 0) {
        std::cout << "[vis_poset] poset1 moebius coloring: no nodes.\n";
        return;
    }

    // Build upward cover graph for poset1: u->v means u <=1 v (bottom->top)
    auto cover1_up = mob::build_cover_up_from_poset1_nodes(g_nodes_for_ui);

    const int top = 0; // root is index 0

    // show full poset
    g_interval_active = false; g_interval_mask.clear();
    g_focus_active = false;    g_focus_mask.clear();

    for (int i = 0; i < n; ++i) {
        if (g_poset_meshes_2d[i]) g_poset_meshes_2d[i]->setEnabled(g_show_poset_2d);
        if (g_poset_meshes_3d[i]) g_poset_meshes_3d[i]->setEnabled(g_show_poset_3d);
    }

    // rebuild edges for full poset
    rebuild_downflip_network_filtered();

    int cnt0=0, cntP=0, cntN=0, cntOther=0;

    for (int t = 0; t < n; ++t) {
        // [0,T] convention: T is below 0, so compute mu(T,0)
        std::int64_t mu = mob::mobius_xy_from_cover(cover1_up, t, top);

        glm::vec3 c;
        if (mu == 0) { c = glm::vec3(0.f,0.f,0.f); cnt0++; }
        else if (mu == 1) { c = glm::vec3(1.f,0.f,0.f); cntP++; }
        else if (mu == -1) { c = glm::vec3(0.f,0.f,1.f); cntN++; }
        else { c = (mu > 0) ? glm::vec3(1.f,0.f,0.f) : glm::vec3(0.f,0.f,1.f); cntOther++; }

        if (g_poset_meshes_2d[t]) g_poset_meshes_2d[t]->setSurfaceColor(c);
        if (g_poset_meshes_3d[t]) g_poset_meshes_3d[t]->setSurfaceColor(c);
    }

    std::cout << "[moebius] poset1 colored nodes by mu(0,T): "
              << "mu=0 black=" << cnt0 << ", mu=+1 red=" << cntP << ", mu=-1 blue=" << cntN;
    if (cntOther) std::cout << " (|mu|!=1: " << cntOther << ")";
    std::cout << "\n";
}



// reset all node colors to default
// wrapper for moebius coloring
void reset_node_coloring() {
    reset_all_node_colors();
    std::cout << "[vis_poset] reset node coloring to default\n";
}











} // namespace viz_poset

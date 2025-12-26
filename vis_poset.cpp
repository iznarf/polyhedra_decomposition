#include "vis_poset.h"
#include "poset.h"
#include "visualization.h"      
#include "geometry_utils.h"   
#include "poset.h"
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

using glm::vec3;
using df::vertex_id;

namespace {

    // faces as local index triples, like in visualization.cpp 
    static std::vector<std::array<int,3>>
    faces_from_triangles(const df::Tri2& t,
                        const std::unordered_map<vertex_id,int>& to_local)
    {
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

    // attach "global id" scalar quantity to a Polyscope mesh
    static void add_global_id_quantity(polyscope::SurfaceMesh* mesh,
                                    const std::vector<vertex_id>& ids)
    {
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
//      otherwise start new "block" below everything.
static std::vector<int>
compute_levels_longest_root_with_unreachable(const std::vector<pst::Node>& nodes)
{
    const int n = (int)nodes.size();
    std::vector<int> level(n, -1);
    if (n == 0) return level;

    auto is_down = [](df::StepKind k) {
        return k == df::StepKind::EdgeFlip_down
            || k == df::StepKind::VertexInsertion_down
            || k == df::StepKind::VertexDeletion_down;
    };

    // Build DOWN adjacency + indegrees for topo sort
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

    // Kahn topo order on DOWN subgraph
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

    // (A) Longest-path-from-root DP over topo
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

    // (B) Place unreachable nodes in reverse topo (children processed before parents)
    // Rule: if any child has a level, put u at (min_child_level - 1).
    // If u has no children with levels (sink in its component), start a new block below.
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
            // new unreachable sink block
            level[u] = ++max_level;
        }

        max_level = std::max(max_level, level[u]);
    }

    // Safety: clamp anything still -1 (shouldn't happen) to 0
    for (int& lv : level) if (lv < 0) lv = 0;

    return level;
}





    // make planar vertex positions for a given triangulation centered at (cx,cz)
    // in the X-Z plane, with a small uniform scale.
    static std::vector<vec3>
    make_planar_poset_vertices(const std::vector<vertex_id>& ids,
                            const std::vector<df::P2>& points2d,
                            float cx, float cz,
                            float scale)
    {
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
    static std::vector<vec3>
    make_lifted_poset_vertices(const std::vector<vertex_id>& ids,
                            const std::vector<df::P2>& points2d,
                            float cx, float cz,
                            float scale_xy,
                            float scale_z)
    {
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

            // Note: Polyscope convention (x,y,z)
            V.emplace_back(X, Y, Z);
        }

        return V;
    }



    // global storage of registered poset meshes for UI toggling
    std::vector<polyscope::SurfaceMesh*> g_poset_meshes_2d;
    std::vector<polyscope::SurfaceMesh*> g_poset_meshes_3d;

    // down_flip edges visualization
    std::vector<glm::vec3> g_node_centers; // center position for each node
    polyscope::CurveNetwork* g_downflip_network = nullptr;

    bool g_show_poset_2d = true;
    bool g_show_poset_3d = true;
    bool g_show_down_flips = true;

    polyscope::CurveNetwork* g_poset2_cover_network = nullptr;
    bool g_show_poset2_covers = true;



    // --- global state so we can rebuild local posets from the UI ---

    // keep input + nodes so the UI can rebuild local neighborhoods
    const df::InputData*       g_input_for_poset = nullptr;
    std::vector<pst::Node>     g_nodes_for_ui;
    int                        g_center_node_idx = 0;
    int                        g_max_local_depth = 5;


} // anonymous namespace




namespace viz_poset {

void register_poset(const df::InputData& D, const std::vector<pst::Node>& nodes) {
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
        return;
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

    // For each node, reconstruct triangulation and register meshes
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

            g_poset_meshes_2d.push_back(m2);
            g_poset_meshes_3d.push_back(m3);
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
}

void poset_ui() {
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

    /*
    ImGui::Separator();
    ImGui::Text("local down flip poset from chosen root");

    if (!g_input_for_poset || g_nodes_for_ui.empty()) {
        ImGui::Text("no poset data available.");
        return;
    }

    int max_index = static_cast<int>(g_nodes_for_ui.size()) - 1;
    if (max_index < 0) return;

    if (g_center_node_idx > max_index) g_center_node_idx = max_index;
    if (g_center_node_idx < 0)        g_center_node_idx = 0;

    ImGui::SliderInt("center node index",
                     &g_center_node_idx,
                     0, max_index);

    ImGui::SliderInt("max local depth",
                     &g_max_local_depth,
                     1, 10);

    if (ImGui::Button("rebuild local poset")) {
        const auto& center_hist =
            g_nodes_for_ui[g_center_node_idx].history;

        std::vector<pst::Node> local_nodes;
        pst::build_poset_local_down_from_history(
            *g_input_for_poset,
            center_hist,
            static_cast<std::size_t>(g_max_local_depth),
            local_nodes,
            600
        );

        pst::debug_print_local_poset_histories(
        local_nodes,
        center_hist.size()  // center_history_len
        );

        // show this local poset instead
        register_poset(*g_input_for_poset, local_nodes);
        
    }
    */
}

void register_poset2_cover_edges(const std::vector<std::vector<int>>& cover_out)
{
    const int n = (int)cover_out.size();
    if (n == 0) return;

    // node centers
    if ((int)g_node_centers.size() != n) {
        std::cerr << "[vis_poset] register_poset2_cover_edges: "
                  << "node centers not initialized or size mismatch.\n";
        return;
    }

    // remove 
    if (g_poset2_cover_network) {
        g_poset2_cover_network->remove();
        g_poset2_cover_network = nullptr;
    }

    // build edges 
    // cover_out: u -> v means u <=2 v (bottom->top)
    // for visualization like poset1 (top->down), reverse: v -> u
    std::vector<glm::uvec2> edges;
    edges.reserve(n);

    for (int u = 0; u < n; ++u) {
        for (int v : cover_out[u]) {
            if (v < 0 || v >= n) continue;
            edges.emplace_back((unsigned)v, (unsigned)u); // reversed for display
        }
    }

    if (edges.empty()) return;

    g_poset2_cover_network = polyscope::registerCurveNetwork(
        "poset2 covers", g_node_centers, edges);

    g_poset2_cover_network->setEnabled(g_show_poset2_covers);
    g_poset2_cover_network->setRadius(0.00037f, true);           
    g_poset2_cover_network->setColor(glm::vec3(0.0f, 0.8f, 0.0f)); // green
}


} // namespace viz_poset

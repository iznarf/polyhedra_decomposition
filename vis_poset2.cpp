#include "vis_poset2.h"

#include "visualization.h"
#include "geometry_utils.h"
#include "replay.h"

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
#include <limits>

using glm::vec3;
using df::vertex_id;

namespace {

using UEdge = glm::uvec2;

// --- small helpers copied from vis_poset.cpp ---

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

static void add_global_id_quantity(polyscope::SurfaceMesh* mesh,
                                   const std::vector<vertex_id>& ids)
{
    std::vector<double> values;
    values.reserve(ids.size());
    for (auto id : ids) values.push_back((double)id);
    mesh->addVertexScalarQuantity("global id", values);
}

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
        minx = std::min(minx, x); maxx = std::max(maxx, x);
        miny = std::min(miny, y); maxy = std::max(maxy, y);
    }

    double cx_local = 0.5 * (minx + maxx);
    double cy_local = 0.5 * (miny + maxy);
    if (!(maxx > minx)) maxx = minx + 1.0;
    if (!(maxy > miny)) maxy = miny + 1.0;

    std::vector<vec3> V;
    V.reserve(ids.size());

    for (auto id : ids) {
        const auto& p = points2d[id];
        double x = CGAL::to_double(p.x());
        double y = CGAL::to_double(p.y());
        float X = cx + (float)((x - cx_local) * scale);
        float Z = cz + (float)((y - cy_local) * scale);
        V.emplace_back(X, 0.0f, Z);
    }
    return V;
}

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
        minx = std::min(minx, x); maxx = std::max(maxx, x);
        miny = std::min(miny, y); maxy = std::max(maxy, y);
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

        float X = cx + (float)((x - cx_local) * scale_xy);
        float Z = cz + (float)((y - cy_local) * scale_xy);
        float Y = (float)(z * scale_z);
        V.emplace_back(X, Y, Z);
    }
    return V;
}

// Build "down edges" for layout & drawing: we want top->down.
// But P2.cover_out stores bottom->top (u<=2 v).
// So we reverse: v -> u becomes a down edge.
static void build_down_edges_from_cover(const std::vector<std::vector<int>>& cover_out,
                                       std::vector<std::vector<int>>& down_adj,
                                       std::vector<UEdge>& down_edges)
{
    int N = (int)cover_out.size();
    down_adj.assign(N, {});
    down_edges.clear();
    down_edges.reserve(N);

    for (int u = 0; u < N; ++u) {
        for (int v : cover_out[u]) {
            if (v < 0 || v >= N) continue;
            // cover_out: u -> v (u <=2 v), so reverse for drawing/layout:
            // v -> u (v is above u)
            down_adj[v].push_back(u);
            down_edges.emplace_back((unsigned)v, (unsigned)u);
        }
    }

    for (auto& nbrs : down_adj) {
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }
}

// Levels = longest path length from TOP node 0 along down edges.
// (If node 0 is not the unique top, this still gives a consistent layout;
// disconnected components get placed after.)
static std::vector<int>
compute_levels_longest_from_root0(const std::vector<std::vector<int>>& down_adj)
{
    int N = (int)down_adj.size();
    std::vector<int> indeg(N, 0);
    for (int u = 0; u < N; ++u) for (int v : down_adj[u]) indeg[v]++;

    std::queue<int> q;
    auto indeg_work = indeg;
    for (int i = 0; i < N; ++i) if (indeg_work[i] == 0) q.push(i);

    std::vector<int> topo;
    topo.reserve(N);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        topo.push_back(u);
        for (int v : down_adj[u]) if (--indeg_work[v] == 0) q.push(v);
    }

    if ((int)topo.size() != N) {
        std::cerr << "[poset2_vis] WARNING: cover graph has a cycle; layout fallback.\n";
        return std::vector<int>(N, 0);
    }

    const int NEG = std::numeric_limits<int>::min();
    std::vector<int> dist(N, NEG);
    if (N > 0) dist[0] = 0;

    for (int u : topo) {
        if (dist[u] == NEG) continue;
        for (int v : down_adj[u]) dist[v] = std::max(dist[v], dist[u] + 1);
    }

    // Place unreachable components below everything (simple strategy)
    int maxlv = 0;
    std::vector<int> level(N, -1);
    for (int i = 0; i < N; ++i) {
        if (dist[i] != NEG) {
            level[i] = dist[i];
            maxlv = std::max(maxlv, level[i]);
        }
    }
    for (int i = 0; i < N; ++i) {
        if (level[i] == -1) level[i] = ++maxlv;
    }
    return level;
}

// --- global UI state (same pattern as vis_poset.cpp) ---
std::vector<polyscope::SurfaceMesh*> g_meshes2d;
std::vector<polyscope::SurfaceMesh*> g_meshes3d;

std::vector<glm::vec3> g_node_centers;
polyscope::CurveNetwork* g_cover_network = nullptr;

bool g_show_2d = true;
bool g_show_3d = true;
bool g_show_covers = true;

} // namespace


namespace viz_poset2 {

void register_poset2(const df::InputData& D,
                     const std::vector<pst::Node>& nodes,
                     const pst2::Poset2& P2)
{
    // remove old
    for (auto* m : g_meshes2d) if (m) m->remove();
    for (auto* m : g_meshes3d) if (m) m->remove();
    g_meshes2d.clear();
    g_meshes3d.clear();

    if (g_cover_network) { g_cover_network->remove(); g_cover_network = nullptr; }

    if (nodes.empty()) {
        std::cout << "[poset2_vis] no nodes.\n";
        return;
    }
    if ((int)P2.cover_out.size() != (int)nodes.size()) {
        std::cerr << "[poset2_vis] ERROR: cover_out size != nodes size\n";
        return;
    }

    const int N = (int)nodes.size();

    // Build top->down adjacency + edge list from cover_out
    std::vector<std::vector<int>> down_adj;
    std::vector<UEdge> down_edges;
    build_down_edges_from_cover(P2.cover_out, down_adj, down_edges);

    // Compute levels for grid layout
    auto level = compute_levels_longest_from_root0(down_adj);
    int max_level = 0;
    for (int lv : level) max_level = std::max(max_level, lv);

    std::vector<std::vector<int>> byLevel(max_level + 1);
    for (int i = 0; i < N; ++i) byLevel[level[i]].push_back(i);

    // Layout params (feel free to match vis_poset exactly)
    const float LEVEL_SPACING   = 4.0f;
    const float NODE_SPACING    = 4.0f;
    const float TRI_SCALE_PLAN  = 0.7f;
    const float TRI_SCALE_LIFT  = 0.7f;
    const float LIFT_HEIGHT_SCL = 0.5f;

    g_node_centers.assign(N, glm::vec3(0.0f));

    // Register meshes per node (same reconstruction as poset1)
    for (int lv = 0; lv <= max_level; ++lv) {
        auto& list = byLevel[lv];
        int k = (int)list.size();
        if (k == 0) continue;

        for (int j = 0; j < k; ++j) {
            int node_idx = list[j];
            float row_z  = (float)lv * LEVEL_SPACING;
            float base_x = (j - 0.5f * (k - 1)) * NODE_SPACING;

            g_node_centers[node_idx] = glm::vec3(base_x, 0.0f, row_z);

            df::Tri2 tri = D.tri_poset;
            pst::replay_history_poset(tri, nodes[node_idx].history, D);

            auto ids      = viz::present_ids(tri);
            auto to_local = viz::make_local_index(ids);
            auto faces    = faces_from_triangles(tri, to_local);

            auto V2 = make_planar_poset_vertices(ids, D.points2d, base_x, row_z, TRI_SCALE_PLAN);
            auto V3 = make_lifted_poset_vertices(ids, D.points2d, base_x, row_z, TRI_SCALE_LIFT, LIFT_HEIGHT_SCL);

            std::string name2d = "poset2 node " + std::to_string(node_idx) + " 2D";
            std::string name3d = "poset2 node " + std::to_string(node_idx) + " lifted";

            auto* m2 = polyscope::registerSurfaceMesh(name2d, V2, faces);
            add_global_id_quantity(m2, ids);
            m2->setEnabled(g_show_2d);
            m2->setEdgeWidth(1.0f);
            m2->setEdgeColor(glm::vec3(0,0,0));

            auto* m3 = polyscope::registerSurfaceMesh(name3d, V3, faces);
            add_global_id_quantity(m3, ids);
            m3->setEnabled(g_show_3d);
            m3->setTransparency(0.6f);
            m3->setEdgeWidth(1.0f);
            m3->setEdgeColor(glm::vec3(0,0,0));

            g_meshes2d.push_back(m2);
            g_meshes3d.push_back(m3);
        }
    }

    // Register cover edges as a curve network (top->down edges)
    if (!down_edges.empty()) {
        g_cover_network = polyscope::registerCurveNetwork("poset2 covers", g_node_centers, down_edges);
        g_cover_network->setEnabled(g_show_covers);
        g_cover_network->setRadius(0.00037f, true); // constant screen size
        g_cover_network->setColor(glm::vec3(0.0f, 0.8f, 0.0f)); // green
    }

    std::cout << "[poset2_vis] registered " << N << " nodes and "
              << down_edges.size() << " cover edges.\n";
}



void poset2_ui()
{
    if (g_meshes2d.empty() && g_meshes3d.empty()) {
        ImGui::Text("poset2: no meshes registered");
        return;
    }

    if (ImGui::Checkbox("show poset2 2D meshes", &g_show_2d)) {
        for (auto* m : g_meshes2d) if (m) m->setEnabled(g_show_2d);
    }
    if (ImGui::Checkbox("show poset2 lifted meshes", &g_show_3d)) {
        for (auto* m : g_meshes3d) if (m) m->setEnabled(g_show_3d);
    }
    if (ImGui::Checkbox("show poset2 cover edges", &g_show_covers)) {
        if (g_cover_network) g_cover_network->setEnabled(g_show_covers);
    }
}

} // namespace viz_poset2

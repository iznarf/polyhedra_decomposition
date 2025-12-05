#include "vis_poset.h"

#include "visualization.h"      // present_ids / make_local_index
#include "geometry_utils.h"   
#include "poset.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
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

    // Faces as local index triples, like in visualization.cpp (but reimplemented)
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

    // Attach "global id" scalar quantity to a Polyscope mesh
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

// Compute a simple BFS level (distance from root) for each node in the poset
static std::vector<int>
compute_levels(const std::vector<pst::Node>& nodes)
{
    const int n = static_cast<int>(nodes.size());
    std::vector<int> level(n, -1);
    std::queue<int> q;

    level[0] = 0;
    q.push(0);

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        int lu = level[u];

        for (int c : nodes[u].children) {
            if (c < 0 || c >= n) continue;
            if (level[c] == -1 || level[c] > lu + 1) {
                level[c] = lu + 1;
                q.push(c);
            }
        }
    }

    // Some nodes might not be reachable from root (shouldn't happen ideally),
    // but if they are, just set them to level 0.
    for (int i = 0; i < n; ++i) {
        if (level[i] < 0) level[i] = 0;
    }

    return level;
}

// Make planar vertex positions for a given triangulation centered at (cx,cz)
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

// Make lifted vertices above same grid center (cx,cz)
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

// Global storage of registered poset meshes for UI toggling
std::vector<polyscope::SurfaceMesh*> g_poset_meshes_2d;
std::vector<polyscope::SurfaceMesh*> g_poset_meshes_3d;

bool g_show_poset_2d = true;
bool g_show_poset_3d = true;

} // anonymous namespace


namespace viz_poset {

void register_poset(const df::InputData& D,
                    const std::vector<pst::Node>& nodes)
{
    // clear previous poset meshes if any
    for (auto* m : g_poset_meshes_2d) {
        if (m) m->remove();
    }
    for (auto* m : g_poset_meshes_3d) {
        if (m) m->remove();
    }
    g_poset_meshes_2d.clear();
    g_poset_meshes_3d.clear();

    if (nodes.empty()) {
        std::cout << "[vis_poset] no nodes to visualize.\n";
        return;
    }

    // Compute BFS level for each node (distance from root)
    auto level = compute_levels(nodes);
    int n = static_cast<int>(nodes.size());

    int max_level = 0;
    for (int lv : level) max_level = std::max(max_level, lv);

    // Build list of nodes per level
    std::vector<std::vector<int>> byLevel(max_level + 1);
    for (int i = 0; i < n; ++i) {
        int lv = level[i];
        if (lv < 0) lv = 0;
        byLevel[lv].push_back(i);
    }

    // Grid layout parameters
    const float LEVEL_SPACING   = 4.0f;  // vertical spacing between levels (in Z)
    const float NODE_SPACING    = 4.0f;  // horizontal spacing between nodes (in X)
    const float TRI_SCALE_PLAN  = 0.7f;  // scale for 2D shape
    const float TRI_SCALE_LIFT  = 0.7f;  // horizontal scale for 3D
    const float LIFT_HEIGHT_SCL = 0.5f;  // vertical scale for lift

    // For each node, reconstruct triangulation and register meshes
    for (int lv = 0; lv <= max_level; ++lv) {
        auto& nodesAtLevel = byLevel[lv];
        int k = static_cast<int>(nodesAtLevel.size());
        if (k == 0) continue;

        // x positions centered around 0
        for (int idx_in_level = 0; idx_in_level < k; ++idx_in_level) {
            int node_idx = nodesAtLevel[idx_in_level];
            const auto& node = nodes[node_idx];

            float row_z  = static_cast<float>(lv) * LEVEL_SPACING;
            float base_x = (idx_in_level - 0.5f * (k - 1)) * NODE_SPACING;



            // reconstruct triangulation at this node
            df::Tri2 tri = D.tri_poset;            // root is upper triangulation
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
            m2->setEdgeColor(glm::vec3(0, 0, 0));  // black edges, or whatever

            auto* m3 = polyscope::registerSurfaceMesh(name3d, V3, faces);
            add_global_id_quantity(m3, ids);
            m3->setEnabled(g_show_poset_3d);
            glm::vec3 col3(0.2f, 0.4f, 0.8f);
            m3->setSurfaceColor(col3);
            m3->setTransparency(0.6f);

            m3->setEdgeWidth(1.0f);                // non-zero => edges visible
            m3->setEdgeColor(glm::vec3(0, 0, 0));  // black edges, or whatever


            g_poset_meshes_2d.push_back(m2);
            g_poset_meshes_3d.push_back(m3);
        }
    }

    std::cout << "[vis_poset] registered " << g_poset_meshes_2d.size()
              << " poset nodes (2D+3D meshes).\n";
}


void poset_ui()
{
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
}

} // namespace viz_poset

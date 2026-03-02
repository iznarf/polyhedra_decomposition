#include "poset_ear_star_vis.h"
#include "visualization.h"  

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <glm/glm.hpp>

#include <cmath>
#include <sstream>
#include <unordered_map>
#include <set>
#include <algorithm>

namespace pst_es_viz {

using glm::vec3;
using df::vertex_id;


static void register_triangulation_mesh_with_offset(
    const df::Tri2& tri,
    const std::vector<df::P2>& points2d,
    const std::string& name,
    const glm::vec3& offset
) {
    auto ids = viz::present_ids(tri);
    std::unordered_map<vertex_id, int> to_local = viz::make_local_index(ids);

    // vertices (XZ plane) + offset
    std::vector<vec3> V;
    V.reserve(ids.size());
    for (auto id : ids) {
        const auto& p = points2d.at(id);
        V.push_back(offset + vec3((float)p.x(), 0.0f, (float)p.y()));
    }

    // faces (finite only), dedup + skip degenerates
    std::vector<std::array<int,3>> F;
    F.reserve(tri.number_of_faces());

    std::set<std::array<int,3>> seen;

    for (auto f = tri.finite_faces_begin(); f != tri.finite_faces_end(); ++f) {
        int a = to_local.at(f->vertex(0)->info());
        int b = to_local.at(f->vertex(1)->info());
        int c = to_local.at(f->vertex(2)->info());

        if (a == b || b == c || a == c) continue;

        std::array<int,3> key = {a,b,c};
        std::sort(key.begin(), key.end());
        if (!seen.insert(key).second) continue;

        F.push_back({a,b,c});
    }

    if (polyscope::hasSurfaceMesh(name)) polyscope::removeSurfaceMesh(name);
    auto* m = polyscope::registerSurfaceMesh(name, V, F);
    m->setEnabled(true);

    // optional scalar = global id
    std::vector<double> gid;
    gid.reserve(ids.size());
    for (auto id : ids) gid.push_back((double)id);
    m->addVertexScalarQuantity("global id", gid);

    m->setSurfaceColor(glm::vec3(0.6f, 0.8f, 1.0f));
    m->setSmoothShade(false);
    m->setEdgeWidth(1.0);
}

void register_poset_triangulations_grid(
    const df::Tri2& tri_start,
    const pst_es::FlipPoset& P,
    const std::vector<df::P2>& points2d,
    df::vertex_id star_id,
    const df::P2& star_point,
    const std::string& base_name,
    double spacing
) {
    const int N = (int)P.nodes.size();
    if (N == 0) return;

    // grid
    int cols = (int)std::ceil(std::sqrt((double)N));
    int rows = (int)std::ceil((double)N / (double)cols);

    for (int idx = 0; idx < N; ++idx) {
        // reconstruct triangulation at node idx
        df::Tri2 tri;
        bool ok = pst_es::reconstruct_triangulation(
            tri_start, P, idx, star_id, star_point, tri
        );
        if (!ok) continue;

        int r = idx / cols;
        int c = idx % cols;

        // center grid around origin
        double xoff = spacing * (c - 0.5 * (cols - 1));
        double zoff = spacing * (r - 0.5 * (rows - 1));

        glm::vec3 offset((float)xoff, 0.0f, (float)zoff);

        std::ostringstream oss;
        oss << base_name << "_" << idx;

        register_triangulation_mesh_with_offset(tri, points2d, oss.str(), offset);
    }

    polyscope::view::resetCameraToHomeView();
}

} // namespace pst_es_viz
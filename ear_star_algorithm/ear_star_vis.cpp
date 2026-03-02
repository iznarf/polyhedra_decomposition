#include "visualization.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <array>
#include <set>
#include <algorithm>

using glm::vec3;
using df::vertex_id;

void register_triangulation_as_mesh_ear_star(
    const df::Tri2& tri,
    const std::vector<df::P2>& points2d,
    const std::string& name
) {
    auto ids = viz::present_ids(tri);
    std::unordered_map<vertex_id, int> to_local = viz::make_local_index(ids);

    //vertices in XZ-plane (y = 0)
    std::vector<vec3> V;
    V.reserve(ids.size());
    for (auto id : ids) {
        const auto& p = points2d.at(id);
        V.push_back(vec3((float)p.x(), 0.0f, (float)p.y()));
    }

    //faces (finite only), dedup + skip degenerates
    std::vector<std::array<int, 3>> F;
    F.reserve(tri.number_of_faces());

    std::set<std::array<int, 3>> seen;

    for (auto f = tri.finite_faces_begin(); f != tri.finite_faces_end(); ++f) {
        int a = to_local.at(f->vertex(0)->info());
        int b = to_local.at(f->vertex(1)->info());
        int c = to_local.at(f->vertex(2)->info());

        if (a == b || b == c || a == c) continue;

        std::array<int, 3> key = {a, b, c};
        std::sort(key.begin(), key.end());
        if (!seen.insert(key).second) continue;

        F.push_back({a, b, c});
    }

    std::cout << "ids.size()=" << ids.size()
          << " V.size()=" << V.size()
          << " F.size()=" << F.size()
          << " finite_faces=" << tri.number_of_faces()
          << std::endl;

    // avoid name collisions
    if (polyscope::hasSurfaceMesh(name)) {
        polyscope::removeSurfaceMesh(name);
    }

    auto* m2 = polyscope::registerSurfaceMesh(name, V, F);

    // global id scalar quantity ( V order = ids order)
    std::vector<double> gid;
    gid.reserve(ids.size());
    for (auto id : ids) gid.push_back((double)id);
    m2->addVertexScalarQuantity("global id", gid);

    m2->setEnabled(true);
    m2->setSurfaceColor(glm::vec3(0.6f, 0.8f, 1.0f));
}


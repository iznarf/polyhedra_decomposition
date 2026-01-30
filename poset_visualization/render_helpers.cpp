#include "render_helpers.h"

#include "geometry_utils.h"           
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h> 

#include <algorithm>
#include <limits>



using glm::vec3;

namespace viz_helpers {

    std::vector<std::array<int,3>>
    faces_from_triangles(const df::Tri2& t,
                        const std::unordered_map<df::vertex_id,int>& to_local) {
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

    // gives vertices the global id in polyscope as scalar quantity
    void add_global_id_quantity(polyscope::SurfaceMesh* mesh,
                                const std::vector<df::vertex_id>& ids) {
        std::vector<double> values;
        values.reserve(ids.size());
        for (auto id : ids) {
            values.push_back(static_cast<double>(id));
        }
        mesh->addVertexScalarQuantity("global id", values);
    }

    std::vector<vec3>
    make_planar_poset_vertices(const std::vector<df::vertex_id>& ids,
                            const std::vector<df::P2>& points2d,
                            float cx, float cz, float scale) {
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

    std::vector<vec3>
    make_lifted_poset_vertices(const std::vector<df::vertex_id>& ids,
                            const std::vector<df::P2>& points2d,
                            float cx, float cz, float scale_xy, float scale_z) {
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


    // we need this for bounding box computation 
    // bbox is for triangle mesh visualization 
    glm::vec3 bbox_center(const std::vector<glm::vec3>& V) {
        if (V.empty()) return glm::vec3(0.f);
        glm::vec3 lo = V[0], hi = V[0];
        for (const auto& p : V) {
            lo = glm::min(lo, p);
            hi = glm::max(hi, p);
        }
        return 0.5f * (lo + hi);
    }

    void recenter(std::vector<glm::vec3>& V, const glm::vec3& c) {
        for (auto& p : V) p -= c;
    }

    // gives each mesh the node index as a scalar quantity 
    void add_node_id_quantity(polyscope::SurfaceMesh* mesh, int nodeIndex) {
        const size_t nV = mesh->nVertices();
        std::vector<double> values(nV, (double)nodeIndex);
        mesh->addVertexScalarQuantity("poset node", values);
    }



} // namespace viz_helpers

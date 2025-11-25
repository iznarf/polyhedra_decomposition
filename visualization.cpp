#include "visualization.h"
#include "input.h"
#include "debug.h"

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <glm/glm.hpp>
#include <unordered_map>
#include <vector>
#include <array>
#include <imgui.h>

using df::DebugTetrahedron;
using df::DebugTetKind;

using df::vertex_id;
using glm::vec3;

namespace  {



// vertcies for planar triangulation z = 0
std::vector<vec3> points_planar(const std::vector<vertex_id>& ids, const std::vector<df::P2>& P) {
    std::vector<vec3> V; V.reserve(ids.size());
    for (auto id : ids) {
        const auto& p = P[id];
        V.push_back(vec3((float)p.x(), 0.0f, (float)p.y()));
    }
    return V;
}

// vertices for lifted triangulation z = x^2 + y^2 
std::vector<vec3> points_lifted(const std::vector<vertex_id>& ids, const std::vector<df::P2>& P) {
    std::vector<vec3> V; V.reserve(ids.size());
    for (auto id : ids) {
        const auto& p = P[id];
        float z = (float)(p.x()*p.x() + p.y()*p.y());
        V.push_back(vec3((float)p.x(), z, (float)p.y()));
    }
    return V;
}

// faces as local index triples
std::vector<std::array<int,3>> faces_from_triangles(const df::Tri2& t, const std::unordered_map<vertex_id,int>& to_local) {
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


// we use this function to give each vertex the global index as a scalar quantity in polyscope
static void add_global_id_quantity(polyscope::SurfaceMesh* mesh,
                                   const std::vector<df::vertex_id>& ids) {
    // polyscope scalar quantities are doubles, so just cast
    std::vector<double> values;
    values.reserve(ids.size());
    for (auto id : ids) {
        values.push_back(static_cast<double>(id));
    }

    mesh->addVertexScalarQuantity("global id", values);
}

// one debug mesh for the tetrahedra with geometry + ids + type
struct DebugTetGeom {
    std::array<vertex_id,4> ids;
    std::array<vec3,4>      pos;
    DebugTetKind            kind;
};

std::vector<DebugTetGeom> g_debug_tets;
int g_debug_idx = 0;
polyscope::SurfaceMesh* g_debug_mesh = nullptr;

// faces of a tetrahedron in canonical order
const std::vector<std::array<size_t,3>> kTetFaces = {
    {0,1,2},
    {0,3,1},
    {1,3,2},
    {0,2,3}
};


// update/create the single debug mesh for the current index
void update_debug_tet_mesh() {
    if (g_debug_tets.empty())
        return;

    const DebugTetGeom& dt = g_debug_tets[g_debug_idx];

    std::vector<vec3> V = {
        dt.pos[0], dt.pos[1], dt.pos[2], dt.pos[3]
    };

    if (!g_debug_mesh) {
        g_debug_mesh = polyscope::registerSurfaceMesh("debug tetra", V, kTetFaces);
        //g_debug_mesh->setTransparency(0.5f); // make it transparent
    } else {
        g_debug_mesh->updateVertexPositions(V);
    }

    // debug tets are yellow 
    glm::vec3 ochre(0.87f, 0.69f, 0.17f);
;

    g_debug_mesh->setSurfaceColor(ochre);


}

} // anonymous namespace





namespace viz {

static inline glm::vec3 lift_paraboloid(const df::P2& p) {
  double x = CGAL::to_double(p.x());
  double y = CGAL::to_double(p.y());
  return glm::vec3((float)x, (float)(x*x + y*y), (float)y); // y-up
}

// collect the global indices present in a triangulation
std::vector<vertex_id> present_ids(const df::Tri2& t) {
    std::vector<vertex_id> ids;
    ids.reserve(t.number_of_vertices());
    for (auto v = t.finite_vertices_begin(); v != t.finite_vertices_end(); ++v)
        ids.push_back(v->info());
    return ids;
}

// map global index to local compact index [0,..,V-1] for polyscope
std::unordered_map<vertex_id,int> make_local_index(const std::vector<vertex_id>& ids) {
    std::unordered_map<vertex_id,int> m;
    m.reserve(ids.size());
    for (int i = 0; i < (int)ids.size(); ++i) m[ids[i]] = i;
    return m;
}

// register one triangulation as (2D mesh + lifted mesh)
// we use this function to register both lower and upper triangulations at startup
void register_triangulation_as_mesh(const df::Tri2& tri,
                         const std::vector<df::P2>& points2d,
                         const std::string& name_planar,
                         const std::string& name_lifted) {
    auto ids      = viz::present_ids(tri);
    auto to_local = viz::make_local_index(ids);

    auto vertices_2d = points_planar(ids, points2d);
    auto faces       = faces_from_triangles(tri, to_local);



    glm::vec3 lightBlue(0.6f, 0.8f, 1.0f);  

    // planar mesh
    auto* m2 = polyscope::registerSurfaceMesh(name_planar, vertices_2d, faces);
    add_global_id_quantity(m2, ids);
    // do not show 2D mesh by default
    m2->setEnabled(false); 
    m2->setSurfaceColor(lightBlue);
   


    // lifted mesh
    auto vertices_3d = points_lifted(ids, points2d);
    auto* m3 = polyscope::registerSurfaceMesh(name_lifted, vertices_3d, faces);
    add_global_id_quantity(m3, ids);
    m3->setSurfaceColor(lightBlue);

    //m3->setTransparency(0.5); 
}





// we use this function to register and to update the current triangulation
void show_or_update_current(const df::InputData& D) {
    // build fresh buffers from the current triangulation
    const auto& tri = D.tri_current;

    const auto ids      = present_ids(tri);
    const auto to_local = make_local_index(ids);
    const auto V2       = points_planar(ids,  D.points2d);
    const auto V3       = points_lifted(ids,  D.points2d);
    const auto F        = faces_from_triangles(tri, to_local);

    // remove old meshes only if they exist
    if (polyscope::hasSurfaceMesh("current 2D")) {
        polyscope::getSurfaceMesh("current 2D")->remove();
    }
    if (polyscope::hasSurfaceMesh("current lifted")) {
        polyscope::getSurfaceMesh("current lifted")->remove();
    }

    glm::vec3 aubergine(0.549f, 0.227f, 0.459f);


    // (re)register meshes
    auto* m2 =polyscope::registerSurfaceMesh("current 2D", V2, F);
    add_global_id_quantity(m2, ids);
    // do not show 2D mesh by default
    m2->setEnabled(false); 
    m2->setSurfaceColor(aubergine);

    auto* m3 = polyscope::registerSurfaceMesh("current lifted", V3, F);
    add_global_id_quantity(m3, ids);
    m3->setSurfaceColor(aubergine);
    m3->setTransparency(0.6f);
}






// --------------------------------
// debug tetrahedra visualization
// --------------------------------

void load_debug_tetrahedra(const df::InputData& D,
                           const std::vector<DebugTetrahedron>& tets)
{
    g_debug_tets.clear();
    g_debug_idx  = 0;

    
    if (polyscope::hasSurfaceMesh("debug tetra")) {
        polyscope::getSurfaceMesh("debug tetra")->remove();
        g_debug_mesh = nullptr;
    }

    g_debug_tets.reserve(tets.size());

    for (const auto& t : tets) {
        DebugTetGeom g;
        g.ids  = t.verts;
        g.kind = t.kind;

        // turn 4 vertex ids into 4 lifted positions
        for (int i = 0; i < 4; ++i) {
            vertex_id vid = t.verts[i];
            const df::P2& p2 = D.points2d[vid];
            g.pos[i] = lift_paraboloid(p2);
        }

        g_debug_tets.push_back(g);
    }

    if (!g_debug_tets.empty()) {
        update_debug_tet_mesh();
    } else {
        std::cout << "load_debug_tetrahedra: no tets to visualize\n";
    }
}

void debug_tet_ui() {
    if (g_debug_tets.empty()) {
        ImGui::Text("no debug tetrahedra");
        return;
    }

    int n = static_cast<int>(g_debug_tets.size());
    ImGui::Text("debug tetrahedra (yellow): %d", n);
    ImGui::Text("current index: %d", g_debug_idx);

    // show vertex ids of the current tet
    const auto& cur = g_debug_tets[g_debug_idx];
    ImGui::Text("vertices: (%d, %d, %d, %d)",
                (int)cur.ids[0], (int)cur.ids[1],
                (int)cur.ids[2], (int)cur.ids[3]);

    ImGui::Text("kind: %s",
        cur.kind == DebugTetKind::EdgeFlip ? "2-2 flip" : "1-3 flip");

    // Buttons
    if (ImGui::Button("prev") && g_debug_idx > 0) {
        --g_debug_idx;
        update_debug_tet_mesh();
    }
    ImGui::SameLine();
    if (ImGui::Button("next") && g_debug_idx < n - 1) {
        ++g_debug_idx;
        update_debug_tet_mesh();
    }

}
// --------------------------------
// end debug tetrahedra visualization


// --------------------------------
// replay visualization (needs fixing)
// --------------------------------

// update visualization of replay mesh
void show_or_update_replay(const df::InputData& D) {
    const auto& tri = D.tri_replay;

    const auto ids      = present_ids(tri);
    const auto to_local = make_local_index(ids);
    const auto V2       = points_planar(ids,  D.points2d);
    const auto V3       = points_lifted(ids,  D.points2d);
    const auto F        = faces_from_triangles(tri, to_local);

    // remove old replay meshes if they exist
    if (polyscope::hasSurfaceMesh("replay 2D")) {
        polyscope::getSurfaceMesh("replay 2D")->remove();
    }
    if (polyscope::hasSurfaceMesh("replay lifted")) {
        polyscope::getSurfaceMesh("replay lifted")->remove();
    }

    // i need neon green for replay meshes
    glm::vec3 neon_green(0.0f, 1.0f, 0.0f);
    

    auto* m2 = polyscope::registerSurfaceMesh("replay 2D", V2, F);
    add_global_id_quantity(m2, ids);
   
    m2->setSurfaceColor(neon_green);
    auto* m3 = polyscope::registerSurfaceMesh("replay lifted", V3, F);
    add_global_id_quantity(m3, ids);
  
    m3->setSurfaceColor(neon_green);
    m3->setTransparency(0.5f);
}




} // namespace viz



#include "envelope.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/Env_surface_data_traits_3.h>
#include <CGAL/envelope_3.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>

#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>


#include <polyscope/polyscope.h>
#include <polyscope/curve_network.h>

#include <glm/glm.hpp>
#include <array>
#include <map>
#include <vector>
#include <string>




#include <map>
#include <queue>
#include <stdexcept>

namespace env_max {

namespace {

using EK = CGAL::Exact_predicates_exact_constructions_kernel;

// ---------------------------
// Envelope traits with surface data
// ---------------------------
using BaseTraits = CGAL::Env_triangle_traits_3<EK>;
using Traits     = CGAL::Env_surface_data_traits_3<BaseTraits, TriTag, TriTag>; 
using Surface_3   = Traits::Surface_3;
using Xy_surface_3 = Traits::Xy_monotone_surface_3;
using Envelope_diagram_2 = CGAL::Envelope_diagram_2<Traits>;

using P2 = EK::Point_2;
using P3 = EK::Point_3;

// ---------------------------
// helper to triangulate polygon
// ---------------------------
struct FaceInfo {
  int nesting_level = -1;
  bool in_domain() const { return nesting_level % 2 == 1; }
};

using Vb  = CGAL::Triangulation_vertex_base_2<EK>;
using Fb0 = CGAL::Constrained_triangulation_face_base_2<EK>;              
using Fb  = CGAL::Triangulation_face_base_with_info_2<FaceInfo, EK, Fb0>; 
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using Itag = CGAL::Exact_intersections_tag;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<EK, Tds, Itag>;

// mark faces of constrained delaunay triangulation which are inside the boundary polygon 
static void mark_domains(CDT& cdt){
    for (auto f = cdt.all_faces_begin(); f != cdt.all_faces_end(); ++f) {
        f->info().nesting_level = -1;
    }

    std::queue<CDT::Face_handle> q;

    // start from infinite face (outside)
    CDT::Face_handle inf = cdt.infinite_face();
    inf->info().nesting_level = 0;
    q.push(inf);

    while (!q.empty()) {
        CDT::Face_handle fh = q.front();
        q.pop();

        for (int i = 0; i < 3; ++i) {
        CDT::Face_handle n = fh->neighbor(i);
        if (n->info().nesting_level != -1) continue;

        // crossing a constrained edge increases nesting level
        const bool constrained = cdt.is_constrained(CDT::Edge(fh, i));
        n->info().nesting_level = fh->info().nesting_level + (constrained ? 1 : 0);
        q.push(n);
        }
    }
}

// extract the outer boundary as ccw point cycle from an envelope face
// returns empty if unbounded or no outer boundary
static std::vector<P2> face_outer_polygon(const Envelope_diagram_2::Face& f){
    if (f.is_unbounded() || !f.has_outer_ccb()) return {};

    std::vector<P2> poly;
    auto circ = f.outer_ccb();
    if (circ == nullptr) return {};

    auto start = circ;
    do {
        // source vertex point of the halfedge
        const auto& vp = circ->source()->point();
        poly.push_back(vp);
        ++circ;
    } while (circ != start);

    // remove duplicate last==first if it happens 
    if (poly.size() >= 2 && poly.front() == poly.back()) poly.pop_back();
    return poly;
}

static P3 eval_on_triangle_plane(const EK::Triangle_3& tri, const P2& p){
  // terrain assumption: triangle is not vertical, so plane has non-zero z coefficient
  EK::Plane_3 pl(tri.vertex(0), tri.vertex(1), tri.vertex(2));

  const EK::FT a = pl.a();
  const EK::FT b = pl.b();
  const EK::FT c = pl.c();
  const EK::FT d = pl.d();

  if (c == EK::FT(0)) {
    // not a terrain patch; fall back to projecting to closest point in plane at z=0 (not correct for envelope).
    // for our use case (lifted planar triangulations) this should not happen
    throw std::runtime_error("envelope: vertical triangle plane encountered (c==0)");
  }

  const EK::FT x = p.x();
  const EK::FT y = p.y();
  const EK::FT z = (-d - a * x - b * y) / c;

  return P3(x, y, z);
}

// Map 2D exact point -> output vertex index to use them as keys in map 
struct Less_xy {
  EK::Less_xy_2 less;
  bool operator()(const P2& a, const P2& b) const { return less(a, b); }
};



static int get_or_add_vertex(const P2& p2,
                            const EK::Triangle_3& tri,
                            std::map<P2, int, Less_xy>& vmap,
                            std::vector<std::array<double, 3>>& V){
    auto it = vmap.find(p2);
    if (it != vmap.end()) return it->second;

    const P3 p3 = eval_on_triangle_plane(tri, p2);

    const int idx = (int)V.size();
    V.push_back({CGAL::to_double(p3.x()), CGAL::to_double(p3.y()), CGAL::to_double(p3.z())});
    vmap.emplace(p2, idx);
    return idx;
}


// Helpers to visualize envelope results in polyscope
static int get_or_add_2d_point_vertex(
    const P2& p,
    std::map<P2, int, Less_xy>& vmap,
    std::vector<glm::vec3>& V,
    float y_plane,
    float x_offset) {
    auto it = vmap.find(p);
    if (it != vmap.end()) return it->second;

    const float x = (float)CGAL::to_double(p.x());
    const float z = (float)CGAL::to_double(p.y());

    const int idx = (int)V.size();
    V.push_back(glm::vec3(x + x_offset, y_plane, z));
    vmap.emplace(p, idx);
    return idx;
}


void register_envelope_diagram_edges_2d(
    const std::string& name,
    const Envelope_diagram_2& diag,
    glm::vec3 color,
    float radius,
    float y_plane,
    float x_offset,
    float z_offset,
    float scale)
{
    std::vector<glm::vec3> V;
    std::vector<std::array<int, 2>> E;
    std::map<P2, int, Less_xy> vmap;

    auto get_vid = [&](const P2& p) -> int {
        auto it = vmap.find(p);
        if (it != vmap.end()) return it->second;

        const double x = CGAL::to_double(p.x());
        const double z = CGAL::to_double(p.y()); // (x,y) in diagram == (x,z) in world

        const int idx = (int)V.size();
        V.push_back(glm::vec3((float)x, y_plane, (float)z)); // store UNTRANSFORMED first
        vmap.emplace(p, idx);
        return idx;
    };

    // halfedges -> edges
    for (auto he = diag.halfedges_begin(); he != diag.halfedges_end(); ++he) {
        if (he->is_fictitious()) continue;
        if (!(he < he->twin())) continue; // keep one direction

        const P2 s = he->source()->point();
        const P2 t = he->target()->point();

        const int a = get_vid(s);
        const int b = get_vid(t);
        if (a != b) E.push_back({a, b});
    }

    // center in XZ, then scale, then translate to (x_offset, z_offset)
    if (!V.empty()) {
        glm::vec2 c(0.f);
        for (auto& p : V) c += glm::vec2(p.x, p.z);
        c /= (float)V.size();

        for (auto& p : V) {
            glm::vec2 q = glm::vec2(p.x, p.z) - c;
            q *= scale;
            p.x = q.x + x_offset;
            p.z = q.y + z_offset;
            p.y = y_plane;
        }
    }

    if (polyscope::hasCurveNetwork(name)) polyscope::removeStructure(name);
    auto* cn = polyscope::registerCurveNetwork(name, V, E);
    cn->setColor(color);
    cn->setRadius(radius, false);
    cn->setEnabled(true);
}

}



bool compute_and_register_lower_envelope_diagram_edges_2d(
    const std::string& name,
    const std::vector<InputTriangle>& tris,
    glm::vec3 color,
    float radius,
    float y_plane,
    float x_offset,
    float z_offset,
    float scale
) {
    if (tris.empty()) {
        if (polyscope::hasCurveNetwork(name)) polyscope::removeStructure(name);
        return true;
    }

    try {
        // 1) build CGAL surfaces with tags
        std::vector<Surface_3> surfaces;
        surfaces.reserve(tris.size());

        for (const auto& t : tris) {
            const P3 a(EK::FT(t.p0[0]), EK::FT(t.p0[1]), EK::FT(t.p0[2]));
            const P3 b(EK::FT(t.p1[0]), EK::FT(t.p1[1]), EK::FT(t.p1[2]));
            const P3 c(EK::FT(t.p2[0]), EK::FT(t.p2[1]), EK::FT(t.p2[2]));
            EK::Triangle_3 tri(a, b, c);
            surfaces.emplace_back(tri, t.tag);
        }

        // 2) compute lower envelope diagram (planar subdivision)
        Envelope_diagram_2 diag;
        CGAL::lower_envelope_3(surfaces.begin(), surfaces.end(), diag);

        // 2.5) compute centroid of all diagram edge endpoints (in diagram coords)
        double cx = 0.0, cy = 0.0;
        int cnt = 0;
        for (auto he = diag.halfedges_begin(); he != diag.halfedges_end(); ++he) {
            if (!he->is_fictitious() && he < he->twin()) {
                const P2 s = he->source()->point();
                const P2 t = he->target()->point();

                cx += CGAL::to_double(s.x()); cy += CGAL::to_double(s.y()); cnt++;
                cx += CGAL::to_double(t.x()); cy += CGAL::to_double(t.y()); cnt++;
            }
        }
        if (cnt > 0) { cx /= (double)cnt; cy /= (double)cnt; }

        // 3) convert diagram edges to a Polyscope curve network in the XZ plane
        std::vector<glm::vec3> V;
        std::vector<std::array<int, 2>> E;
        std::map<P2, int, Less_xy> vmap;

        auto get_or_add_vertex = [&](const P2& p) -> int {
            auto it = vmap.find(p);
            if (it != vmap.end()) return it->second;

            const double px = CGAL::to_double(p.x());
            const double py = CGAL::to_double(p.y());

            // center -> scale -> translate into grid cell
            const float x = (float)((px - cx) * (double)scale) + x_offset;
            const float z = (float)((py - cy) * (double)scale) + z_offset;

            const int idx = (int)V.size();
            V.emplace_back(x, y_plane, z);
            vmap.emplace(p, idx);
            return idx;
        };

        for (auto he = diag.halfedges_begin(); he != diag.halfedges_end(); ++he) {
            if (!he->is_fictitious() && he < he->twin()) {
                const P2 s = he->source()->point();
                const P2 t = he->target()->point();

                const int i0 = get_or_add_vertex(s);
                const int i1 = get_or_add_vertex(t);
                if (i0 != i1) E.push_back({ i0, i1 });
            }
        }

        // replace old
        if (polyscope::hasCurveNetwork(name)) polyscope::removeStructure(name);

        auto* cn = polyscope::registerCurveNetwork(name, V, E);
        cn->setColor(color);
        cn->setRadius(radius, false);
        cn->setEnabled(true);

        return true;
    } catch (...) {
        return false;
    }
}







// output is triangle mesh with out.V = list of 3D vertices 
// out.F = list of triangle faces as indices into out.V
// out.tri_tag = per triangle face tag {mesh_id, face_id} indicating from which input triangle it was induced

bool compute_lower_envelope(const std::vector<InputTriangle>& tris, Mesh& out){
    out = Mesh{};
    if (tris.empty()) return true;

    try {
        // 1) convert every triangle into CGAL EK::Triangle_3 with attached triangle tag = {mesh_id, face_id}
        std::vector<Surface_3> surfaces;
        surfaces.reserve(tris.size());


        for (const auto& t : tris) {
        const P3 a(EK::FT(t.p0[0]), EK::FT(t.p0[1]), EK::FT(t.p0[2]));
        const P3 b(EK::FT(t.p1[0]), EK::FT(t.p1[1]), EK::FT(t.p1[2]));
        const P3 c(EK::FT(t.p2[0]), EK::FT(t.p2[1]), EK::FT(t.p2[2]));
        EK::Triangle_3 tri(a, b, c);

        // BaseTraits::Surface_3 is EK::Triangle_3
        // we wtap triangle_3 into surface_3 (with data) because envelope_3 works with surfaces
        const Surface_3 s(tri, t.tag);
        surfaces.push_back(s);
        }

        // 2) compute envelope diagram
        // what is the diagram: project all triangle surfaces plane and compute where they intersect
        // result is a planar subdivision whose faces are labeled by the triangle that has maximum height on that region
        // diagram now contains the planar subdivision of the upper envelope projection with all faces labeled by their inducing triangle surface
        Envelope_diagram_2 diag;
        CGAL::lower_envelope_3(surfaces.begin(), surfaces.end(), diag);


        // show planar subdivision edges
        register_envelope_diagram_edges_2d("envelope diagram", diag,glm::vec3(1,0,0),0.005f,0.0f, 0.0f, 0.0f, 1.0f);


        // 3) For each bounded face: triangulate its 2D region and lift on inducing triangle
        // iterate over all envelope regions, take boundary polygon and read which input triangle it is 

        // prepare vertex map 
        std::map<P2, int, Less_xy> vmap;

        // loop over faces in diagram
        for (auto fit = diag.faces_begin(); fit != diag.faces_end(); ++fit) {
        const auto& f = *fit;
        if (f.is_unbounded()) continue;
        if (f.number_of_surfaces() == 0) continue; // no height winner stored 
        
        // get outer boundary polygon of face
        // important here: faces are not triangles but polygons
        std::vector<P2> poly = face_outer_polygon(f); // extract boundary polygon as ccw point cycle 
        // we skip it degenerate 
        if (poly.size() < 3) continue;

        
        // we get the inducing surface: we want the trianguation surface patch tag to know from where it comes 
        const Xy_surface_3 xys = f.surface(); // triangle surface patch that is highest on this face
        const TriTag tag = xys.data(); // data tells us which triangulation/node and which face in it


        // xys inherits BaseTraits xy-surface, which represents triangle patch
        // in Env_triangle_traits_3 base surface type is Triangle_3 and xy-surface can be converted to it
        // we store it by static_cast to the base
        const BaseTraits::Xy_monotone_surface_3& base_xys = xys;
        const EK::Triangle_3 base_tri = base_xys; // relies on conversion provided by traits




        // triangulate polygon via constrained delaunay triangulation 
        // CHANGE: use Tri2, not CDT 
        CDT cdt; // constrained delaunay triangulation
        std::vector<CDT::Vertex_handle> vhs; 
        vhs.reserve(poly.size());

        // vertices are inserted in order, we get them from boundary polygons
        for (const auto& p : poly) {
            vhs.push_back(cdt.insert(p));
        }

        for (size_t i = 0; i < vhs.size(); ++i) {
            cdt.insert_constraint(vhs[i], vhs[(i + 1) % vhs.size()]);
        }

        // otherwise cdt would fill whole plane, we only want ínterior region of polygon 
        mark_domains(cdt);

        // extract triangles from cdt that are inside the polygon
        for (auto tf = cdt.finite_faces_begin(); tf != cdt.finite_faces_end(); ++tf) {
            if (!tf->info().in_domain()) continue;

            // read triangles P2 
            const P2 p0 = tf->vertex(0)->point();
            const P2 p1 = tf->vertex(1)->point();
            const P2 p2 = tf->vertex(2)->point();

            // add verticees or get existing onees to avoid duplicates
            const int i0 = get_or_add_vertex(p0, base_tri, vmap, out.V);
            const int i1 = get_or_add_vertex(p1, base_tri, vmap, out.V);
            const int i2 = get_or_add_vertex(p2, base_tri, vmap, out.V);

            // add face with triangle tag
            out.F.push_back({i0, i1, i2});
            out.tri_tag.push_back({tag.mesh_id, tag.face_id});
        }
        }

        return true;
    } catch (...) {
        return false;
    }
}


} // namespace env_max
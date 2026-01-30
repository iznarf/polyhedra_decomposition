#include "poset_tri_mesh.h"

#include "input.h"
#include "poset.h"
#include "replay.h"

#include "render_helpers.h"
#include "poset_vis_helpers.h"
#include "poset_utils.h" 
#include "visualization.h"

#include <polyscope/polyscope.h>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

namespace pst_vis {

void remove_triangulation_vis(TriVisResult& R) {
  for (auto* m : R.meshes2d) if (m) m->remove();
  for (auto* m : R.meshes3d) if (m) m->remove();
  R.meshes2d.clear();
  R.meshes3d.clear();
  R.centers.clear();
}

void register_poset_as_triangulation_P1(const std::string& name_prefix,
                                        const df::InputData& D,
                                        const pst::Poset1& P1,
                                        TriVisResult& out)
{
  remove_triangulation_vis(out);

  const int n = (int)P1.nodes.size();
  if (n <= 0) return;

  out.meshes2d.assign(n, nullptr);
  out.meshes3d.assign(n, nullptr);
  out.centers.assign(n, glm::vec3(0.f));

  out.pivot2.assign(n, glm::vec3(0.f));
  out.pivot3.assign(n, glm::vec3(0.f));


  const float TRI_SCALE_PLAN  = 0.7f;
  const float TRI_SCALE_LIFT  = 0.7f;
  const float LIFT_HEIGHT_SCL = 0.5f;

  for (int node_idx = 0; node_idx < n; ++node_idx) {
    // build triangulation for node
    df::Tri2 tri = D.tri_poset;
    pst::replay_history_poset(tri, P1.nodes[node_idx].history, D);

    // convert to polyscope mesh
    auto ids      = viz::present_ids(tri);
    auto to_local = viz::make_local_index(ids);
    auto faces    = viz_helpers::faces_from_triangles(tri, to_local);

    auto V2 = viz_helpers::make_planar_poset_vertices(ids, D.points2d, 0.f, 0.f, TRI_SCALE_PLAN);
    auto V3 = viz_helpers::make_lifted_poset_vertices(ids, D.points2d, 0.f, 0.f, TRI_SCALE_LIFT, LIFT_HEIGHT_SCL);

    glm::vec3 c3 = viz_helpers::bbox_center(V3);
    // store pivot (center of this mesh in local coords)
    out.pivot3[node_idx] = viz_helpers::bbox_center(V3);
    out.pivot2[node_idx] = viz_helpers::bbox_center(V2);

    std::string name2d = name_prefix + std::to_string(node_idx) + " 2D";
    std::string name3d = name_prefix + std::to_string(node_idx) + " 3D";

    auto* m2 = polyscope::registerSurfaceMesh(name2d, V2, faces);
    viz_helpers::add_node_id_quantity(m2, node_idx); // gives mesh the node index as scalar quantity
    viz_helpers::add_global_id_quantity(m2, ids); // gives mesh vertices the global id as scalar quantity
    m2->setEnabled(true);
    m2->setSurfaceColor(glm::vec3(0.6f, 0.8f, 1.0f));
    m2->setEdgeWidth(1.0f);
    m2->setEdgeColor(glm::vec3(0, 0, 0));

    auto* m3 = polyscope::registerSurfaceMesh(name3d, V3, faces);
    viz_helpers::add_node_id_quantity(m3, node_idx); // gives mesh the node index as scalar quantity
    viz_helpers::add_global_id_quantity(m3, ids); // gives mesh vertices the global id as scalar quantity
    m3->setEnabled(true);
    m3->setSurfaceColor(glm::vec3(0.2f, 0.4f, 0.8f));
    m3->setTransparency(0.6f);
    m3->setEdgeWidth(1.0f);
    m3->setEdgeColor(glm::vec3(0, 0, 0));

    out.meshes2d[node_idx] = m2;
    out.meshes3d[node_idx] = m3;
  }
}
void apply_triangulation_centers(TriVisResult& R,
                                 const std::vector<glm::vec3>& centers,
                                 float uniformScale)
{
    const float liftY = 0.05f; // visual offset for the 3D mesh only (UP axis!)

    if (centers.size() != R.meshes3d.size()) {
        std::cerr << "[poset_tri_mesh] centers size mismatch: centers="
                  << centers.size() << " nodes=" << R.meshes3d.size() << "\n";
        return;
    }
    if (R.pivot2.size() != R.meshes3d.size()) {
        std::cerr << "[poset_tri_mesh] pivot2 size mismatch: pivot2="
                  << R.pivot2.size() << " nodes=" << R.meshes3d.size() << "\n";
        return;
    }

    R.centers = centers;

    for (int i = 0; i < (int)centers.size(); ++i) {

        // Use planar pivot so the 3D lift stays "above" the plane visually
        glm::mat4 S    = glm::scale(glm::mat4(1.f), glm::vec3(uniformScale));
        glm::mat4 Tpiv = glm::translate(glm::mat4(1.f), -R.pivot2[i]);

        // 2D at center
        if (R.meshes2d[i]) {
            glm::mat4 Tpos2 = glm::translate(glm::mat4(1.f), centers[i]);
            R.meshes2d[i]->setTransform(Tpos2 * S * Tpiv);
        }

        // 3D lifted in Y (UP)
        if (R.meshes3d[i]) {
            glm::vec3 liftedCenter = centers[i] + glm::vec3(0.f, liftY, 0.f);
            glm::mat4 Tpos3 = glm::translate(glm::mat4(1.f), liftedCenter);
            R.meshes3d[i]->setTransform(Tpos3 * S * Tpiv);
        }
    }
}





} // namespace pst_vis





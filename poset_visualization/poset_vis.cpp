#include "poset_vis.h"

#include "poset2.h"
#include "poset_utils.h"
#include "poset_vis_helpers.h"
#include "poset_tri_mesh.h"

#include "dedekind_poset.h"


#include <polyscope/polyscope.h>
#include <polyscope/point_cloud.h>
#include <polyscope/curve_network.h>

#include <glm/glm.hpp>
#include <vector>
#include <string>
#include <iostream>



namespace pst_vis {

pst_vis::TriVisResult g_triVis_P1;   //global for poset1 visualizations
pst_vis::TriVisResult g_triVis_P2;  //global for poset2 visualizations


std::vector<glm::vec3> g_gridPos_P1; // global for poset1 grid positions
std::vector<glm::vec3> g_gridPos_P2; // global for poset2 grid positions

float g_gridXSpacing_P1 = 0.2f;
float g_gridZSpacing_P1 = 0.2f;
float g_gridXSpacing_P2 = 0.2f;
float g_gridZSpacing_P2 = 0.2f;




void register_cover_up_as_edges(const std::string& name, const std::vector<glm::vec3>& pos, const std::vector<std::vector<int>>& cover_up, glm::vec3 color) {
    const int n = (int)pos.size();
    if ((int)cover_up.size() != n) {
        std::cerr << "[poset_vis] cover_up.size() != pos.size()\n";
        return;
    }

    // build edge list (u -> v)
    std::vector<std::array<int, 2>> edges;
    edges.reserve(n * 2);

    for (int u = 0; u < n; ++u) {
        for (int v : cover_up[u]) {
            if (0 <= v && v < n) {
                edges.push_back({u, v});
            }
        }
    }

    // replace old network if exists
    if (polyscope::hasCurveNetwork(name)) {
        polyscope::removeStructure(name);
    }

    auto* cn = polyscope::registerCurveNetwork(name, pos, edges);
    cn->setRadius(0.002f, false);   // tweak later
    cn->setColor(color);
    cn->setEnabled(true);
}




std::vector<glm::vec3> compute_grid_pos(const std::vector<int>& levels,glm::vec3 center,float xSpacing,float zSpacing){
    auto byLevel = pst_vis::group_nodes_by_level(levels);
    return pst_vis::grid_for_poset_nodes(byLevel, (int)levels.size(), center, xSpacing, zSpacing);
}

// register poset nodes as a point cloud laid out on a grid according to levels
// give points in point cloud node ids as a scalar quantity, as well as levels
void register_poset_as_grid(const std::string& name,
                            int number_of_nodes,
                            const std::vector<int>& levels,
                            glm::vec3 center,
                            float xSpacing,
                            float zSpacing,
                            glm::vec3 color)
{
    if ((int)levels.size() != number_of_nodes) {
        std::cerr << "[poset_vis] levels.size() != number_of_nodes\n";
        return;
    }


    auto pos = compute_grid_pos(levels, center, xSpacing, zSpacing);



    if (polyscope::hasPointCloud(name))
        polyscope::removeStructure(name);

    auto* pc = polyscope::registerPointCloud(name, pos);
    pc->setPointRadius(0.01f, false);
    pc->setPointColor(color);

    // level quantity (int -> double is fine, polyscope stores scalars as double internally)
    pc->addScalarQuantity("level", levels);

    // node id quantity (0..n-1)
    std::vector<double> nodeIds;
    nodeIds.reserve(number_of_nodes);
    for (int i = 0; i < number_of_nodes; ++i) nodeIds.push_back((double)i);
    pc->addScalarQuantity("node id", nodeIds);
}


void register_poset1_as_triangulations(const df::InputData& D, const pst::Poset1& P1, const glm::vec3& center, float xSpacing, float zSpacing){
    remove_triangulation_vis(g_triVis_P1);

    // 1) build meshes
    register_poset_as_triangulation_P1("P1 node ", D, P1, g_triVis_P1);

    const int n = (int)P1.nodes.size();
    if (n == 0) return;

    // 2) grid layout based on levels
    auto pos = compute_grid_pos(P1.levels, center, xSpacing, zSpacing);

    // 3) store grid positions globally
    pst_vis::g_gridXSpacing_P1 = xSpacing;
    pst_vis::g_gridZSpacing_P1 = zSpacing;


    // 4) apply transforms from the grid 
    float meshScale = 0.08f;
    apply_triangulation_centers(g_triVis_P1, pos, meshScale);

    std::cout << "[poset_vis] registered " << n << " P1 triangulation nodes placed on grid\n";

}


void register_poset2_as_triangulations(const df::InputData& D, const pst::Poset1& P1, const pst2::Poset2& P2, const glm::vec3& center, float xSpacing, float zSpacing){
        remove_triangulation_vis(g_triVis_P2);

    // build meshes (independent polyscope objects via different names)
    register_poset_as_triangulation_P1("P2 node ", D, P1, g_triVis_P2);

    const int n = (int)P1.nodes.size();
    if (n == 0) return;

    // grid layout based on levels of P2
    auto pos = compute_grid_pos(P2.levels, center, xSpacing, zSpacing);

    // store grid positions globally
    pst_vis::g_gridXSpacing_P2 = xSpacing;
    pst_vis::g_gridZSpacing_P2 = zSpacing;

    // mesh scale scales the size of the triangulation meshes 
    float meshScale = 0.08f;
    apply_triangulation_centers(g_triVis_P2, pos, meshScale);

    std::cout << "[poset_vis] registered " << n << " P2 triangulation nodes placed on grid\n";
}


// function to show/hide all trivis meshes in a TriVisResult
void set_trivis_enabled(const pst_vis::TriVisResult& R, bool show2d, bool show3d){
    for (auto* m : R.meshes2d) if (m) m->setEnabled(show2d);
    for (auto* m : R.meshes3d) if (m) m->setEnabled(show3d);
}

void register_poset1(const df::InputData& D, const pst::Poset1& P1,const glm::vec3& center, float xSpacing, float zSpacing, glm::vec3 color){
    std::cout << "[poset_vis] register_poset1 called. P1.levels.size()=" << P1.levels.size() << "\n";
    
    // grid first (for camera framing + toggles)
    // red color for P1
    register_poset_as_grid("P1 grid nodes", (int)P1.levels.size(), P1.levels, center, xSpacing, zSpacing, color);

    // triangulations
    register_poset1_as_triangulations(D, P1, center, xSpacing, zSpacing);

    // edges
    auto pos = compute_grid_pos(P1.levels, center, xSpacing, zSpacing);
    pst_vis::g_gridPos_P1 = pos; // store grid positions globally
    std::cout << "[poset_vis] g_gridPos_P1 set. size=" << pst_vis::g_gridPos_P1.size()
          << "  addr=" << (void*)&pst_vis::g_gridPos_P1 << "\n";
    register_cover_up_as_edges("P1 edges", pos, P1.cover_up, color);
}

void register_poset2(const df::InputData& D, const pst::Poset1& P1, const pst2::Poset2& P2,const glm::vec3& center, float xSpacing, float zSpacing, glm::vec3 color){
    std::cout << "[poset_vis] register_poset2 called. P2.levels.size()=" << P2.levels.size() << "\n";

    
    // grid first
    // green color for P2
    register_poset_as_grid("P2 grid nodes", (int)P2.levels.size(), P2.levels, center, xSpacing, zSpacing, color);

    register_poset2_as_triangulations(D, P1, P2, center, xSpacing, zSpacing);

    auto pos = compute_grid_pos(P2.levels, center, xSpacing, zSpacing);
    pst_vis::g_gridPos_P2 = pos; // store grid positions globally
    std::cout << "[poset_vis] g_gridPos_P2 set. size=" << pst_vis::g_gridPos_P2.size()
          << "  addr=" << (void*)&pst_vis::g_gridPos_P2 << "\n";

    register_cover_up_as_edges("P2 edges", pos, P2.cover_up, color);
}

// ---------------------------------------------------------------------
// register poset as a colored grid (per-node colors), returns computed positions
// this is for the completion: new nodes are dark blue, old nodes are in the color used elsewhere
std::vector<glm::vec3> register_poset_as_grid_colored(
    const std::string& name,
    int number_of_nodes,
    const std::vector<int>& levels,
    const glm::vec3& center,
    float xSpacing,
    float zSpacing,
    const std::vector<glm::vec3>& colors,
    float pointRadius
) {
    if ((int)levels.size() != number_of_nodes) {
        std::cerr << "[poset_vis] levels.size() != number_of_nodes\n";
        return {};
    }
    if ((int)colors.size() != number_of_nodes) {
        std::cerr << "[poset_vis] colors.size() != number_of_nodes\n";
        return {};
    }

    // same placement logic as everywhere else:
    // levels -> group by level -> grid positions
    auto byLevel = group_nodes_by_level(levels);
    auto pos = grid_for_poset_nodes(byLevel, number_of_nodes, center, xSpacing, zSpacing);

    if (polyscope::hasPointCloud(name))
        polyscope::removeStructure(name);

    auto* pc = polyscope::registerPointCloud(name, pos);
    pc->setPointRadius(pointRadius, false);

    // scalar quantities
    pc->addScalarQuantity("level", levels);

    std::vector<double> nodeIds;
    nodeIds.reserve(number_of_nodes);
    for (int i = 0; i < number_of_nodes; ++i) nodeIds.push_back((double)i);
    pc->addScalarQuantity("node id", nodeIds);

    // per-node colors
    auto* cq = pc->addColorQuantity("old/new", colors);
    cq->setEnabled(true);

    return pos;
}





} // namespace pst_vis













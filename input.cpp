#include "input.h"
#include <CGAL/Convex_hull_traits_adapter_2.h> // for convex hull with property map
#include <CGAL/property_map.h> // for make_property_map
#include <CGAL/enum.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <algorithm>
#include <numeric>    // std::iota
#include <random>
#include <cmath>
#include <cassert>


namespace df {

// checks if three points are collinear
bool are_collinear(const P2& a, const P2& b, const P2& c) {
    // use CGAL orientation test
    return CGAL::orientation(a, b, c) == CGAL::COLLINEAR;
}

// checks the hull point set and removes collinear points
// modifies hull_ids in place
static void prune_collinear_hull(const std::vector<df::P2>& P, std::vector<std::size_t>& hull_ids) {
    if (hull_ids.size() < 3) return;
        std::vector<std::size_t> out; out.reserve(hull_ids.size());
        auto n = hull_ids.size();
    //iterate over hull ids and check if the previous, current, 
    // next are collinear, modulo because we want a closed loop
    for (size_t i = 0; i < n; ++i) {
        const auto ia = hull_ids[(i + n - 1) % n];
        const auto ib = hull_ids[i];
        const auto ic = hull_ids[(i + 1) % n];
        // check if ia, ib, ic are collinear
        // if not, ib is a true corner, keept it
        // if collinear, ib lies on a straight line between ia and ic, skip it
        if (!are_collinear(P[ia], P[ib], P[ic])) out.push_back(ib);
    }
    // replace old hull_ids with filtered one
    hull_ids.swap(out);
    // if hull contains less then 3 points, choose another input and print warning 
    if (hull_ids.size() < 3) {
        std::cerr << "[Warning] Hull degenerated after pruning: "
                  << hull_ids.size()
                  << " vertices remain (need >= 3). ";
    }
}

// generate n 2d points inside a disk of radius R, uniformly over the area to spread points evenly
// we pick random angle and random radius for polar coordinates and then convert to cartesian coordinates
static std::vector<P2> sample_points_in_disk(int n, double R, std::mt19937& rng) { // // std::mt19937& rng: random number generator
    std::uniform_real_distribution<double> Uang(0.0, 2.0 * PI); // random angle
    std::uniform_real_distribution<double> Urad(0.0, 1.0); // random radius
    std::vector<P2> pts; pts.reserve(n);
    for (int i = 0; i < n; ++i) {
        double ang = Uang(rng); // random angle between 0 and 2pi
        double r   = R * std::sqrt(Urad(rng));     // sqrt of random radius
        pts.emplace_back(r * std::cos(ang), r * std::sin(ang)); // convert to cartesian
    }  
    return pts;
}

// function to build edge list from triangulation
static std::vector<std::array<df::vertex_id, 2>>
build_edge_list(const df::Tri2& tri)
{
    std::vector<std::array<df::vertex_id, 2>> edges;
    edges.reserve(3 * tri.number_of_faces()); // upper bound

    for (auto e = tri.finite_edges_begin(); e != tri.finite_edges_end(); ++e) {
        auto fh = e->first;
        int i   = e->second; // edge opposite vertex i
        auto a = fh->vertex(tri.cw(i))->info();
        auto b = fh->vertex(tri.ccw(i))->info();
        edges.push_back({a, b});   
    }

    return edges;
}

// function to build face list from triangulation
static std::vector<std::array<df::vertex_id, 3>> build_face_list(const df::Tri2& tri) {
    std::vector<std::array<df::vertex_id, 3>> F;
    F.reserve(tri.number_of_faces());
    for (auto f = tri.finite_faces_begin(); f != tri.finite_faces_end(); ++f) {
        auto v0 = f->vertex(0)->info();
        auto v1 = f->vertex(1)->info();
        auto v2 = f->vertex(2)->info();
        // enforce a canonical ordering
        F.push_back({ v0, v1, v2 });
    }
    return F;
}


// create random input: points, hull triangulation, full triangulation
// n_points: number of vertices 
// seed: starting value for the random number generator -> if seed stays the same, we can get the same random points
InputData make_random_input(int n_points, unsigned seed) {
    assert(n_points >= 3); // need at least 3 points to form a triangulation
    std::mt19937 rng(seed); 

    InputData D;
    // 1) random points with interior points
    D.points2d = sample_points_in_disk(n_points, /*R=*/1.0, rng);

    // 2) make an global index array [0,...,n-1]
    std::vector<std::size_t> global_indices(D.points2d.size());
    std::iota(global_indices.begin(), global_indices.end(), 0);

    // 3) build convex hull: We use the property map to then give CGAL just the global indices for the convex hull
    // since we have the property map, CGAL can then build the convex hull based on the indices
    // build property map: global index -> point: for integer k it gives back points2d[k]
    // in O(1) time
    auto point_map = CGAL::make_property_map(D.points2d);

    // 4) ask CGAL for the hull *indices*
    std::vector<std::size_t> hull_ids;
    hull_ids.reserve(D.points2d.size());

    // 5) convex hull as global indices in ccw order 
    using hull_adapter = CGAL::Convex_hull_traits_adapter_2<K, decltype(point_map)>;
    // this gives back the global indices of the hull points in ccw order
    CGAL::convex_hull_2(global_indices.begin(), global_indices.end(),
                        std::back_inserter(hull_ids),
                        hull_adapter(point_map));

    // prune collinear points from hull
    prune_collinear_hull(D.points2d, hull_ids);


    // 6) build (point, global index) pairs for the hull
    std::vector<std::pair<P2, vertex_id>> hull_pairs;
    hull_pairs.reserve(hull_ids.size());
    for (auto id : hull_ids)
        hull_pairs.emplace_back(D.points2d[id], id);
    
    // 7) insert into triangulation; each vertex’s info() becomes the global index
    D.tri_upper.clear();
    D.tri_upper.insert(hull_pairs.begin(), hull_pairs.end());
    D.tri_current.insert(hull_pairs.begin(), hull_pairs.end());

  
    using vertex_handle = Tri2::Vertex_handle;

    // global index -> vertex handle map for upper triangulation
    D.index_to_vertex_handle_upper.assign(D.points2d.size(), vertex_handle{});
    for (auto v = D.tri_upper.finite_vertices_begin();
        v != D.tri_upper.finite_vertices_end(); ++v)
    D.index_to_vertex_handle_upper[v->info()] = v; 

    // global index -> vertex handle map for current triangulation
    D.index_to_vertex_handle_current.assign(D.points2d.size(), vertex_handle{});
    for (auto v = D.tri_current.finite_vertices_begin();
        v != D.tri_current.finite_vertices_end(); ++v)
    D.index_to_vertex_handle_current[v->info()] = v;
    
    // indices of all points for shuffling
    // we shuffle the order of insertion to generate random triangulations
    std::vector<std::size_t> indices(D.points2d.size());
    std::iota(indices.begin(), indices.end(), 0);

    // shuffle the indices (not the points)
    std::shuffle(indices.begin(), indices.end(), rng);

    // now insert (point, id) pairs using the shuffled order
    std::vector<std::pair<P2, std::size_t>> shuffled_pairs;
    shuffled_pairs.reserve(D.points2d.size());
    for (auto id : indices)
        shuffled_pairs.emplace_back(D.points2d[id], id);

    D.tri_lower.clear();
    D.tri_lower.insert(shuffled_pairs.begin(), shuffled_pairs.end());
 
    // global index -> vertex handle map for lower triangulation
    D.index_to_vertex_handle_lower.assign(D.points2d.size(), vertex_handle{});
    for (auto v = D.tri_lower.finite_vertices_begin();
        v != D.tri_lower.finite_vertices_end(); ++v)
    D.index_to_vertex_handle_lower[v->info()] = v; 


    return D;
}

} // namespace df



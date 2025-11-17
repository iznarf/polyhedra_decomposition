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
#include <iostream>


namespace df {

    using DT2 = CGAL::Delaunay_triangulation_2<K, Tri_ds>;



// check if "p_new" is collinear with ANY pair of points in "pts"
static bool forms_collinear_triple(const std::vector<P2>& pts, const P2& p_new) {
    std::size_t n = pts.size();
    if (n < 2) return false; // need at least 2 existing points to form a triple

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (CGAL::orientation(pts[i], pts[j], p_new) == CGAL::COLLINEAR) {
                return true;
            }
        }
    }
    return false;
}




// generate n 2d points inside a disk of radius R, uniformly over the area,
// with the extra constraint: no three points are collinear
static std::vector<P2> sample_points_in_disk(int n, double R, std::mt19937& rng) {
    std::uniform_real_distribution<double> Uang(0.0, 2.0 * PI); // random angle
    std::uniform_real_distribution<double> Urad(0.0, 1.0);      // random radius

    std::vector<P2> pts;
    pts.reserve(n);

    for (int k = 0; k < n; ++k) {
        P2 p;
        // rejection sampling until p does not form a collinear triple
        while (true) {
            double ang = Uang(rng);                 // angle in [0, 2pi)
            double r   = R * std::sqrt(Urad(rng));  // radius with area-uniform sampling
            p = P2(r * std::cos(ang), r * std::sin(ang));

            if (!forms_collinear_triple(pts, p)) {
                break; // accept this point
            }
            // else: try again
        }

        pts.push_back(p);
    }

    return pts;
}






// create random input: points, hull triangulation, full triangulation
// n_points: number of vertices 
// seed: starting value for the random number generator -> if seed stays the same, we can get the same random points
InputData make_random_input(int n_points, unsigned seed) {
    assert(n_points >= 3); // need at least 3 points to form a triangulation
    std::mt19937 rng(seed); 

    InputData D;
    // clear flip history from earlier runs
    D.flip_history.clear(); 
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

    return D;
}

} // namespace df



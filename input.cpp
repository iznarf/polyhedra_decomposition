#include "input.h"
#include <CGAL/Convex_hull_traits_adapter_2.h> // for convex hull with property map
#include <CGAL/property_map.h> // for make_property_map
#include <CGAL/enum.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Object.h>

#include <algorithm>
#include <numeric>    
#include <random>
#include <cmath>
#include <cassert>
#include <iostream>
#include <CGAL/intersections.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Segment_3.h>



namespace df {


// check if the point p_new is collinear with any pair of points in pts
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

    // 4) ask CGAL for the hull indices
    std::vector<std::size_t> hull_ids;
    // build hull edge list
    hull_ids.reserve(D.points2d.size());

    // 5) convex hull as global indices in ccw order 
    using hull_adapter = CGAL::Convex_hull_traits_adapter_2<K, decltype(point_map)>;
    // this gives back the global indices of the hull points in ccw order
    CGAL::convex_hull_2(global_indices.begin(), global_indices.end(),
                        std::back_inserter(hull_ids),
                        hull_adapter(point_map));

    // build hull edge list for later use
    int H = hull_ids.size();
    D.hull_edges.clear();
    for(int i = 0; i < H; ++i) {
        vertex_id a = hull_ids[i];
        vertex_id b = hull_ids[(i + 1) % H]; // circular
            // this needs fixing, we do not want to order the indices here
            if (a > b) std::swap(a, b);          
            D.hull_edges.emplace_back(a, b);

    }

    // 6) build (point, global index) pairs for the hull
    std::vector<std::pair<P2, vertex_id>> hull_pairs;
    hull_pairs.reserve(hull_ids.size());
    for (auto id : hull_ids)
        hull_pairs.emplace_back(D.points2d[id], id);
    
    // 7) insert into triangulation; each vertex’s info() becomes the global index
    D.tri_upper.clear();
    D.tri_upper.insert(hull_pairs.begin(), hull_pairs.end());
    D.tri_current.insert(hull_pairs.begin(), hull_pairs.end());

    // we need this for applying then all recorded flips
    D.tri_replay.clear();
    D.tri_replay.insert(hull_pairs.begin(), hull_pairs.end());

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


// fix: put it in geometry helper
static inline P3 lift(const P2& p) {
    return P3(p.x(), p.y(), p.x() * p.x() + p.y() * p.y());
}

bool is_hull_edge(const InputData& D, vertex_id u, vertex_id v)
{
    if (u > v) std::swap(u, v);

    for (auto& e : D.hull_edges) {
        if (e.first == u && e.second == v)
            return true;
    }
    return false;
}

// checks if the lifted surfaces of lower and upper triangulations intersect
// returns true if they do intersect since this is invalid input for the algorithm
bool lifted_triangulations_intersect(const InputData& D)
{
    const Tri2& upper = D.tri_upper;
    const Tri2& lower = D.tri_lower;

    using Triangle3 = CGAL::Triangle_3<K>;
    using Segment3  = CGAL::Segment_3<K>;

    // Loop over all finite faces in upper
    for (auto fu = upper.finite_faces_begin(); fu != upper.finite_faces_end(); ++fu) {

        df::vertex_id au = fu->vertex(0)->info();
        df::vertex_id bu = fu->vertex(1)->info();
        df::vertex_id cu = fu->vertex(2)->info();

        const P3 a3u = lift(D.points2d[au]);
        const P3 b3u = lift(D.points2d[bu]);
        const P3 c3u = lift(D.points2d[cu]);

        Triangle3 triangle_upper(a3u, b3u, c3u);

        std::array<df::vertex_id,3> up_ids = { au, bu, cu };

        // Loop over all finite faces in lower
        for (auto fl = lower.finite_faces_begin(); fl != lower.finite_faces_end(); ++fl) {

            df::vertex_id al = fl->vertex(0)->info();
            df::vertex_id bl = fl->vertex(1)->info();
            df::vertex_id cl = fl->vertex(2)->info();

            // make lifted triangle out of it
            const P3 a3l = lift(D.points2d[al]);
            const P3 b3l = lift(D.points2d[bl]);
            const P3 c3l = lift(D.points2d[cl]);
            Triangle3 trianle_lower(a3l, b3l, c3l);

            std::array<df::vertex_id,3> low_ids = { al, bl, cl };
           
            // Count common vertices between these two faces (in 2D).
            int common = 0;
            std::vector<df::vertex_id> common_ids;
            for (auto iu : up_ids) {
                for (auto il : low_ids) {
                    if (iu == il) {
                        ++common;
                        common_ids.push_back(iu);
                    }
                }
            }

            // we still have to check the CGAL intersection type
            auto intersection_object = CGAL::intersection(triangle_upper, trianle_lower);
            if (!intersection_object) {
                continue; // no intersection
            }

            if (intersection_object) {
                if  (const CGAL::Point_3<K>* s = std::get_if<CGAL::Point_3<K>>(&*intersection_object)) {
                    if(common == 1) {
                    // intersection is a point -> allowed
                    continue;
                    }
                }
                if (const CGAL::Segment_3<K>* s = std::get_if<CGAL::Segment_3<K>>(&*intersection_object)) {
                    if (common == 1) {
                        // we have to check if point is mislabeled as segment
                        // check distance of source and target of segment to the 3 lifted points of upper triangle
                        const P3& p0 = s->source();
                        const P3& p1 = s->target();
                        double d0 = std::min({ CGAL::squared_distance(p0, a3u),
                                            CGAL::squared_distance(p0, b3u),
                                            CGAL::squared_distance(p0, c3u) });
                        double d1 = std::min({ CGAL::squared_distance(p1, a3u),
                                            CGAL::squared_distance(p1, b3u),
                                            CGAL::squared_distance(p1, c3u) });
                        const double eps = 1e-12;
                        if (d0 < eps && d1 < eps) {
                            // both endpoints are very close to lifted upper triangle vertices
                            // treat as point intersection
                            continue;
                        }
                        else {
                            // intersection is a segment but triangles only share one vertex -> invalid
                            std::cout << "[intersect] lifted triangles intersect in segment but only share one vertex\n";
                            return true;
                        }
                    }
                    if (common == 2){ 
                        // check if the common points are forming a hull edge
                        df::vertex_id u = common_ids[0];
                        df::vertex_id v = common_ids[1];
                        // now check if they form a hull edge
                        if (!is_hull_edge(D, u, v)) {
                            return true;
                        }
                        else {
                            // intersection is a hull edge -> allowed
                            continue;
                        }  
                    }    
                    if (common == 0 || common == 3) {
                        // intersection is a segment but triangles do not share an edge -> invalid
                        std::cout << "[intersect] lifted triangles intersect in segment but do not share an edge\n";
                        return true;
                    }
                }
                if (const Triangle3* t = std::get_if<Triangle3>(&*intersection_object)) {
                    if(common != 3){
                        std::cout << "[intersect] lifted triangles share an intersection face\n";
                        return true; // allowed if they are the same triangle
                    }
                }
                // we have to check if intersection object is a vector of points (coplanar triangles)
                if (const std::vector<CGAL::Point_3<K>>* pts = std::get_if<std::vector<CGAL::Point_3<K>>>(&*intersection_object)) {
                    return true; // intersection is a face or edge or multiple points -> invalid   
                }  
            }
        }
    }
    // no interior intersections found
    return false;
}

InputData make_random_valid_input(int n_points, unsigned seed_start)
{
    unsigned seed = seed_start;

    while (true) {
        InputData D = make_random_input(n_points, seed);

        if (!lifted_triangulations_intersect(D)) {
            std::cout << "[input] using seed " << seed
                      << " (lifted triangulations are non-intersecting)\n";
            return D;
        }


        std::cout << "[input] seed " << seed
                  << " gives intersecting lifts, trying next seed...\n";
        ++seed;
    }
}

} // namespace df





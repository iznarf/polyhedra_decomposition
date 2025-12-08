#include "input.h"
#include "geometry_utils.h"
#include <CGAL/Convex_hull_traits_adapter_2.h> // for convex hull with property map
#include <CGAL/property_map.h> // for make_property_map
#include <CGAL/enum.h>
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

struct Convert_vertex_RT_to_Tri2
{
    mutable bool first_vertex = true;  // the first one is the infinite vertex

    df::Tri2::Vertex operator()(const df::Tri2Regular::Vertex&) const {
        return df::Tri2::Vertex();
    }

    void operator()(const df::Tri2Regular::Vertex& src,
                    df::Tri2::Vertex& tgt) const
    {
        if (first_vertex) {
            // this is the infinite vertex: no point/info to set
            first_vertex = false;
            return;
        }

        // src.point() is a Weighted_point_2<K>, get its center point()
        tgt.point() = src.point().point();   // type: df::P2
        tgt.info()  = src.info();            // keep vertex_id
    }
};

struct Convert_face_RT_to_Tri2
{
    df::Tri2::Face operator()(const df::Tri2Regular::Face&) const {
        return df::Tri2::Face();
    }

    void operator()(const df::Tri2Regular::Face&,
                    df::Tri2::Face&) const
    {
        // nothing to copy, faces don't store geometry
    }
};






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



// return true iff all lifted points are vertices of their convex hull
// ensure that all points are in convex position in 3D
bool vertices_in_convex_position(const InputData& D)
{
    const std::size_t n = D.points2d.size();
    if (n == 0) return true;
    if (n <= 3) return true; // any less or equal than 3 points in 3D are trivially in convex position

    // 1) lift all points to paraboloid 
    std::vector<P3> lifted;
    lifted.reserve(n);
    for (const P2& p : D.points2d) {
        lifted.push_back(lift(p));
    }

    // 2) compute convex hull in 3D
    Polyhedron_3 hull;
    CGAL::convex_hull_3(lifted.begin(), lifted.end(), hull);

    // 3) convex hull vertices are a subset of the input points
    // if every input point is extreme, these two cardinalities are the same 
    const std::size_t n_hull_vertices = hull.size_of_vertices();

    return n_hull_vertices == n;
}



// create random input: points, hull triangulation, full triangulation
// n_points: number of vertices 
// seed: starting value for the random number generator -> if seed stays the same, we can get the same random points
InputData make_random_input(int n_points, unsigned seed) {
    assert(n_points >= 3); // need at least 3 points to form a triangulation
    std::mt19937 rng(seed); 

    InputData D;
    // 1) random points with interior points

    
    //D.points2d = sample_points_in_disk(n_points, 1.0, rng);
    /*
    //print 2d points list
    std::cout << "Generated " << n_points << " random 2D points:\n";
    for (std::size_t i = 0; i < D.points2d.size(); ++i) {
        const P2& p = D.points2d[i];
        std::cout << " " << p.x() << " " << p.y() << "\n";
    }
    */

    
    
    

    
    
    //example points 1

    D.points2d = {
        P2(-2,1), //lift to z = 2
        P2(2,1),   //lift to z = 2
        P2(0,4),    //lift to z = 2
        P2(0.1,3),    //lift to z = 1
        P2(-0.8,1.7),   //lift to z = 1
        P2(1,1.5)     //lift to z = 1
    };

    
    

    // 2) make an global index array [0,...,n-1]
    std::vector<std::size_t> global_indices(D.points2d.size());
    std::iota(global_indices.begin(), global_indices.end(), 0);

    // build weighted points for regular triangulation
    D.points2d_weighted.clear();
    for (std::size_t i = 0; i < D.points2d.size(); ++i) {
        const P2& p = D.points2d[i];
        K::FT w = -(p.x() * p.x() + p.y() * p.y());
        P2_weighted wp(p, w);

        auto vh = D.tri_regular.insert(wp);
        vh->info() = i; // set global index as info
        D.points2d_weighted.push_back(wp); 
    }

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

    
    // 7) insert into triangulation; each vertex info() becomes global index
    D.tri_upper.clear();
    D.tri_upper.insert(hull_pairs.begin(), hull_pairs.end());

    D.tri_current.clear();
    D.tri_current.insert(hull_pairs.begin(), hull_pairs.end());

    // we need this for applying then all recorded flips
    D.tri_replay.clear();
    D.tri_replay.insert(hull_pairs.begin(), hull_pairs.end());

    // we need this for building the poset
    D.tri_poset.clear();
    D.tri_poset.insert(hull_pairs.begin(), hull_pairs.end());

    using vertex_handle = Tri2::Vertex_handle;
    
    // indices of all points 
    std::vector<std::size_t> indices(D.points2d.size());
    std::iota(indices.begin(), indices.end(), 0);

    // 8) build full triangulation with all points
    
    std::vector<std::pair<P2, std::size_t>> lower_pairs;
    lower_pairs.reserve(D.points2d.size());
    for (auto id : indices)
        lower_pairs.emplace_back(D.points2d[id], id);
    
    
    D.tri_lower.clear();
    //D.tri_lower.insert(lower_pairs.begin(), lower_pairs.end());

    
    
    // make lower triangulation be the (unweighted) copy of the regular triangulation
    df::Convert_vertex_RT_to_Tri2 cv;
    df::Convert_face_RT_to_Tri2   cf;

    auto inf_v =
        D.tri_lower.tds().copy_tds(
            D.tri_regular.tds(),
            D.tri_regular.infinite_vertex(),
            cv, cf
        );

    D.tri_lower.set_infinite_vertex(inf_v);
    CGAL_assertion(D.tri_lower.is_valid());
    

    
    return D;
}



// checks if an edge is part of the convex hull in 2D
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
// general fix here: we have to avoid the exact intersection cases, we want to check as much as possible using geometry predicates
bool lifted_triangulations_intersect(const InputData& D)
{
    const Tri2& upper = D.tri_upper;
    const Tri2& lower = D.tri_lower;

   
    using Segment3  = CGAL::Segment_3<K>;

    // Loop over all finite faces in upper
    for (auto fu = upper.finite_faces_begin(); fu != upper.finite_faces_end(); ++fu) {

        df::vertex_id au = fu->vertex(0)->info();
        df::vertex_id bu = fu->vertex(1)->info();
        df::vertex_id cu = fu->vertex(2)->info();

        const P3 a3u = df::lift(D.points2d[au]);
        const P3 b3u = df::lift(D.points2d[bu]);
        const P3 c3u = df::lift(D.points2d[cu]);

        CGAL::Triangle_3<K> triangle_upper(a3u, b3u, c3u);

        std::array<df::vertex_id,3> up_ids = { au, bu, cu };

        // Loop over all finite faces in lower
        for (auto fl = lower.finite_faces_begin(); fl != lower.finite_faces_end(); ++fl) {

            df::vertex_id al = fl->vertex(0)->info();
            df::vertex_id bl = fl->vertex(1)->info();
            df::vertex_id cl = fl->vertex(2)->info();

            // make lifted triangle out of it
            const P3 a3l = df::lift(D.points2d[al]);
            const P3 b3l = df::lift(D.points2d[bl]);
            const P3 c3l = df::lift(D.points2d[cl]);
            CGAL::Triangle_3<K> triangle_lower(a3l, b3l, c3l);

            std::array<df::vertex_id,3> low_ids = { al, bl, cl };
           
            // count common vertices between these two faces (in 2D)
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
            auto intersection_object = CGAL::intersection(triangle_upper, triangle_lower);
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
                if (const CGAL::Triangle_3<K>* t = std::get_if<CGAL::Triangle_3<K>>(&*intersection_object)) {
                    std::cout << "[intersect] lifted triangles share an intersection face\n";
                    return true; // intersection is a face -> invalid
                }
                // we have to check if intersection object is a vector of points (coplanar triangles)
                if (const std::vector<CGAL::Point_3<K>>* pts = std::get_if<std::vector<CGAL::Point_3<K>>>(&*intersection_object)) {
                    return true; // intersection is a face or edge or multiple points -> invalid   
                }  
            }
        }
    }
    // loop over edges in upper triangulation
    for (auto eu = upper.finite_edges_begin(); eu != upper.finite_edges_end(); ++eu) {
        Tri2::Face_handle fu = eu->first;
        int idxu = eu->second;

        df::vertex_id au = fu->vertex(upper.cw(idxu))->info();
        df::vertex_id bu = fu->vertex(upper.ccw(idxu))->info();

        const P3 a3u = df::lift(D.points2d[au]);
        const P3 b3u = df::lift(D.points2d[bu]);

        Segment3 seg_upper(a3u, b3u);

        // loop over edges in lower triangulation
        for (auto el = lower.finite_edges_begin(); el != lower.finite_edges_end(); ++el) {
            Tri2::Face_handle fl = el->first;
            int idxl = el->second;

            df::vertex_id al = fl->vertex(lower.cw(idxl))->info();
            df::vertex_id bl = fl->vertex(lower.ccw(idxl))->info();

            const P3 a3l = df::lift(D.points2d[al]);
            const P3 b3l = df::lift(D.points2d[bl]);

            Segment3 seg_lower(a3l, b3l);



            // check intersection of these two segments
            if (CGAL::do_intersect(seg_upper, seg_lower)) {
                // check if they share endpoints and how many and which ones since we have to check if the edge is hull edge 
                int counter = 0;
                std::vector<df::vertex_id> common_ids_edge;
                
                if (au == al || au == bl) {
                    ++counter;
                    common_ids_edge.push_back(au);
                }
                if (bu == al || bu == bl) {
                    ++counter;
                    common_ids_edge.push_back(bu);
                }
                // check if the edge is hull edge if they share two endpoints
                if (counter == 2) {
                    df::vertex_id u = common_ids_edge[0];
                    df::vertex_id v = common_ids_edge[1];
                    // now check if they form a hull edge
                    if (!is_hull_edge(D, u, v)) {
                        return true;
                    }
                    else {
                        // intersection is a hull edge -> allowed
                        continue;
                    }
                }
                if (counter == 1) {
                    // intersection is a point -> allowed
                    continue;
                }
                else {
                    // invalid intersection
                    std::cout << "[intersect] lifted edges intersect\n";
                    return true;
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

        if (!df::vertices_in_convex_position(D)) {
        std::cout << "[check] vertices are NOT in convex position\n";
        } else {
            std::cout << "[check] vertices ARE in convex position\n";
        }

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





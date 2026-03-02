#include "input_ear_star.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>



std::vector<LabeledPoint> make_convex_ngon_with_star(int n, double R = 1.0) {
    if (n < 3) throw std::runtime_error("n must be >= 3");

    std::vector<LabeledPoint> out;
    out.reserve(n + 1);
    df::vertex_id id = 0;

    // boundary points on a circle: strictly convex, CCW order by construction
    const double two_pi = 2.0 * std::acos(-1.0);
    for (int i = 0; i < n; ++i) {
        id = i; // boundary vertices get ids 0..n-1
        double ang = two_pi * double(i) / double(n);
        double x = R * std::cos(ang);
        double y = R * std::sin(ang);
        out.push_back({id, df::P2(x, y), std::to_string(id)});
    }

    // interior point, hardcoded position for now 
    id = n; // star vertex gets id -1
    out.push_back({id, df::P2(-0.3, 0.2), "*"});

    return out;
}

df::Tri2 make_triangulation(const std::vector<LabeledPoint>& points) {
    df::Tri2 tri;
    for (const auto& lp : points) {
        auto vh = tri.insert(lp.p);
        vh->info() = lp.id; // set vertex info to the global id of this point
    }
    for (auto v = tri.finite_vertices_begin(); v != tri.finite_vertices_end(); ++v) {
    std::cout << "vertex info=" << v->info() << "\n";
}
    return tri;
}



InputData_ear_star make_input_ear_star(int n, double R) {
    InputData_ear_star out;
    out.points = make_convex_ngon_with_star(n, R);
    out.tri_start = make_triangulation(out.points);
    return out;
}

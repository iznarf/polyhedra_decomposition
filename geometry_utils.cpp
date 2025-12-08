#include "geometry_utils.h"
#include <CGAL/enum.h>

namespace df {

bool quad_strictly_convex(const P2& a, const P2& b,
                          const P2& c, const P2& d)
{
    // 1) C and D on opposite sides of AB
    auto oc = CGAL::orientation(a, b, c);
    auto od = CGAL::orientation(a, b, d);

    if (oc == CGAL::COLLINEAR || od == CGAL::COLLINEAR || oc == od)
        return false;

    // 2) A and B on opposite sides of CD
    auto oa = CGAL::orientation(c, d, a);
    auto ob = CGAL::orientation(c, d, b);

    if (oa == CGAL::COLLINEAR || ob == CGAL::COLLINEAR || oa == ob)
        return false;

    // if both conditions hold, the diagonals intersect in the interior
    // and the quadrilateral is strictly convex 
    return true;
}




// we take the triangle (a,b,c) in 2D and make sure it is ccw oriented
// then we compute the 3D orientation of (a3,b3,c3,d3)
// if it is positive then d3 lies in the normal direction of (a3,b3,c3): convetion "above plane"
// if it is negative then d3 lies in the opposite direction of the normal "below plane"


CGAL::Orientation oriented_height_sign(
    const P2& a2, const P2& b2, const P2& c2,
    const P3& a3, const P3& b3, const P3& c3, const P3& d3)
{
    auto o2 = CGAL::orientation(a2, b2, c2);

    // make a,b,c ccw in 2D 
    if (o2 == CGAL::RIGHT_TURN) {
        // swap b and c in 3D
        return CGAL::orientation(a3, c3, b3, d3);
    }
    if (o2 == CGAL::COLLINEAR) {
        return CGAL::COLLINEAR;
    }
    // ccw case
    return CGAL::orientation(a3, b3, c3, d3);
}



/*
static inline double debug_height(double x, double y) {
    constexpr double eps = 1e-9;
    auto eq = [eps](double a, double b) { return std::abs(a - b) < eps; };

    // z = 2  first 3 points
    if (eq(x, -2.0) && eq(y, 1.0)) return 2.0;
    if (eq(x,  2.0) && eq(y, 1.0)) return 2.0;
    if (eq(x,  0.0) && eq(y, 4.0)) return 2.0;

    // z = 1 last 3 points
    if (eq(x,  0.1) && eq(y, 3.0)) return 1.0;
    if (eq(x, -0.8) && eq(y, 1.7)) return 1.0;
    if (eq(x,  1.0) && eq(y, 1.5)) return 1.0;

    return 0.0; // fallback
}


P3 lift(const P2& p) {
    double x = CGAL::to_double(p.x());
    double y = CGAL::to_double(p.y());
    double z = debug_height(x, y);
    return P3(x, y, z);
}

P3 lift_regular(const P2_weighted& p) {
    double x = CGAL::to_double(p.x());
    double y = CGAL::to_double(p.y());
    double z = debug_height(x, y);
    return P3(x, y, z);
}


*/
P3 lift(const P2& p) {
    return P3(p.x(), p.y(), p.x()*p.x() + p.y()*p.y());
}

P3 lift_regular(const P2_weighted& p) {
    return P3(p.x(), p.y(), p.x()*p.x() + p.y()*p.y());
}









}




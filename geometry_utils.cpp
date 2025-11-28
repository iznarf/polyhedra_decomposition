#include "geometry_utils.h"
#include <CGAL/enum.h>

namespace df {

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


inline P3 lift(const P2& p) {
    return P3(p.x(), p.y(), p.x()*p.x() + p.y()*p.y());
}

inline P3 lift_regular(const P2_weighted& p) {
    return P3(p.x(), p.y(), p.x()*p.x() + p.y()*p.y());
}


}




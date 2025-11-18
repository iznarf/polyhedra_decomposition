#pragma once

#include "input.h"         


namespace df {

// Lifted height sign: returns orientation of (a', b', c', d')
// after ensuring ccw(a,b,c) in 2D.
CGAL::Orientation oriented_height_sign(
    const P2& a2, const P2& b2, const P2& c2,
    const P3& a3, const P3& b3, const P3& c3, const P3& d3);


} // namespace df

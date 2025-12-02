#pragma once
#include "input.h"  

namespace df { 

// Compare heights of two lifted segments (c,d) and (u,v) at their 2D intersection point with CGAL intersection
int compare_heights_at_intersection(const df::P2& c2, const df::P2& d2,
                                    const df::P2& u2, const df::P2& v2,
                                    const df::P3& c3, const df::P3& d3,
                                    const df::P3& u3, const df::P3& v3);




bool height_test_orientation_based(const P2& c2, const P2& d2,
                                   const P2& u2, const P2& v2,
                                   const P3& c3, const P3& d3,
                                   const P3& u3, const P3& v3);
} // namespace df


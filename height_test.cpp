#include "height_test.h"
#include "input.h"

#include <CGAL/Segment_2.h>
#include <CGAL/intersections.h>
#include <CGAL/number_utils.h>   
#include <CGAL/squared_distance_2.h>

namespace df {

using FT       = K::FT; // number type used by kernel 
using Segment2 = CGAL::Segment_2<K>;



// computation of parameter t such that p = (1 - t) * a + t * b in 2D
static inline FT compute_param_t(const P2& p, const P2& a, const P2& b) {
    // we can compute t using x or y coordinate, depending on which is more stable
    const FT dx = b.x() - a.x();
    const FT dy = b.y() - a.y();

    // choose the more stable coordinate to divide by
    // for vertical segments dx can be very close to 0
    // for horizontal segments dy can be very close to 0
    if (CGAL::abs(dx) >= CGAL::abs(dy)) {
        // if the segments intersect, dx != 0 in this branch
        return (p.x() - a.x()) / dx;
    } else {
        // likewise, dy != 0 here
        return (p.y() - a.y()) / dy;
    }
}

// linear interpolation of z along a 3D segment for parameter t
static inline FT interpolate_z(const P3& a3, const P3& b3, const FT& t)
{
    return (FT(1) - t) * a3.z() + t * b3.z();
}

// approximate equality of two 2D points
static inline bool approx_equal_p2(const P2& a, const P2& b)
{
    // epsilon on squared distance
    const FT eps_sq = FT(1e-12);
    return CGAL::squared_distance(a, b) < eps_sq;
}



// Compare heights of two lifted segments (c,d) and (u,v) at their 2D intersection point
// return -1 if (c,d) is below (u,v)
// return +1 if (c,d) is above (u,v)
// return 0 if equal height / coplanar / overlapping segment in 2D (should not happen here)
int compare_heights_at_intersection(const P2& c2, const P2& d2, const P2& u2, const P2& v2, const P3& c3, const P3& d3, const P3& u3, const P3& v3) {
    Segment2 seg_cd(c2, d2);
    Segment2 seg_uv(u2, v2);

    auto inter = CGAL::intersection(seg_cd, seg_uv);
    const FT eps_height = FT(1e-12);

  
    // intersection is a single 2D point
    if (const P2* p = std::get_if<P2>(&*inter)) {
        const P2& ip = *p;

        FT t_cd = compute_param_t(ip, c2, d2);
        FT t_uv = compute_param_t(ip, u2, v2);

        FT z_cd = interpolate_z(c3, d3, t_cd);
        FT z_uv = interpolate_z(u3, v3, t_uv);

        if (z_cd + eps_height < z_uv) return -1; // (c,d) below (u,v)
        if (z_cd > z_uv + eps_height) return +1; // (c,d) above (u,v)
        return 0;                                // equal within tolerance
    }

   
    // intersection is a 2D segment 
    // sometimes CGAL encodes a single point as a segment
    if (const Segment2* s = std::get_if<Segment2>(&*inter)) {
        const P2& p_src = s->source();
        const P2& p_tgt = s->target();

        // if source == target, treat as single point
        if (approx_equal_p2(p_src, p_tgt)) {
            const P2& ip = p_src;

            FT t_cd = compute_param_t(ip, c2, d2);
            FT t_uv = compute_param_t(ip, u2, v2);

            FT z_cd = interpolate_z(c3, d3, t_cd);
            FT z_uv = interpolate_z(u3, v3, t_uv);

            if (z_cd + eps_height < z_uv) return -1; // (c,d) below (u,v)
            if (z_cd > z_uv + eps_height) return +1; // (c,d) above (u,v)
            return 0;
        }

        // true overlapping 2D segment
        return 0;
    }

    // should not happen if we only call this when segments intersect in 2D
    return 0;
}














// orientation-based height test
// returns true if (c,d) is above (u,v), false if (c,d) is below or coplanar with (u,v)
bool height_test_orientation_based(const P2& c2, const P2& d2,
                                   const P2& u2, const P2& v2,
                                   const P3& c3, const P3& d3,
                                   const P3& u3, const P3& v3)
{
    // determine consistent ordering of (u,v) relative to (c,d)
    auto su = CGAL::orientation(c2, d2, u2); // >0: u is LEFT of c->d, <0: RIGHT
    P2 u2_ordered = u2;
    P2 v2_ordered = v2;
    P3 u3_ordered = u3;
    P3 v3_ordered = v3;

    if (su == CGAL::COLLINEAR) {
        // segments are collinear in 2D -> treat as coplanar / overlapping
        // should not happen since we checked for that before
        return false;
    }
    if (su == CGAL::NEGATIVE) {
        std::swap(u2_ordered, v2_ordered);
        std::swap(u3_ordered, v3_ordered);
    }

    // now u is left of (c,d) and v is right of (c,d)

    const auto o = CGAL::orientation(u3_ordered, v3_ordered, c3, d3);
   
    
    if (o == CGAL::NEGATIVE) {
        // (c,d) above (u,v)
        return true;
    } else if (o == CGAL::POSITIVE) {
        // (c,d) below (u,v)
        return false;
    } else {
        // coplanar
        return false;
    }
}



} // namespace df::reg

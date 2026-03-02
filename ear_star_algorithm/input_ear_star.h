#pragma once

#include "input.h"  // for CGAL types 
#include <vector>
#include <string>




struct LabeledPoint {
    df::vertex_id id;         // 0..n-1 for boundary, -1 for star
    df::P2 p;                  // 2D point GCAL type
    std::string label;         // "0","1",...,"n-1","*"
};



struct InputData_ear_star {
    std::vector<LabeledPoint> points;   // all random points of planar point set A 
    df::Tri2            tri_start;          // triangulation of the input point set, this is the starting point for building all triangultions of this point set via flips and insertions
};


InputData_ear_star make_input_ear_star(int n, double R);
#pragma once

#include <string>
#include <vector>

#include "poset.h"  
#include "input.h" 

namespace regcheck {

struct NodeRegularity {
    int  node_idx = -1;
    int  num_vertices = 0;
    int  num_faces = 0;
    bool is_regular = false;   // only valid if error.empty()
    std::string error;         // if non-empty: something failed
};

// checks regularity for every node triangulation in P1 by calling:
//   wsl M2 -q <script.m2>
// prints a table of non-regular nodes (and errors)
//
// - scale: coordinates are exported as integers round(scale * coord)
// - if points are integer already, set scale=1
// - wsl_cmd: "wsl" by default, we could also use "wsl.exe"
std::vector<NodeRegularity> check_poset_regularity_wsl(
    const df::InputData& D,
    const pst::Poset1& P1,
    double scale = 1e6,
    const std::string& wsl_cmd = "wsl",
    bool print_table = true
);


} // namespace regcheck


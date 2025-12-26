#pragma once
#include "input.h"
#include "poset.h"
#include <vector>


namespace pst2 {

// returns true iff T(t1) <=_2 T(t2): the characteristic section of t1
// lies everywhere below (or equal to) that of t2

bool compare(int t1, int t2, const df::InputData& D, const std::vector<pst::Node>& nodes);

}

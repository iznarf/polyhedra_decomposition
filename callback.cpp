#include "callback.h"
#include "visualization.h"

namespace cb {

VisualizerCallback::VisualizerCallback(const df::OmegaQuad& omega,
                                       const std::string& planar_name,
                                       const std::string& lifted_name)
  : omega_(omega), planar_name_(planar_name), lifted_name_(lifted_name) {}

void VisualizerCallback::on_init(const df::Tri2& tri) {
  viz::register_planar_triangulation(tri, planar_name_);
  viz::register_lifted_triangulation(tri, omega_, lifted_name_);
}

void VisualizerCallback::on_finish() {
  // no-op for now
}

} // namespace cb

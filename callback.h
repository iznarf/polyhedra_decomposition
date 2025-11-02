#pragma once
#include <string>
#include <vector>
#include "input.h"

namespace cb {

// Minimal callback interface — for now we only need init + finish.
struct Callback {
  virtual ~Callback() = default;

  // Called when the initial triangulation is ready (step 0).
  virtual void on_init(const df::Tri2& tri) = 0;

  // Called when you’re done (optional for now).
  virtual void on_finish() {}
};

// A concrete callback that visualizes via Polyscope.
// It registers: (1) planar (z=0) and (2) lifted (z=ω) meshes.
class VisualizerCallback : public Callback {
public:
  VisualizerCallback(const df::OmegaQuad& omega,
                     const std::string& planar_name  = "planar",
                     const std::string& lifted_name  = "lifted");

  void on_init(const df::Tri2& tri) override;
  void on_finish() override; // no-op (kept for symmetry)

private:
  const df::OmegaQuad& omega_;
  std::string planar_name_;
  std::string lifted_name_;
};

} // namespace cb

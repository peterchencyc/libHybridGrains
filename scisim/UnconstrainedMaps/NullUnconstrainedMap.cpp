#include "NullUnconstrainedMap.h"

void NullUnconstrainedMap::flow(const VectorXs &q0, const VectorXs &v0,
                                FlowableSystem &fsys, const unsigned iteration,
                                const scalar &dt, VectorXs &q1, VectorXs &v1) {
  q1 = q0;
  v1 = v0;
}

std::string NullUnconstrainedMap::name() const {
  return "null_unconstrained_map";
}

void NullUnconstrainedMap::serialize(std::ostream &output_stream) const {
  // No state to serialize
}

std::unique_ptr<UnconstrainedMap> NullUnconstrainedMap::clone() const {
  return std::unique_ptr<UnconstrainedMap>{new NullUnconstrainedMap};
}

#ifndef NULL_UNCONSTRAINED_MAP_H
#define NULL_UNCONSTRAINED_MAP_H

#include "UnconstrainedMap.h"

class NullUnconstrainedMap final : public UnconstrainedMap {

public:
  NullUnconstrainedMap() = default;

  virtual ~NullUnconstrainedMap() override = default;

  virtual void flow(const VectorXs &q0, const VectorXs &v0,
                    FlowableSystem &fsys, const unsigned iteration,
                    const scalar &dt, VectorXs &q1, VectorXs &v1) override;

  virtual std::string name() const override;

  virtual void serialize(std::ostream &output_stream) const override;

  virtual std::unique_ptr<UnconstrainedMap> clone() const override;
};

#endif

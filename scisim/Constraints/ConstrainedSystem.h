#ifndef CONSTRAINED_SYSTEM_H
#define CONSTRAINED_SYSTEM_H

#include <memory>

#include "scisim/Math/MathDefines.h"

class Constraint;

class ConstrainedSystem {

public:
  virtual void
  computeActiveSet(const VectorXs &q0, const VectorXs &qp, const VectorXs &v,
                   const bool reduce_bandwidth,
                   std::vector<std::unique_ptr<Constraint>> &active_set) = 0;
  virtual void
  computeImpactBases(const VectorXs &q,
                     const std::vector<std::unique_ptr<Constraint>> &active_set,
                     MatrixXXsc &impact_bases) const = 0;
  virtual void computeContactBases(
      const VectorXs &q, const VectorXs &v,
      const std::vector<std::unique_ptr<Constraint>> &active_set,
      MatrixXXsc &contact_bases) const = 0;

  virtual void clearConstraintCache() = 0;
  virtual void cacheConstraint(const Constraint &constraint,
                               const VectorXs &r) = 0;
  virtual void getCachedConstraintImpulse(const Constraint &constraint,
                                          VectorXs &r) const = 0;
  virtual bool constraintCacheEmpty() const = 0;

  virtual void computeActiveSetNew(
      const VectorXs &q, const VectorXs &v,
      std::vector<std::vector<std::unique_ptr<Constraint>>> &con_table);

protected:
  ConstrainedSystem() = default;
  ConstrainedSystem(const ConstrainedSystem &) = default;
  ConstrainedSystem(ConstrainedSystem &&) noexcept = default;
  ConstrainedSystem &operator=(const ConstrainedSystem &other) = default;
  ConstrainedSystem &operator=(ConstrainedSystem &&other) noexcept = default;
  virtual ~ConstrainedSystem() = 0;
};

#endif

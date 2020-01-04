#ifndef NULL_FRICTION_SOLVER_H
#define NULL_FRICTION_SOLVER_H

#include "FrictionSolver.h"

class NullFrictionSolver final : public FrictionSolver {

public:
  NullFrictionSolver() = default;

  virtual ~NullFrictionSolver() override = default;

  virtual void solve(const unsigned iteration, const scalar &dt,
                     const FlowableSystem &fsys, const SparseMatrixsc &M,
                     const SparseMatrixsc &Minv, const VectorXs &CoR,
                     const VectorXs &mu, const VectorXs &q0, const VectorXs &v0,
                     std::vector<std::unique_ptr<Constraint>> &active_set,
                     const MatrixXXsc &contact_bases,
                     const VectorXs &nrel_extra, const VectorXs &drel_extra,
                     VectorXs &f, VectorXs &alpha, VectorXs &beta,
                     VectorXs &vout, bool &solve_succeeded,
                     scalar &error) override;

  virtual unsigned numFrictionImpulsesPerNormal(
      const unsigned ambient_space_dimensions) const override;

  virtual void serialize(std::ostream &output_stream) const override;

  virtual std::unique_ptr<FrictionSolver> clone() const override;

  virtual std::string name() const override;
};

#endif

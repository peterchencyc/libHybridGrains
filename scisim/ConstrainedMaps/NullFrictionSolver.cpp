#include "NullFrictionSolver.h"

void NullFrictionSolver::solve(
    const unsigned iteration, const scalar &dt, const FlowableSystem &fsys,
    const SparseMatrixsc &M, const SparseMatrixsc &Minv, const VectorXs &CoR,
    const VectorXs &mu, const VectorXs &q0, const VectorXs &v0,
    std::vector<std::unique_ptr<Constraint>> &active_set,
    const MatrixXXsc &contact_bases, const VectorXs &nrel_extra,
    const VectorXs &drel_extra, VectorXs &f, VectorXs &alpha, VectorXs &beta,
    VectorXs &vout, bool &solve_succeeded, scalar &error) {
  f.setZero();
  alpha.setZero();
  beta.setZero();
  vout = v0;
  solve_succeeded = true;
  error = 0.0;
}

unsigned NullFrictionSolver::numFrictionImpulsesPerNormal(
    const unsigned ambient_space_dimensions) const {
  return ambient_space_dimensions - 1;
}

void NullFrictionSolver::serialize(std::ostream &output_stream) const {
  // No state to serialize
}

std::unique_ptr<FrictionSolver> NullFrictionSolver::clone() const {
  return std::unique_ptr<FrictionSolver>{new NullFrictionSolver};
}

std::string NullFrictionSolver::name() const { return "null_friction_solver"; }

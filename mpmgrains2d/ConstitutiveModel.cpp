//
//  ConstitutiveModel.cpp
//  SCISim
//
//  Created by Yonghao Yue on 11/11/16.
//
//
#include "ConstitutiveModel.h"

#include <Eigen/LU>
#include <iostream>
#include <map>

#include "scisim/Math/MathUtilities.h"

namespace ConstitutiveModel {

// explicit elastic be prediction: be* = ben + dt(L ben + ben L^T)
void elasticPrediction(const Matrix22sc &be_n, const Matrix22sc &L,
                       const double &dt, Matrix22sc &be_star) {
  be_star = be_n + dt * (L * be_n + be_n * L.transpose());
}

//*
// Drucker-Pragar
void plasticCorrection(const Matrix22sc &be_star, const scalar kappa,
                       const scalar mu, const scalar alpha,
                       int &just_homogenized, Matrix22sc &be_next) {
  if (just_homogenized != 0) {
    // remove the shear strain:
    const double J = sqrt(be_star.determinant());
    if (J > 1.0) {
      be_next = J * Matrix2s::Identity();
      return;
    }

    // once the material is packed enough, let the material act as a (usual)
    // dens granular continuum;
    just_homogenized = 0;
  }

  Matrix22sc be_star2 = be_star;

  // numerically, J can be slightly bigger than 1.0 even with the
  // re-normalization above.
  const double J = sqrt(be_star2.determinant());
  const double yield_stress =
      std::max<double>(0.0, 0.5 * alpha * kappa * (1.0 + J) * (1.0 - J));
  // const double yield_stress = 0.5 * alpha * kappa * (1.0 + J) * (1.0 - J);

  const Matrix22sc dev_be_star2 =
      be_star2 - 0.5 * be_star2.trace() * Matrix22sc::Identity();
  const Matrix22sc bar_dev_be_star2 = dev_be_star2 / J;

  if (mu * bar_dev_be_star2.norm() > yield_stress) {
    double lambda2 = yield_stress / (mu * bar_dev_be_star2.norm());

    double lambda1_squared =
        be_star2.determinant() - lambda2 * lambda2 * dev_be_star2.determinant();

    // checking:
    double test = 0.5 * (J * yield_stress / mu) * (J * yield_stress / mu);
    if (fabs(test - (-lambda2 * lambda2 * dev_be_star2.determinant())) >
        1.0e-6) {
      std::cout << "error in lambda_1?" << std::endl;
    }

    if (lambda1_squared < 0.0) {
      std::cout << "lambda1^2 is negative..." << std::endl;
    }

    double lambda1 = sqrt(lambda1_squared);
    be_next = lambda1 * Matrix22sc::Identity() + lambda2 * dev_be_star2;

    if (fabs(be_next.determinant() - be_star2.determinant()) > 1.0e-6) {
      std::cout << "determinants should agree" << std::endl;
    }
  } else {
    be_next = be_star2;
  }
}

void computeTau(const Matrix22sc &be, const scalar kappa, const scalar mu,
                const int just_homogenized, Matrix22sc &tau) {
  if (just_homogenized != 0) {
    tau = Matrix22sc::Zero();
    return;
  }

  const scalar J = sqrt(be.determinant());
  const Matrix22sc be_bar = be / pow(be.determinant(), 1.0 / 2.0);
  const Matrix22sc dev_be_bar =
      be_bar - be_bar.trace() * Matrix22sc::Identity() / 2.0;

  tau = kappa * 0.5 * (J + 1.0) * (J - 1.0) * Matrix22sc::Identity() +
        mu * dev_be_bar;
}

void computeStrainFromCauchyStress(const Matrix22sc &sigma, const scalar kappa,
                                   const scalar mu, scalar &J,
                                   Matrix22sc &be_bar) {
  const scalar T = sigma.trace();
  J = 0.5 * (T / kappa + sqrt((T * T) / (kappa * kappa) + 4.0));

  const Matrix22sc dev_sigma =
      sigma - sigma.trace() * Matrix22sc::Identity() / 2.0;
  const Matrix22sc dev_be_bar = J * dev_sigma / mu;

  if (1.0 - dev_be_bar.determinant() < 0.0) {
    std::cout << "ERROR: cannot compute strain from stress... sigma = " << sigma
              << std::endl;
    std::cerr << "ERROR: cannot compute strain from stress... sigma = " << sigma
              << std::endl;
  }
  const scalar t = sqrt(1.0 - dev_be_bar.determinant());
  be_bar = dev_be_bar + t * Matrix22sc::Identity();
}
} // namespace ConstitutiveModel

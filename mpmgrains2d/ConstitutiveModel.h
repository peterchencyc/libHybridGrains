//
//  ConstitutiveModel.h
//  SCISim
//
//  Created by Yonghao Yue on 11/11/16.
//
//

#ifndef ConstitutiveModel_2D_h
#define ConstitutiveModel_2D_h

#include <memory>

#include "scisim/Math/MathDefines.h"

namespace ConstitutiveModel {

void computeBeNext(const Matrix22sc &be_n, const Matrix22sc &L,
                   const double &dt, const scalar kappa, const scalar mu,
                   const scalar alpha, int &just_homogenized,
                   Matrix22sc &be_next);

void elasticPrediction(const Matrix22sc &be_n, const Matrix22sc &L,
                       const double &dt, Matrix22sc &be_star);

void plasticCorrection(const Matrix22sc &be_star, const scalar kappa,
                       const scalar mu, const scalar alpha,
                       int &just_homogenized, Matrix22sc &be_next);

// void computeBeGradient(const Matrix22sc& be_n, const Matrix22sc& L, const
// scalar kappa, const scalar mu, const scalar alpha, const double& dt,
// Eigen::Matrix4d& grad);

void computeTau(const Matrix22sc &be, const scalar kappa, const scalar mu,
                const int just_homogenized, Matrix22sc &tau);

// void computeTauGradient(const Matrix22sc& be, const scalar kappa, const
// scalar mu, Eigen::Matrix4d& grad);

void computeStrainFromCauchyStress(const Matrix22sc &sigma, const scalar kappa,
                                   const scalar mu, scalar &J,
                                   Matrix22sc &be_bar);

} // namespace ConstitutiveModel

#endif /* ConstitutiveModel_h */

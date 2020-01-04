#ifndef MATERIAL_POINTS_2D_H
#define MATERIAL_POINTS_2D_H

#include <Eigen/StdVector>

#include "scisim/Math/MathDefines.h"

#include "BasisFunctions.h"

struct PhysicsGrid;
struct MPMStaticPlane;

// TODO : Move this struct out of here
struct InitialMaterialPoint final {

  InitialMaterialPoint(const Vector2s &position0, const Matrix22sc &be_bar0,
                       const scalar &J0, const Vector2s &velocity0)
      : x(position0), be_bar(be_bar0), J(J0), v(velocity0),
        m(std::numeric_limits<scalar>::signaling_NaN()),
        volume(std::numeric_limits<scalar>::signaling_NaN()) {}

  Vector2s x;
  Matrix22sc be_bar;
  scalar J;
  Vector2s v;
  scalar m;
  scalar volume;

  // the particle extent (half width)
  // this parameter might be replaced by [Matrix22sc extent] in the future for
  // CPDI
  scalar hl;

  bool always_fixed;
  bool kinematically_scripted;
};

struct MaterialPoints final {

  MaterialPoints();

  void clear();

  void serialize(std::ostream &strm) const;
  void deserialize(std::istream &input_stream);

  bool empty() const;

  // Computes an 'easy' (not tight) bounding sphere centered at the average
  // position
  void computeBoundingSphere(scalar &radius, Vector2s &center) const;

  void
  setPoints(const std::vector<InitialMaterialPoint,
                              Eigen::aligned_allocator<InitialMaterialPoint>>
                &points);
  void conservativeResize(const unsigned n);

  scalar computeTotalEnergy() const;
  Vector2s computeTotalMomentum() const;
  scalar computeTotalAngularMomentum() const;

  void computeHyperelasticCauchyStress(const scalar &shear_modulus,
                                       const scalar &bulk_modulus);

  // TODO: Make Ando-style a separate method
  // if viscosity is zero, snow paper style blending (of PIC and FLIP velocities
  // based on valpha) is used; if viscosity is nonzero, Ando 2013 style blending
  // is used. for "zero" viscosity, set both viscosity to zero and valpha to
  // one.
  // void updateVelocity( const std::unique_ptr<BasisFunctions>& in_SF, const
  // PhysicsGrid& in_Grid, const scalar& valpha, const scalar& viscosity, const
  // scalar& dt );
  void updateVelocities(const std::unique_ptr<BasisFunctions> &basis_funcs,
                        const PhysicsGrid &grid, const scalar &alpha,
                        const scalar &dt);
  void updatePositions(const std::unique_ptr<BasisFunctions> &basis_funcs,
                       const PhysicsGrid &grid, const scalar &dt);

  // TODO: Rename kappa and mu
  void
  updateDeformationGradient(const std::unique_ptr<BasisFunctions> &basis_funcs,
                            const PhysicsGrid &grid, const scalar &dt,
                            const scalar &kappa, const scalar &mu);
  // const scalar& alpha, const scalar& beta, // for Drucker Prager
  // (non-associated) model const bool lowDensityStrainFreeModel, const scalar&
  // criticalDensity, const bool simo92, const bool exponentialUpdate, const
  // bool densityUpdate, const scalar& initialParticleVolume )

  void computeVelGrad(const std::unique_ptr<BasisFunctions> &basis_funcs,
                      const PhysicsGrid &grid, const Matrix2Xsc &vel_star);

  void elasticPrediction(
      const scalar &dt,
      const std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>
          &vel_Grad);

  void plasticCorrection(const scalar &kappa, const scalar &mu,
                         const scalar &alpha);

  void
  resolvePlaneCollisions(const std::vector<MPMStaticPlane> &plane_obstacles);

  void addMaterialPoint(const Vector2s &_q0, const Vector2s &_q,
                        const Vector2s &_v, const scalar &_m,
                        const scalar &_volume, const scalar &_hl,
                        const Matrix22sc &_be_bar, const scalar &_J,
                        const int &_just_homogenized);
  void deleteMaterialPoint(const std::vector<unsigned> &idx_to_del,
                           std::vector<int> &old_to_new_idx_ref_map);

  void updateWeights(std::function<scalar(const VectorXs &)> &weight_func);
  void updateMassAndVolumeUsingWeights();

  scalar totalMass() const;

  unsigned npoints;
  // Material point rest locations
  Matrix2Xsc q0;
  // Material point locations
  Matrix2Xsc q;
  // Mateiral point velocities
  Matrix2Xsc v;
  // Material point masses
  VectorXs m;
  // Material point volumes
  VectorXs volume;

  // Initial material point masses
  VectorXs m0;
  // Initial material point volumes
  VectorXs volume0;
  // mass volume weight
  VectorXs weight;

  // the particle extent (half width)
  // this parameter might be replaced by [Matrix9Xsc extent] in the future for
  // CPDI
  VectorXs hl;

  // Volume-preserving left Cauchy-Green strain
  std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>> be_bar;
  // Determinant of the deformation gradient
  VectorXs J;
  // Cauchy stress
  std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>> sigma;
  // Velocity gradient
  std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>> vel_grad;

  VectorXi just_homogenized;

  std::vector<bool> always_fixed;
  std::vector<bool> kinematically_scripted;

  // VectorXs homog_factor;
};

#endif

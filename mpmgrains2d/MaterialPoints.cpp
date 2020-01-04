#include "MaterialPoints.h"

#include "PhysicsGrid.h"
#include "StaticPlane.h"

#include "scisim/Math/MathUtilities.h"
#include "scisim/Utilities.h"

#include "ConstitutiveModel.h"

#include <iostream>

MaterialPoints::MaterialPoints() : npoints(0) {}

void MaterialPoints::clear() {
  npoints = 0;
  q0.resize(2, 0);
  q.resize(2, 0);
  v.resize(2, 0);
  m.resize(0);
  volume.resize(0);
  m0.resize(0);
  volume0.resize(0);
  weight.resize(0);
  hl.resize(0);
  be_bar.clear();
  J.resize(0);
  sigma.clear();
  vel_grad.clear();
  just_homogenized.resize(0);
  always_fixed.clear();
  kinematically_scripted.clear();
}

void MaterialPoints::serialize(std::ostream &strm) const {
  Utilities::serializeBuiltInType(npoints, strm);
  MathUtilities::serialize(q0, strm);
  MathUtilities::serialize(q, strm);
  MathUtilities::serialize(v, strm);
  MathUtilities::serialize(m, strm);
  MathUtilities::serialize(volume, strm);
  MathUtilities::serialize(m0, strm);
  MathUtilities::serialize(volume0, strm);
  MathUtilities::serialize(weight, strm);
  MathUtilities::serialize(hl, strm);
  MathUtilities::serialize(be_bar, strm);
  MathUtilities::serialize(J, strm);
  MathUtilities::serialize(sigma, strm);
  MathUtilities::serialize(vel_grad, strm);
  MathUtilities::serialize(just_homogenized, strm);
  Utilities::serializeVectorBuiltInType(always_fixed, strm);
  Utilities::serializeVectorBuiltInType(kinematically_scripted, strm);
}

void MaterialPoints::deserialize(std::istream &input_stream) {
  npoints = Utilities::deserialize<unsigned>(input_stream);
  q0 = MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  q = MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  v = MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  m = MathUtilities::deserialize<VectorXs>(input_stream);
  volume = MathUtilities::deserialize<VectorXs>(input_stream);
  m0 = MathUtilities::deserialize<VectorXs>(input_stream);
  volume0 = MathUtilities::deserialize<VectorXs>(input_stream);
  weight = MathUtilities::deserialize<VectorXs>(input_stream);
  hl = MathUtilities::deserialize<VectorXs>(input_stream);
  be_bar = MathUtilities::deserialize<
      std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>>(
      input_stream);
  J = MathUtilities::deserialize<VectorXs>(input_stream);
  sigma = MathUtilities::deserialize<
      std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>>(
      input_stream);
  vel_grad = MathUtilities::deserialize<
      std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>>(
      input_stream);
  just_homogenized = MathUtilities::deserialize<VectorXi>(input_stream);
  Utilities::deserializeVectorBuiltInType<bool>(always_fixed, input_stream);
  Utilities::deserializeVectorBuiltInType<bool>(kinematically_scripted,
                                                input_stream);
}

bool MaterialPoints::empty() const { return npoints == 0; }

void MaterialPoints::computeBoundingSphere(scalar &radius,
                                           Vector2s &center) const {
  assert(q.size() != 0);
  center = q.rowwise().mean();
  radius = 0.0;
  for (unsigned pnt_idx = 0; pnt_idx < npoints; pnt_idx++) {
    using std::max;
    radius = max(radius, (q.col(pnt_idx) - center).squaredNorm());
  }
  using std::sqrt;
  radius = sqrt(radius);
}

// TODO: Pull the sanity checks into their own function
void MaterialPoints::setPoints(
    const std::vector<InitialMaterialPoint,
                      Eigen::aligned_allocator<InitialMaterialPoint>> &points) {
  conservativeResize(unsigned(points.size()));
  for (std::vector<InitialMaterialPoint>::size_type pnt_idx = 0;
       pnt_idx < points.size(); pnt_idx++) {
    VectorXs x = points[pnt_idx].x;
    q0.col(long(pnt_idx)) = x;
    q.col(long(pnt_idx)) = x;
    v.col(long(pnt_idx)) = points[pnt_idx].v;
    m0(long(pnt_idx)) = points[pnt_idx].m;
    m(long(pnt_idx)) = points[pnt_idx].m;
    assert(m(long(pnt_idx)) > 0.0);
    volume0(long(pnt_idx)) = points[pnt_idx].volume;
    volume(long(pnt_idx)) = points[pnt_idx].volume;
    assert(volume(long(pnt_idx)) > 0.0);
    weight(long(pnt_idx)) = 1.0;
    hl(long(pnt_idx)) = points[pnt_idx].hl;
    assert(hl(long(pnt_idx)) > 0.0);
    MatrixXXsc be_bar_l = points[pnt_idx].be_bar;
    be_bar[pnt_idx] = be_bar_l;
    assert((be_bar[pnt_idx] - be_bar[pnt_idx].transpose())
               .lpNorm<Eigen::Infinity>() <= 1.0e-9);
    assert(fabs(be_bar[pnt_idx].determinant() - 1.0) <= 1.0e-9);
    J(long(pnt_idx)) = points[pnt_idx].J;
    sigma[long(pnt_idx)].setZero();
    vel_grad[long(pnt_idx)].setZero();
    just_homogenized(pnt_idx) = 0;
    always_fixed[pnt_idx] = points[pnt_idx].always_fixed;
    kinematically_scripted[pnt_idx] = points[pnt_idx].kinematically_scripted;
  }
}

void MaterialPoints::conservativeResize(const unsigned n) {
  npoints = n;
  q0.conservativeResize(2, n);
  q.conservativeResize(2, n);
  v.conservativeResize(2, n);
  m.conservativeResize(n);
  volume.conservativeResize(n);
  m0.conservativeResize(n);
  volume0.conservativeResize(n);
  weight.conservativeResize(n);
  hl.conservativeResize(n);
  be_bar.resize(n);
  J.conservativeResize(n);
  sigma.resize(n);
  vel_grad.resize(n);
  just_homogenized.conservativeResize(n);
  always_fixed.resize(n);
  kinematically_scripted.resize(n);
}

scalar MaterialPoints::computeTotalEnergy() const {
  static bool warning_printed = false;
  if (!warning_printed) {
    std::cerr << "Warning, material point energy computation not coded up yet."
              << std::endl;
    warning_printed = true;
  }
  return 0.0;
}

Vector2s MaterialPoints::computeTotalMomentum() const {
  assert(npoints == m.size());
  assert(npoints == v.cols());
  Vector2s total_momentum = Vector2s::Zero();
  for (unsigned pnt_idx = 0; pnt_idx < npoints; pnt_idx++) {
    total_momentum += m(pnt_idx) * v.col(pnt_idx);
  }
  return total_momentum;
}

scalar MaterialPoints::computeTotalAngularMomentum() const {
  assert(npoints == m.size());
  assert(npoints == q.cols());
  assert(npoints == v.cols());
  scalar total_angular_momentum = 0.0;
  for (unsigned pnt_idx = 0; pnt_idx < npoints; pnt_idx++) {
    total_angular_momentum +=
        m(pnt_idx) * MathUtilities::cross(q.col(pnt_idx), v.col(pnt_idx));
  }
  return total_angular_momentum;
}

void MaterialPoints::computeHyperelasticCauchyStress(
    const scalar &shear_modulus, const scalar &bulk_modulus) {
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < npoints; pnt_idx++) {
    assert((be_bar[pnt_idx] - be_bar[pnt_idx].transpose())
               .lpNorm<Eigen::Infinity>() <= 1.0e-9);
    assert(fabs(be_bar[pnt_idx].determinant() - 1.0) <= 1.0e-9);

    Matrix22sc tau{Matrix22sc::Zero()};
    if (J(pnt_idx) <= 1.0) {
      ConstitutiveModel::computeTau(be_bar[pnt_idx] * J(pnt_idx), bulk_modulus,
                                    shear_modulus, just_homogenized(pnt_idx),
                                    tau);
    }
    assert(J(pnt_idx) != 0.0);
    sigma[pnt_idx] = tau / J(pnt_idx);
    assert((sigma[pnt_idx] - sigma[pnt_idx].transpose())
               .lpNorm<Eigen::Infinity>() <= 1.0e-9);
  }
}

void MaterialPoints::updateVelocities(
    const std::unique_ptr<BasisFunctions> &basis_funcs, const PhysicsGrid &grid,
    const scalar &alpha, const scalar &dt) {

#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < npoints; pnt_idx++) {
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, *this, grid)};
    assert((stencil.first < stencil.second).all());

    Vector2s vpic{Vector2s::Zero()};
    Vector2s accel{Vector2s::Zero()};

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(grid.nodeIndicesValid(x_idx, y_idx));

        const scalar wght{
            basis_funcs->weight(pnt_idx, *this, {x_idx, y_idx}, grid)};
#ifndef NDEBUG
        weight_sum += wght;
#endif

        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};

        vpic += wght * grid.velocity.col(flat_node_idx);
        accel += wght * grid.acceleration.col(flat_node_idx);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);

    const Vector2s vflip{v.col(pnt_idx) + accel * dt};

    if (!kinematically_scripted[pnt_idx])
      v.col(pnt_idx) = (1.0 - alpha) * vpic + alpha * vflip;
  }
}

void MaterialPoints::updatePositions(
    const std::unique_ptr<BasisFunctions> &basis_funcs, const PhysicsGrid &grid,
    const scalar &dt) {
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < npoints; pnt_idx++) {
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, *this, grid)};
    assert((stencil.first < stencil.second).all());

    Vector2s gvel{0.0, 0.0};

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(grid.nodeIndicesValid(x_idx, y_idx));

        // Switching to the new interface, for the compatibility with uGIMP
        // bases. const scalar weight{ basis_funcs->weight( px, { x_idx, y_idx
        // }, grid ) };
        const scalar wght{
            basis_funcs->weight(pnt_idx, *this, {x_idx, y_idx}, grid)};
#ifndef NDEBUG
        weight_sum += wght;
#endif

        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};

        gvel += wght * grid.velocity.col(flat_node_idx);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);

    if (!kinematically_scripted[pnt_idx])
      q.col(pnt_idx) += dt * gvel;
    else
      q.col(pnt_idx) += dt * v.col(pnt_idx);
  }
}

static Matrix22sc normalizedDeterminant(const Matrix22sc &A) {
  assert(A.determinant() != 0.0);
  using std::pow;
  return A / pow(A.determinant(), 1.0 / 2.0);
}

void MaterialPoints::computeVelGrad(
    const std::unique_ptr<BasisFunctions> &basis_funcs, const PhysicsGrid &grid,
    const Matrix2Xsc &vel_star) {
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < npoints; ++pnt_idx) {

    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, *this, grid)};

    // Compute the gradient of the velocity
    Matrix22sc grad_v{Matrix22sc::Zero()};

    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        const Vector2s grad_weight{
            basis_funcs->weightGrad(pnt_idx, *this, {x_idx, y_idx}, grid)};
        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};

        grad_v += vel_star.col(flat_node_idx) * grad_weight.transpose();
      }
    }

    vel_grad[pnt_idx] = grad_v;
  }
}

void MaterialPoints::elasticPrediction(
    const scalar &in_dt,
    const std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>
        &vel_Grad) {
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < npoints; ++pnt_idx) {
    // Get a local mutable copy of be_bar
    const Matrix22sc be_bar_map{be_bar[pnt_idx]};
    assert((be_bar_map - be_bar_map.transpose()).lpNorm<Eigen::Infinity>() <=
           1.0e-6);
    assert(fabs(be_bar_map.determinant() - 1.0) <= 1.0e-6);

    // Eigen::Map<Matrix22sr> L{ m_velGrad.col(pnt_idx).data() };
    const Matrix22sc L{vel_Grad[pnt_idx]};
    Matrix22sc be_star;
    ConstitutiveModel::elasticPrediction(be_bar_map * J(pnt_idx), L, in_dt,
                                         be_star);

    const scalar J_star{sqrt(be_star.determinant())};
    assert(J_star > 0.0);

    J(pnt_idx) = J_star;
    be_bar[pnt_idx] = normalizedDeterminant(be_star);
  }
}

void MaterialPoints::plasticCorrection(const scalar &kappa, const scalar &mu,
                                       const scalar &alpha) {
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < npoints; ++pnt_idx) {
    const Matrix22sc be_bar_map{be_bar[pnt_idx]};

    Matrix22sc be_next;
    ConstitutiveModel::plasticCorrection(be_bar_map * J(pnt_idx), kappa, mu,
                                         alpha, just_homogenized(pnt_idx),
                                         be_next);

    // Finalize the deformation gradient
    J(pnt_idx) = sqrt(be_next.determinant());

    // Finalize the normalized strain
    be_bar[pnt_idx] = normalizedDeterminant(be_next);

    assert(fabs(be_bar[pnt_idx].determinant() - 1.0) <= 1.0e-6);
  }
}

void MaterialPoints::resolvePlaneCollisions(
    const std::vector<MPMStaticPlane> &plane_obstacles) {
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < npoints; pnt_idx++) {
    for (const MPMStaticPlane &plane_obstacle : plane_obstacles) {
      // if below the lower bound of the plane
      if (q.col(pnt_idx)(1) < plane_obstacle.lower_bound) {
        continue;
      }

      // Compute the signed distance to the plane
      assert(fabs(plane_obstacle.n.norm() - 1.0) <= 1.0e-6);
      const scalar d{plane_obstacle.n.dot(q.col(pnt_idx) - plane_obstacle.x)};

      // If the objects are not colliding, move on to the next particle
      if (d > 0.0) {
        continue;
      }

      // Move particle out of obstacle if stuck there; position based collision
      // response
      q.col(pnt_idx) += -d * plane_obstacle.n;

      switch (plane_obstacle.boundary_behavior) {
      case MPMStaticPlane::BoundaryBehavior::SLIDING: {
        // Compute the normal relative velocity
        assert(fabs(plane_obstacle.n.norm() - 1.0) <= 1.0e-6);
        const scalar vnormal{plane_obstacle.n.dot(v.col(pnt_idx))};

        // If the bodies are separating, no response
        if (vnormal > 0.0) {
          continue;
        }

        // Kill off the normal relative velocity
        const Vector2s no_normal_vel{v.col(pnt_idx) -
                                     vnormal * plane_obstacle.n};

        // Update the final velocity
        v.col(pnt_idx) = no_normal_vel;

        break;
      }
      case MPMStaticPlane::BoundaryBehavior::STICKING: {
        // Fully sticking
        v.col(pnt_idx).setZero();

        break;
      }
      }
    }
  }
}

void MaterialPoints::updateWeights(
    std::function<scalar(const VectorXs &)> &weight_func) {
  for (unsigned idx = 0; idx < npoints; idx++) {
    weight(idx) = weight_func(q.col(idx));
  }
}

void MaterialPoints::updateMassAndVolumeUsingWeights() {
  for (unsigned idx = 0; idx < npoints; idx++) {
    m(idx) = m0(idx) * weight(idx);
    volume(idx) = volume0(idx) * weight(idx);
  }
}

void MaterialPoints::addMaterialPoint(const Vector2s &_q0, const Vector2s &_q,
                                      const Vector2s &_v, const scalar &_m,
                                      const scalar &_volume, const scalar &_hl,
                                      const Matrix22sc &_be_bar,
                                      const scalar &_J,
                                      const int &_just_homogenized) {
  const int idx = npoints++;
  conservativeResize(npoints);

  q0.col(idx) = _q0;
  q.col(idx) = _q;
  v.col(idx) = _v;
  m(idx) = _m;
  volume(idx) = _volume;
  m0(idx) = _m;
  volume0(idx) = _volume;
  weight(idx) = 1.0;

  hl(idx) = _hl;

  be_bar[idx] = _be_bar;
  J(idx) = _J;
  sigma[idx] = Matrix22sc::Zero();
  vel_grad[idx] = Matrix22sc::Zero();
  just_homogenized[idx] = _just_homogenized;
  always_fixed[idx] = false;
  kinematically_scripted[idx] = false;
  // homog_factor[ idx ] = 1.0;
}

void MaterialPoints::deleteMaterialPoint(
    const std::vector<unsigned> &idx_to_del,
    std::vector<int> &old_to_new_idx_ref_map) {
  old_to_new_idx_ref_map.clear();
  unsigned pt = 0;
  unsigned dest = 0;

  int nPointsBeforeDelete = npoints;

  for (unsigned idx = 0; idx < static_cast<unsigned>(nPointsBeforeDelete);
       idx++) {
    bool erase = false;

    if (pt < idx_to_del.size()) {
      if (idx_to_del[pt] == idx) {
        erase = true;
        pt++;
      }
    }

    if (!erase) {
      if (idx != dest) {
        // physically copy src(idx) to dest;
        q0.block(0, dest, 2, 1) = q0.block(0, idx, 2, 1);
        q.block(0, dest, 2, 1) = q.block(0, idx, 2, 1);
        v.block(0, dest, 2, 1) = v.block(0, idx, 2, 1);
        m.segment(dest, 1) = m.segment(idx, 1);
        volume.segment(dest, 1) = volume.segment(idx, 1);
        m0.segment(dest, 1) = m0.segment(idx, 1);
        volume0.segment(dest, 1) = volume0.segment(idx, 1);
        weight.segment(dest, 1) = weight.segment(idx, 1);

        hl.segment(dest, 1) = hl.segment(idx, 1);

        be_bar[dest] = be_bar[idx];
        J.segment(dest, 1) = J.segment(idx, 1);
        sigma[dest] = sigma[idx];
        vel_grad[dest] = vel_grad[idx];
        just_homogenized[dest] = just_homogenized[idx];
        always_fixed[dest] = always_fixed[idx];
        kinematically_scripted[dest] = kinematically_scripted[idx];
      }
      old_to_new_idx_ref_map.push_back(dest);
      dest++;
    } else {
      old_to_new_idx_ref_map.push_back(-1);
      --npoints;
    }
  }

  conservativeResize(npoints);
}

scalar MaterialPoints::totalMass() const { return m.sum(); }

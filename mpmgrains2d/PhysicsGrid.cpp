#include "PhysicsGrid.h"

#include "BasisFunctions.h"
#include "MaterialPoints.h"
#include "StaticPlane.h"

#include "scisim/Math/MathUtilities.h"
#include "scisim/Utilities.h"

#ifdef OPENMP_ENABLED
#include <omp.h>
#endif

void PhysicsGrid::setDimensions(const Array2s &grid_min,
                                const Array2s &grid_max, const scalar &cell_w) {
  // Center of the region
  const Array2s center = 0.5 * (grid_min + grid_max);

  // Width of the region
  Array2s width = grid_max - grid_min;

  // Number of grid cells along each dimension
  using std::ceil;
  cell_count = (width / cell_w)
                   .unaryExpr([](const scalar &s) { return ceil(s); })
                   .cast<unsigned>();

  // New region width
  width = cell_count.cast<scalar>() * cell_w;

  // New region
  min = center - 0.5 * width;
  max = center + 0.5 * width;
  assert((min <= grid_min).all());
  assert((max >= grid_max).all());

  cell_width = cell_w;

  // start + h * dimensions == end
  assert(
      (max - (min + cell_width * cell_count.cast<scalar>())).abs().maxCoeff() <=
      1.0e-6);

  rasterized_mass.resize(numGridPoints());
  rasterized_momentum.resize(2, numGridPoints());
  rasterized_momentum_new.resize(2, numGridPoints());
  force.resize(2, numGridPoints());
  velocity.resize(2, numGridPoints());
  acceleration.resize(2, numGridPoints());
  coupling_force.setZero(2, numGridPoints());

  homog_factor.setOnes(numGridPoints());

  hybridized.resize(numGridPoints());
}

unsigned PhysicsGrid::numGridPoints() const {
  return (cell_count.x() + 1) * (cell_count.y() + 1);
}

unsigned PhysicsGrid::numGridCells() const {
  return cell_count.x() * cell_count.y();
}

Vector2u PhysicsGrid::computeNodeTripletIndex(const unsigned flat_idx) const {
  const Vector2u point(flat_idx % (cell_count.x() + 1),
                       flat_idx / (cell_count.x() + 1));
  return point;
}

Vector2s PhysicsGrid::computeGridPointLocation(const unsigned flat_idx) const {
  const Vector2u grid_idx = computeNodeTripletIndex(flat_idx);
  return min.matrix() + cell_width * grid_idx.cast<scalar>();
}

unsigned PhysicsGrid::flatNodeIndex(const unsigned xidx,
                                    const unsigned yidx) const {
  assert(nodeIndicesValid(xidx, yidx));
  return yidx * (cell_count.x() + 1) + xidx;
}

// #ifndef NDEBUG
bool PhysicsGrid::nodeIndicesValid(const unsigned xidx,
                                   const unsigned yidx) const {
  if (xidx > cell_count.x()) {
    return false;
  }
  if (yidx > cell_count.y()) {
    return false;
  }
  return true;
}
// #endif

void PhysicsGrid::clear() {
  rasterized_mass.resize(0);
  rasterized_momentum.resize(2, 0);
  rasterized_momentum_new.resize(2, 0);
  force.resize(2, 0);
  velocity.resize(2, 0);
  acceleration.resize(2, 0);
  coupling_force.resize(2, 0);
  homog_factor.resize(0);
  hybridized.resize(0);
}

void PhysicsGrid::clearRasterizedData() {
  rasterized_mass.setZero();
  rasterized_momentum.setZero();
  rasterized_momentum_new.setZero();
  force.setZero();
  velocity.setZero();
  acceleration.setZero();
}

void PhysicsGrid::clearCouplingForce() { coupling_force.setZero(); }

void PhysicsGrid::serialize(std::ostream &strm) const {
  MathUtilities::serialize(min, strm);
  MathUtilities::serialize(max, strm);
  Utilities::serializeBuiltInType(cell_width, strm);
  MathUtilities::serialize(cell_count, strm);
  MathUtilities::serialize(rasterized_mass, strm);
  MathUtilities::serialize(rasterized_momentum, strm);
  MathUtilities::serialize(rasterized_momentum_new, strm);
  MathUtilities::serialize(force, strm);
  MathUtilities::serialize(velocity, strm);
  MathUtilities::serialize(acceleration, strm);
  MathUtilities::serialize(homog_factor, strm);
  MathUtilities::serialize(hybridized, strm);
}

void PhysicsGrid::deserialize(std::istream &input_stream) {
  min = MathUtilities::deserialize<Array2s>(input_stream);
  max = MathUtilities::deserialize<Array2s>(input_stream);
  cell_width = Utilities::deserialize<scalar>(input_stream);
  cell_count = MathUtilities::deserialize<Array2u>(input_stream);
  rasterized_mass = MathUtilities::deserialize<VectorXs>(input_stream);
  rasterized_momentum = MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  rasterized_momentum_new =
      MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  force = MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  velocity = MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  acceleration = MathUtilities::deserialize<Matrix2Xsc>(input_stream);
  homog_factor = MathUtilities::deserialize<VectorXs>(input_stream);
  hybridized = MathUtilities::deserialize<VectorXi>(input_stream);
}

Vector2s PhysicsGrid::computeTotalMomentum() const {
  // TODO: Colwise sum?
  Vector2s total_mom = Vector2s::Zero();
  for (int node_idx = 0; node_idx < rasterized_momentum.cols(); node_idx++) {
    total_mom += rasterized_momentum.col(node_idx);
  }
  return total_mom;
}

scalar PhysicsGrid::computeTotalAngularMomentum() const {
  scalar ang_mom = 0.0;
  for (int node_idx = 0; node_idx < rasterized_momentum.cols(); node_idx++) {
    ang_mom +=
        MathUtilities::cross(computeGridPointLocation(unsigned(node_idx)),
                             rasterized_momentum.col(node_idx));
  }
  return ang_mom;
}

void PhysicsGrid::rasterizePointMasses(
    const MaterialPoints &points,
    const std::unique_ptr<BasisFunctions> &basis_funcs) {
#ifndef OPENMP_ENABLED
  // Precondition: grid mass must be zero
  assert((rasterized_mass.array() == 0.0).all());
  // Precondition: point masses must be positive
  assert((points.m.array() > 0.0).all());

  // For each material point
  for (unsigned pnt_idx = 0; pnt_idx < points.npoints; pnt_idx++) {
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, points, *this)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        const scalar weight{
            basis_funcs->weight(pnt_idx, points, {x_idx, y_idx}, *this)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const unsigned flat_node_idx{flatNodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < rasterized_mass.size());
        rasterized_mass(flat_node_idx) += weight * points.m(pnt_idx);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

  // Postcondition: grid masses are non-negative
  assert((rasterized_mass.array() >= -2.0e-19).all());
  // Postcondition: total grid mass is equal to total point mass
  assert(fabs(points.m.sum() - rasterized_mass.sum()) <= 1.0e-6);
#else
  // Precondition: grid mass must be zero
  assert((rasterized_mass.array() == 0.0).all());
  // Precondition: point masses must be positive
  assert((points.m.array() > 0.0).all());

  // Allocate per-thread space to rasterize masses into
  static std::vector<VectorXs> per_thread_masses;
  if (per_thread_masses.empty()) {
    per_thread_masses.resize(omp_get_max_threads());
  }
  assert(int(per_thread_masses.size()) == omp_get_max_threads());

  // Zero out the storage
#pragma omp parallel for
  for (std::vector<VectorXs>::size_type idx = 0; idx < per_thread_masses.size();
       idx++) {
    per_thread_masses[idx].setZero(rasterized_mass.size());
  }

// For each material point
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < points.npoints; pnt_idx++) {
    const int tid{omp_get_thread_num()};
    assert(tid >= 0);
    assert(tid < int(per_thread_masses.size()));
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, points, *this)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        const scalar weight{
            basis_funcs->weight(pnt_idx, points, {x_idx, y_idx}, *this)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const unsigned flat_node_idx{flatNodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < rasterized_mass.size());
        per_thread_masses[tid](flat_node_idx) += weight * points.m(pnt_idx);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

#pragma omp parallel for
  for (unsigned node_idx = 0; node_idx < rasterized_mass.size(); node_idx++) {
    for (std::vector<VectorXs>::size_type idx = 0;
         idx < per_thread_masses.size(); idx++) {
      rasterized_mass(node_idx) += per_thread_masses[idx](node_idx);
    }
  }

  // Postcondition: grid masses are non-negative
  assert((rasterized_mass.array() >= -2.0e-19).all());
  // Postcondition: total grid mass is equal to total point mass
  assert(fabs(points.m.sum() - rasterized_mass.sum()) <= 1.0e-6);
#endif
}

void PhysicsGrid::rasterizePointMomentum(
    const MaterialPoints &points,
    const std::unique_ptr<BasisFunctions> &basis_funcs) {
#ifndef OPENMP_ENABLED
  // Precondition: grid momenta must be zero
  assert((rasterized_momentum.array() == 0.0).all());
  // Precondition: grid masses must be non zero
  assert((rasterized_mass.array() >= -2.0e-19).all());
  // Precondition: point masses must be positive
  assert((points.m.array() > 0.0).all());

  // For each material point
  for (unsigned pnt_idx = 0; pnt_idx < points.npoints; pnt_idx++) {
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, points, *this)};
    assert((stencil.first < stencil.second).all());

    const Vector2s pnt_mom{points.m(pnt_idx) * points.v.col(pnt_idx)};

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        const scalar weight{
            basis_funcs->weight(pnt_idx, points, {x_idx, y_idx}, *this)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const unsigned flat_node_idx{flatNodeIndex(x_idx, y_idx)};
        rasterized_momentum.col(flat_node_idx) += weight * pnt_mom;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

// Postcondition: grid momentum is equal to point momentum
#ifndef NDEBUG
  {
    const Vector2s p0{points.computeTotalMomentum()};
    const Vector2s p1{computeTotalMomentum()};
    assert((p0 - p1).lpNorm<Eigen::Infinity>() <= 1.0e-6);
  }
#endif
// Postcondition: grid angular momentum is equal to point angular momentum
#ifndef NDEBUG
  {
    const scalar L0{points.computeTotalAngularMomentum()};
    const scalar L1{computeTotalAngularMomentum()};
    assert(fabs(L1 - L0) <= 1.0e-6);
  }
#endif
#else
  // Precondition: grid momenta must be zero
  assert((rasterized_momentum.array() == 0.0).all());
  // Precondition: grid masses must be non zero
  assert((rasterized_mass.array() >= -2.0e-19).all());
  // Precondition: point masses must be positive
  assert((points.m.array() > 0.0).all());

  // Allocate per-thread space to rasterize momentum into
  static std::vector<Matrix2Xsc> per_thread_momentum;
  if (per_thread_momentum.empty()) {
    per_thread_momentum.resize(omp_get_max_threads());
  }
  assert(int(per_thread_momentum.size()) == omp_get_max_threads());

  // Zero out the storage
#pragma omp parallel for
  for (std::vector<Matrix2Xsc>::size_type idx = 0;
       idx < per_thread_momentum.size(); idx++) {
    per_thread_momentum[idx].setZero(rasterized_momentum.rows(),
                                     rasterized_momentum.cols());
  }

// For each material point
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < points.npoints; pnt_idx++) {
    const int tid{omp_get_thread_num()};
    assert(tid >= 0);
    assert(tid < int(per_thread_momentum.size()));
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, points, *this)};
    assert((stencil.first < stencil.second).all());

    const Vector2s pnt_mom{points.m(pnt_idx) * points.v.col(pnt_idx)};

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        const scalar weight{
            basis_funcs->weight(pnt_idx, points, {x_idx, y_idx}, *this)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const unsigned flat_node_idx{flatNodeIndex(x_idx, y_idx)};
        per_thread_momentum[tid].col(flat_node_idx) += weight * pnt_mom;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

  assert(rasterized_mass.size() == rasterized_momentum.cols());
#pragma omp parallel for
  for (unsigned node_idx = 0; node_idx < rasterized_mass.size(); node_idx++) {
    for (std::vector<Matrix2Xsc>::size_type idx = 0;
         idx < per_thread_momentum.size(); idx++) {
      rasterized_momentum.col(node_idx) +=
          per_thread_momentum[idx].col(node_idx);
    }
  }

  // Postcondition: grid momentum is equal to point momentum
#ifndef NDEBUG
  {
    const Vector2s p0{points.computeTotalMomentum()};
    const Vector2s p1{computeTotalMomentum()};
    assert((p0 - p1).lpNorm<Eigen::Infinity>() <= 1.0e-6);
  }
#endif
  // Postcondition: grid angular momentum is equal to point angular momentum
#ifndef NDEBUG
  {
    const scalar L0{points.computeTotalAngularMomentum()};
    const scalar L1{computeTotalAngularMomentum()};
    assert(fabs(L1 - L0) <= 1.0e-6);
  }
#endif
#endif
}

void PhysicsGrid::computeForces(
    const MaterialPoints &points,
    const std::unique_ptr<BasisFunctions> &basis_funcs,
    const Vector2s &near_earth_gravity) {
#ifndef OPENMP_ENABLED
  assert((force.array() == 0.0).all());

  for (unsigned pnt_idx = 0; pnt_idx < points.npoints; pnt_idx++) {
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, points, *this)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif

    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        const scalar weight{
            basis_funcs->weight(pnt_idx, points, {x_idx, y_idx}, *this)};
#ifndef NDEBUG
        weight_sum += weight;
#endif

        if (weight == 0.0) {
          continue;
        }

        const Vector2s wg =
            basis_funcs->weightGrad(pnt_idx, points, {x_idx, y_idx}, *this);

        // f = - V J sigma grad_w
        const Vector2s f_internal = -points.volume(pnt_idx) *
                                    points.J(pnt_idx) *
                                    (points.sigma[pnt_idx] * wg);

        const unsigned flat_node_idx{flatNodeIndex(x_idx, y_idx)};

        // Internal forces
        force.col(flat_node_idx) += f_internal;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

  // Near earth gravity
  assert((near_earth_gravity.array() == near_earth_gravity.array()).all());
  for (unsigned i = 0; i < rasterized_mass.size(); ++i) {
    if (hybridized(i) > 0)
      force.col(i) += 0.5 * rasterized_mass(i) * near_earth_gravity;
    else
      force.col(i) += rasterized_mass(i) * near_earth_gravity;
  }

  // Coupling Force
  for (unsigned i = 0; i < rasterized_mass.size(); ++i) {
    force.col(i) += coupling_force.col(i);
  }
#else
  assert((force.array() == 0.0).all());

  // Allocate per-thread space to compute forces on
  static std::vector<Matrix2Xsc> per_thread_forces;
  if (per_thread_forces.empty()) {
    per_thread_forces.resize(omp_get_max_threads());
  }
  assert(int(per_thread_forces.size()) == omp_get_max_threads());

  // Zero out the storage
#pragma omp parallel for
  for (std::vector<Matrix2Xsc>::size_type idx = 0;
       idx < per_thread_forces.size(); idx++) {
    per_thread_forces[idx].setZero(rasterized_momentum.rows(),
                                   rasterized_momentum.cols());
  }

  // For each material point
#pragma omp parallel for
  for (unsigned pnt_idx = 0; pnt_idx < points.npoints; pnt_idx++) {
    const int tid{omp_get_thread_num()};
    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(pnt_idx, points, *this)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif

    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        const scalar weight{
            basis_funcs->weight(pnt_idx, points, {x_idx, y_idx}, *this)};
#ifndef NDEBUG
        weight_sum += weight;
#endif

        if (weight == 0.0) {
          continue;
        }
        const Vector2s wg =
            basis_funcs->weightGrad(pnt_idx, points, {x_idx, y_idx}, *this);

        // f = - V J sigma grad_w
        const Vector2s f_internal = -points.volume(pnt_idx) *
                                    points.J(pnt_idx) *
                                    (points.sigma[pnt_idx] * wg);

        const unsigned flat_node_idx{flatNodeIndex(x_idx, y_idx)};

        // Internal forces
        per_thread_forces[tid].col(flat_node_idx) += f_internal;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

  assert(rasterized_mass.size() == force.cols());
#pragma omp parallel for
  for (unsigned node_idx = 0; node_idx < rasterized_mass.size(); node_idx++) {
    for (std::vector<Matrix2Xsc>::size_type idx = 0;
         idx < per_thread_forces.size(); idx++) {
      force.col(node_idx) += per_thread_forces[idx].col(node_idx);
    }
  }

  // Internal forces shouldn't introduce momentum
  assert(force.rowwise().sum().lpNorm<Eigen::Infinity>() <= 1.0e-6);

  // Near earth gravity
  assert((near_earth_gravity.array() == near_earth_gravity.array()).all());
#pragma omp parallel for
  for (unsigned i = 0; i < rasterized_mass.size(); ++i) {
    if (hybridized(i) > 0)
      force.col(i) += 0.5 * rasterized_mass(i) * near_earth_gravity;
    else
      force.col(i) += rasterized_mass(i) * near_earth_gravity;
  }

// Coupling Force
#pragma omp parallel for
  for (unsigned i = 0; i < rasterized_mass.size(); ++i) {
    force.col(i) += coupling_force.col(i);
  }
#endif
}

void PhysicsGrid::updateMomentum(const scalar &dt) {
  assert(rasterized_mass.size() == force.cols());
  assert(rasterized_momentum.cols() == force.cols());
  assert(rasterized_momentum_new.cols() == force.cols());

#pragma omp parallel for
  for (int node_idx = 0; node_idx < force.cols(); node_idx++) {
    assert(rasterized_mass(node_idx) >= -2.0e-19);
#ifndef NDEBUG
    if (rasterized_mass(node_idx) == 0.0) {
      assert((rasterized_momentum.col(node_idx).array() == 0.0).all());
      assert((force.col(node_idx).array() == 0.0).all());
    }
#endif
    // TODO: Better to skip this if? Then we don't even need a loop.
    if (rasterized_mass(node_idx) > 0.0) {
      rasterized_momentum_new.col(node_idx) =
          rasterized_momentum.col(node_idx) + dt * force.col(node_idx);
    }
  }
  // TODO: rasterized_momentum_new = rasterized_momentum + dt * force
}

void PhysicsGrid::lumpedMassVelocityAndAccelerationUpdate(const scalar &dt) {
  assert(rasterized_mass.size() == force.cols());
  assert(rasterized_momentum.cols() == force.cols());
  assert(rasterized_momentum_new.cols() == force.cols());
  assert(velocity.cols() == force.cols());
  assert(acceleration.cols() == force.cols());

#pragma omp parallel for
  for (int node_idx = 0; node_idx < force.cols(); node_idx++) {
    assert(rasterized_mass(node_idx) >= -2.0e-19);
    if (rasterized_mass(node_idx) > 0.0) {
      velocity.col(node_idx) =
          rasterized_momentum_new.col(node_idx) / rasterized_mass(node_idx);
      acceleration.col(node_idx) = (rasterized_momentum_new.col(node_idx) -
                                    rasterized_momentum.col(node_idx)) /
                                   (dt * rasterized_mass(node_idx));
    } else {
      velocity.col(node_idx).setZero();
      acceleration.col(node_idx).setZero();
    }
  }
}

void PhysicsGrid::resolvePlaneCollisions(
    const std::vector<MPMStaticPlane> &plane_obstacles) {
#pragma omp parallel for
  for (int node_idx = 0; node_idx < rasterized_momentum_new.cols();
       ++node_idx) {
    // Compute the world space position of this grid point
    const Vector2s px{computeGridPointLocation(node_idx)};

    for (const MPMStaticPlane &plane_obstacle : plane_obstacles) {
      // if below the lower bound of the plane
      if (px(1) < plane_obstacle.lower_bound) {
        continue;
      }

      // Compute the signed distance to the plane
      assert(fabs(plane_obstacle.n.norm() - 1.0) <= 1.0e-6);
      const scalar d{plane_obstacle.n.dot(px - plane_obstacle.x)};

      // If the objects are not colliding, move on to the next particle
      if (d > 0.0) {
        continue;
      }

      switch (plane_obstacle.boundary_behavior) {
      case MPMStaticPlane::BoundaryBehavior::SLIDING: {
        // Skip nodes with no mass
        if (rasterized_mass(node_idx) <= 0.0) {
          continue;
        }

        // Compute the velocity at this grid node
        // TODO: Because momentum is always parallel to velocity, for a material
        // point, can probably skip this step
        Vector2s vel{rasterized_momentum_new.col(node_idx) /
                     rasterized_mass(node_idx)};

        // Compute the normal relative velocity
        assert(fabs(plane_obstacle.n.norm() - 1.0) <= 1.0e-6);
        const scalar vnormal{plane_obstacle.n.dot(vel)};

        // If the bodies are separating, no response
        if (vnormal > 0.0) {
          continue;
        }

        // Kill off the normal relative velocity
        vel -= vnormal * plane_obstacle.n;

        // Update the final momentum
        rasterized_momentum_new.col(node_idx) = rasterized_mass(node_idx) * vel;

        break;
      }
      case MPMStaticPlane::BoundaryBehavior::STICKING: {
        // Fully sticking
        rasterized_momentum_new.col(node_idx).setZero();

        break;
      }
      }
    }
  }
}

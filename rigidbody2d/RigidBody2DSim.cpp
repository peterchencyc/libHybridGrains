// RigidBody2DSim.cpp
//
// Breannan Smith
// Last updated: 12/07/2015

#include "RigidBody2DSim.h"

#include "scisim/CollisionDetection/CollisionDetectionUtilities.h"
#include "scisim/ConstrainedMaps/ImpactFrictionMap.h"
#include "scisim/ConstrainedMaps/ImpactMaps/ImpactMap.h"
#include "scisim/Constraints/ContactGraphTools.h"
#include "scisim/HDF5File.h"
#include "scisim/Math/MathUtilities.h"
#include "scisim/Math/Rational.h"
#include "scisim/UnconstrainedMaps/UnconstrainedMap.h"
#include "scisim/Utilities.h"

#include "AnnulusGeometry.h"
#include "BodyBodyConstraint.h"
#include "BoxBoxTools.h"
#include "BoxGeometry.h"
#include "CircleBoxTools.h"
#include "CircleCircleConstraint.h"
#include "CircleGeometry.h"
#include "KinematicKickCircleCircleConstraint.h"
#include "KinematicObjectCircleConstraint.h"
#include "PythonScripting.h"
#include "SpatialGrid.h"
#include "StateOutput.h"
#include "StaticDrumCircleConstraint.h"
#include "StaticPlaneBodyConstraint.h"
#include "StaticPlaneCircleConstraint.h"
#include "TeleportedCircleCircleConstraint.h"

#include <iostream>

RigidBody2DSim::RigidBody2DSim(const RigidBody2DState &state,
                               const SpatialGrid &spatial_grid)
    : m_state(state), m_spatial_grid(spatial_grid), m_constraint_cache() {}

RigidBody2DSim::RigidBody2DSim()
    : m_state(), m_spatial_grid(1.0), m_constraint_cache() {}

RigidBody2DState &RigidBody2DSim::getState() { return m_state; }

const RigidBody2DState &RigidBody2DSim::getState() const { return m_state; }

RigidBody2DState &RigidBody2DSim::state() { return m_state; }

const RigidBody2DState &RigidBody2DSim::state() const { return m_state; }

SpatialGrid &RigidBody2DSim::grid() { return m_spatial_grid; }

scalar RigidBody2DSim::computeKineticEnergy() const {
  return 0.5 * m_state.v().dot(m_state.M() * m_state.v());
}

scalar RigidBody2DSim::computePotentialEnergy() const {
  scalar U{0.0};
  const std::vector<std::unique_ptr<RigidBody2DForce>> &forces{
      m_state.forces()};
  for (const std::unique_ptr<RigidBody2DForce> &force : forces) {
    U += force->computePotential(m_state.q(), m_state.M());
  }
  return U;
}

scalar RigidBody2DSim::computeTotalEnergy() const {
  return computeKineticEnergy() + computePotentialEnergy();
}

Vector2s RigidBody2DSim::computeTotalMomentum() const {
  VectorXs p;
  computeMomentum(m_state.v(), p);
  assert(p.size() == 2);
  return p;
}

scalar RigidBody2DSim::computeTotalAngularMomentum() const {
  VectorXs L;
  computeAngularMomentum(m_state.v(), L);
  assert(L.size() == 1);
  return L(0);
}

int RigidBody2DSim::nqdofs() const { return int(m_state.q().size()); }

int RigidBody2DSim::nvdofs() const { return int(m_state.v().size()); }

unsigned RigidBody2DSim::numVelDoFsPerBody() const { return 3; }

unsigned RigidBody2DSim::ambientSpaceDimensions() const { return 2; }

bool RigidBody2DSim::isKinematicallyScripted(const int i) const {
  return m_state.fixed(i);
}

void RigidBody2DSim::computeForce(const VectorXs &q, const VectorXs &v,
                                  const scalar &t, VectorXs &F) {
  assert(q.size() % 3 == 0);
  assert(v.size() == q.size());
  assert(v.size() == F.size());
  F.setZero();
  const std::vector<std::unique_ptr<RigidBody2DForce>> &forces{
      m_state.forces()};
  for (const std::unique_ptr<RigidBody2DForce> &force : forces) {
    force->computeForce(q, v, m_state.M(), F);
  }
}

void RigidBody2DSim::linearInertialConfigurationUpdate(const VectorXs &q0,
                                                       const VectorXs &v0,
                                                       const scalar &dt,
                                                       VectorXs &q1) const {
  assert(q0.size() == v0.size());
  assert(q0.size() == q1.size());
  assert(dt > 0.0);
  q1 = q0 + dt * v0;
}

const SparseMatrixsc &RigidBody2DSim::M() const { return m_state.M(); }

const SparseMatrixsc &RigidBody2DSim::Minv() const { return m_state.Minv(); }

const SparseMatrixsc &RigidBody2DSim::M0() const {
  // Mass matrix is invaraint to configuration for this system
  return m_state.M();
}

const SparseMatrixsc &RigidBody2DSim::Minv0() const {
  // Mass matrix is invaraint to configuration for this system
  return m_state.Minv();
}

void RigidBody2DSim::computeMomentum(const VectorXs &v, VectorXs &p) const {
  p = Vector2s::Zero();
  assert(m_state.q().size() % 3 == 0);
  const unsigned nbodies{static_cast<unsigned>(m_state.q().size() / 3)};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    p += m_state.m(bdy_idx) * v.segment<2>(3 * bdy_idx);
  }
}

void RigidBody2DSim::computeAngularMomentum(const VectorXs &v,
                                            VectorXs &L) const {
  L = VectorXs::Zero(1);
  assert(m_state.q().size() % 3 == 0);
  const unsigned nbodies{static_cast<unsigned>(m_state.q().size() / 3)};
  // Contribution from center of mass
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    L(0) += m_state.m(bdy_idx) *
            MathUtilities::cross(m_state.q().segment<2>(3 * bdy_idx),
                                 v.segment<2>(3 * bdy_idx));
  }
  // Contribution from rotation about center of mass
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    L(0) += m_state.I(bdy_idx) * v(3 * bdy_idx + 2);
  }
}

void RigidBody2DSim::boxBoxNarrowPhaseCollision(
    const unsigned idx0, const unsigned idx1, const BoxGeometry &box0,
    const BoxGeometry &box1, const VectorXs &q0, const VectorXs &q1,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  if (isKinematicallyScripted(idx0) || isKinematicallyScripted(idx1)) {
    std::cerr << "Box-Box kinematic collisions not yet supported." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Note: Detection is at q1...
  const Vector2s x0_t1{q1.segment<2>(3 * idx0)};
  const scalar theta0_t1{q1(3 * idx0 + 2)};
  const Vector2s x1_t1{q1.segment<2>(3 * idx1)};
  const scalar theta1_t1{q1(3 * idx1 + 2)};
  Vector2s n;
  std::vector<Vector2s> points;
  BoxBoxTools::isActive(x0_t1, theta0_t1, box0.r(), x1_t1, theta1_t1, box1.r(),
                        n, points);

  // ... but constraint construction is at q0 to conserve angular momentum
  for (const Vector2s &point : points) {
    active_set.emplace_back(new BodyBodyConstraint{idx0, idx1, point, n, q0});
  }
}

void RigidBody2DSim::boxCircleNarrowPhaseCollision(
    const unsigned idx_crcl, const unsigned idx_box,
    const CircleGeometry &circle, const BoxGeometry &box, const VectorXs &q0,
    const VectorXs &q1, const VectorXs &v,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  if (isKinematicallyScripted(idx_crcl)) {
    std::cerr << "Kinematic circle vs. box collisions not yet supported."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Note: Detection is at q1...
  const Vector2s x0_t1{q1.segment<2>(3 * idx_crcl)};
  const Vector2s x1_t1{q1.segment<2>(3 * idx_box)};
  const scalar theta1_t1{q1(3 * idx_box + 2)};
  Vector2s n;
  Vector2s p;
  const bool is_active{CircleBoxTools::isActive(x0_t1, circle.r(), x1_t1,
                                                theta1_t1, box.r(), n, p)};

  // ... but constraint construction is at q0 to conserve angular momentum
  if (is_active) {
    if (!isKinematicallyScripted(idx_box)) {
      if (idx_crcl < idx_box) {
        active_set.emplace_back(
            new BodyBodyConstraint{idx_crcl, idx_box, p, n, q0});
      } else {
        active_set.emplace_back(
            new BodyBodyConstraint{idx_box, idx_crcl, p, -n, q0});
      }
    } else {
      const Vector2s x{q0.segment<2>(3 * idx_box)};
      const Vector2s vel{v.segment<2>(3 * idx_box)};
      const scalar omega{v(3 * idx_box + 2)};
      active_set.emplace_back(new KinematicBoxCircleConstraint{
          idx_crcl, circle.r(), n, idx_box, x, vel, omega, box.r()});
    }
  }
}

void RigidBody2DSim::annulusCircleNarrowPhaseCollision(
    const unsigned idx0, const unsigned idx1, const AnnulusGeometry &annulus,
    const CircleGeometry &circle, const VectorXs &q0, const VectorXs &q1,
    const VectorXs &v,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  assert(!isKinematicallyScripted(idx0));
  assert(!isKinematicallyScripted(idx1));

  const Vector2s q1a{q1.segment<2>(3 * idx0)};
  const Vector2s q1b{q1.segment<2>(3 * idx1)};

  // Distance between the annulus and circle centers
  const scalar d = (q1a - q1b).norm();

  // Radius of the center of the annulus
  const scalar rc = (annulus.r0() + annulus.r1()) / 2.0;

  // A possible collision against the outside of the annulus
  if (d >= rc) {
    if (d <= annulus.r1() + circle.r()) {
      const Vector2s q0a{q0.segment<2>(3 * idx0)};
      const Vector2s q0b{q0.segment<2>(3 * idx1)};
      const Vector2s n{(q0a - q0b).normalized()};
      const scalar ra{annulus.r1()};
      const scalar rb{circle.r()};
      const Vector2s p{q0a + (ra / (ra + rb)) * (q0b - q0a)};
      active_set.emplace_back(
          new CircleCircleConstraint{idx0, idx1, n, p, ra, rb});
    }
  }
  // A possible collision against the inside of the annulus
  else {
    if (d >= annulus.r0() - circle.r()) {
      const Vector2s q0a{q0.segment<2>(3 * idx0)};
      const Vector2s q0b{q0.segment<2>(3 * idx1)};
      const Vector2s n{(q0a - q0b).normalized()};
      std::cout << "n: " << n.transpose() << std::endl;

      std::cout << "Code up the inside collision part!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

void RigidBody2DSim::dispatchNarrowPhaseCollision(
    unsigned idx0, unsigned idx1, const VectorXs &q0, const VectorXs &q1,
    const VectorXs &v,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  assert(q0.size() % 3 == 0);
  assert(q0.size() == q1.size());

  if (isKinematicallyScripted(idx0) && isKinematicallyScripted(idx1)) {
    return;
  }

  // Ensure that the kinematic body is always listed second
  if (isKinematicallyScripted(idx0)) {
    using std::swap;
    swap(idx0, idx1);
  }

  const std::unique_ptr<RigidBody2DGeometry> &geo0{m_state.bodyGeometry(idx0)};
  const std::unique_ptr<RigidBody2DGeometry> &geo1{m_state.bodyGeometry(idx1)};

  switch (geo0->type()) {
  case RigidBody2DGeometryType::CIRCLE: {
    const CircleGeometry &circle_geo0{static_cast<CircleGeometry &>(*geo0)};
    switch (geo1->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      const CircleGeometry &circle_geo1{static_cast<CircleGeometry &>(*geo1)};

      const Vector2s q0a{q0.segment<2>(3 * idx0)};
      const Vector2s q1a{q1.segment<2>(3 * idx0)};
      const scalar ra{circle_geo0.r()};
      const Vector2s q0b{q0.segment<2>(3 * idx1)};
      const Vector2s q1b{q1.segment<2>(3 * idx1)};
      const scalar rb{circle_geo1.r()};

      const std::pair<bool, scalar> ccd_result{
          CollisionDetectionUtilities::ballBallCCDCollisionHappens(
              q0a, q1a, ra, q0b, q1b, rb)};

      if (ccd_result.first) {
#ifndef NDEBUG
        {
          const Vector2s x0{(1.0 - ccd_result.second) * q0a +
                            ccd_result.second * q1a};
          const Vector2s x1{(1.0 - ccd_result.second) * q0b +
                            ccd_result.second * q1b};
          assert(((x0 - x1).squaredNorm() - (ra + rb) * (ra + rb)) <= 1.0e-9);
        }
#endif

        // Creation of constraints at q0 to preserve angular momentum
        const Vector2s n{(q0a - q0b).normalized()};
        assert(!isKinematicallyScripted(idx0));
        if (!isKinematicallyScripted(idx1)) {
          const Vector2s p{q0a + (ra / (ra + rb)) * (q0b - q0a)};
          active_set.emplace_back(
              new CircleCircleConstraint{idx0, idx1, n, p, ra, rb});
        } else {
          const Vector2s vel{v.segment<2>(3 * idx1)};
          const scalar omega{v(3 * idx1 + 2)};
          // active_set.emplace_back( new KinematicObjectCircleConstraint{ idx0,
          // ra, n, idx1, q0b, vel, omega } );
          active_set.emplace_back(new KinematicCircleCircleConstraint(
              idx0, ra, n, idx1, q0b, vel, omega, rb));
        }
      }
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      assert(!isKinematicallyScripted(idx0));
      const BoxGeometry &box_geo1{static_cast<BoxGeometry &>(*geo1)};
      boxCircleNarrowPhaseCollision(idx0, idx1, circle_geo0, box_geo1, q0, q1,
                                    v, active_set);
      break;
    }
    case RigidBody2DGeometryType::ANNULUS: {
      std::cerr << "Code up circle vs. annulus collision detection!"
                << std::endl;
      std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
      break;
#endif
    }
    }
    break;
  }
  case RigidBody2DGeometryType::BOX: {
    const BoxGeometry &box_geo0{static_cast<BoxGeometry &>(*geo0)};
    switch (geo1->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      assert(!isKinematicallyScripted(idx0));
      assert(!isKinematicallyScripted(idx1));
      const CircleGeometry &circle_geo1{static_cast<CircleGeometry &>(*geo1)};
      boxCircleNarrowPhaseCollision(idx1, idx0, circle_geo1, box_geo0, q0, q1,
                                    v, active_set);
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      assert(!isKinematicallyScripted(idx0));
      assert(!isKinematicallyScripted(idx1));
      const BoxGeometry &box_geo1{static_cast<BoxGeometry &>(*geo1)};
      boxBoxNarrowPhaseCollision(idx0, idx1, box_geo0, box_geo1, q0, q1,
                                 active_set);
      break;
    }
    case RigidBody2DGeometryType::ANNULUS: {
      std::cerr << "Code up box vs. annulus collision detection!" << std::endl;
      std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
      break;
#endif
    }
    }
    break;
  }
  case RigidBody2DGeometryType::ANNULUS: {
    const AnnulusGeometry &annulus_geo0{static_cast<AnnulusGeometry &>(*geo0)};
    switch (geo1->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      assert(!isKinematicallyScripted(idx0));
      assert(!isKinematicallyScripted(idx1));
      const CircleGeometry &circle_geo1{static_cast<CircleGeometry &>(*geo1)};
      annulusCircleNarrowPhaseCollision(idx0, idx1, annulus_geo0, circle_geo1,
                                        q0, q1, v, active_set);
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      std::cerr << "Code up annulus vs. box collision detection!" << std::endl;
      std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
      break;
#endif
    }
    case RigidBody2DGeometryType::ANNULUS: {
      std::cerr << "Code up annulus vs. annulus collision detection!"
                << std::endl;
      std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
      break;
#endif
    }
    }
  }
  }
}

static bool
collisionIsActive(const Vector2s &x0, const scalar &theta0,
                  const std::unique_ptr<RigidBody2DGeometry> &geo0,
                  const Vector2s &x1, const scalar &theta1,
                  const std::unique_ptr<RigidBody2DGeometry> &geo1) {
  switch (geo0->type()) {
  case RigidBody2DGeometryType::CIRCLE: {
    const CircleGeometry &circle_geo0{static_cast<CircleGeometry &>(*geo0)};
    switch (geo1->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      const CircleGeometry &circle_geo1{static_cast<CircleGeometry &>(*geo1)};
      if (CircleCircleConstraint::isActive(x0, x1, circle_geo0.r(),
                                           circle_geo1.r())) {
        return true;
      } else {
        return false;
      }
      // break;
    }
    case RigidBody2DGeometryType::BOX: {
      std::cerr << "CIRCLE-BOX case not handled in collisionIsActive"
                << std::endl;
      std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
      break;
#endif
    }
    case RigidBody2DGeometryType::ANNULUS: {
      std::cerr << "CIRCLE-ANNULUS case not handled in collisionIsActive"
                << std::endl;
      std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
      break;
#endif
    }
    }
    break;
  }
  case RigidBody2DGeometryType::BOX: {
    std::cerr << "BOX case not handled in collisionIsActive" << std::endl;
    std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
    break;
#endif
  }
  case RigidBody2DGeometryType::ANNULUS: {
    std::cerr << "ANNULUS case not handled in collisionIsActive" << std::endl;
    std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
    break;
#endif
  }
// GCC and Intel don't realize that we've exhausted all cases and complain about
// no return here.
#ifndef CMAKE_DETECTED_CLANG_COMPILER
  default: {
    std::cerr << "Invalid geometry type in RigidBody2DSim::collisionIsActive, "
                 "this is a bug."
              << std::endl;
    std::exit(EXIT_FAILURE);
#ifndef CMAKE_DETECTED_CLANG_COMPILER
    break;
#endif
  }
#endif
  }
  return false;
}

static bool collisionIsActive(const unsigned idx0, const unsigned idx1,
                              const std::unique_ptr<RigidBody2DGeometry> &geo0,
                              const std::unique_ptr<RigidBody2DGeometry> &geo1,
                              const VectorXs &q) {
  assert(q.size() % 3 == 0);

  const Vector2s x0{q.segment<2>(3 * idx0)};
  const scalar theta0{q(3 * idx0 + 2)};
  const Vector2s x1{q.segment<2>(3 * idx1)};
  const scalar theta1{q(3 * idx1 + 2)};

  return collisionIsActive(x0, theta0, geo0, x1, theta1, geo1);
}

bool RigidBody2DSim::teleportedCollisionIsActive(
    const TeleportedCollision &teleported_collision,
    const std::unique_ptr<RigidBody2DGeometry> &geo0,
    const std::unique_ptr<RigidBody2DGeometry> &geo1, const VectorXs &q) const {
  assert(q.size() % 3 == 0);

  // Get the center of mass of each body after the teleportation for the end of
  // the step
  Vector2s x0;
  Vector2s x1;
  getTeleportedCollisionCenters(q, teleported_collision, x0, x1);

  const scalar theta0{q(3 * teleported_collision.bodyIndex0() + 2)};
  const scalar theta1{q(3 * teleported_collision.bodyIndex1() + 2)};

  return collisionIsActive(x0, theta0, geo0, x1, theta1, geo1);
}

void RigidBody2DSim::getTeleportedCollisionCenter(const unsigned portal_index,
                                                  const bool portal_plane,
                                                  Vector2s &x) const {
  // TODO: Move this if statement into the portal class
  if (portal_plane == 0) {
    m_state.planarPortals()[portal_index].teleportPointThroughPlaneA(x, x);
  } else {
    m_state.planarPortals()[portal_index].teleportPointThroughPlaneB(x, x);
  }
}

void RigidBody2DSim::getTeleportedCollisionCenters(
    const VectorXs &q, const TeleportedCollision &teleported_collision,
    Vector2s &x0, Vector2s &x1) const {
  assert(q.size() % 3 == 0);

#ifndef NDEBUG
  const unsigned nbodies{static_cast<unsigned>(q.size() / 3)};
#endif

  // Indices of the colliding bodies
  const unsigned idx0{teleported_collision.bodyIndex0()};
  assert(idx0 < nbodies);
  const unsigned idx1{teleported_collision.bodyIndex1()};
  assert(idx1 < nbodies);

  // Indices of the portals
  const unsigned prtl_idx0{teleported_collision.portalIndex0()};
  assert(prtl_idx0 < m_state.planarPortals().size() ||
         prtl_idx0 == std::numeric_limits<unsigned>::max());
  const unsigned prtl_idx1{teleported_collision.portalIndex1()};
  assert(prtl_idx1 < m_state.planarPortals().size() ||
         prtl_idx1 == std::numeric_limits<unsigned>::max());

  // If the first object was teleported
  assert(3 * idx0 + 1 < q.size());
  x0 = q.segment<2>(3 * idx0);
  if (prtl_idx0 != std::numeric_limits<unsigned>::max()) {
    getTeleportedCollisionCenter(prtl_idx0, teleported_collision.plane0(), x0);
  }

  // If the second object was teleported
  assert(3 * idx1 + 1 < q.size());
  x1 = q.segment<2>(3 * idx1);
  if (prtl_idx1 != std::numeric_limits<unsigned>::max()) {
    getTeleportedCollisionCenter(prtl_idx1, teleported_collision.plane1(), x1);
  }
}

void RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision(
    const TeleportedCollision &teleported_collision,
    const std::unique_ptr<RigidBody2DGeometry> &geo0,
    const std::unique_ptr<RigidBody2DGeometry> &geo1, const VectorXs &q0,
    const VectorXs &q1,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  assert(q0.size() % 3 == 0);
  assert(q0.size() == q1.size());

  if (isKinematicallyScripted(teleported_collision.bodyIndex0()) ||
      isKinematicallyScripted(teleported_collision.bodyIndex1())) {
    std::cerr
        << "!!!!!!!!!!! Kinematic geometry not supported with periodic "
           "boundary conditions for 2D rigid body sim, skipping this collision."
        << std::endl;
    return;
    // std::exit( EXIT_FAILURE );
  }

  // Get the center of mass of each body after the teleportation for the start
  // and end of the step
  Vector2s x0_t0;
  Vector2s x1_t0;
  getTeleportedCollisionCenters(q0, teleported_collision, x0_t0, x1_t0);
  const Vector2s delta0_t0{
      x0_t0 - q0.segment<2>(3 * teleported_collision.bodyIndex0())};
  const Vector2s delta1_t0{
      x1_t0 - q0.segment<2>(3 * teleported_collision.bodyIndex1())};

  Vector2s x0_t1;
  Vector2s x1_t1;
  getTeleportedCollisionCenters(q1, teleported_collision, x0_t1, x1_t1);
// TODO: If Lees-Edwards conditions are updated to have different locations at
// start and end of step
//       (instead of same at each, as now) these tests will no longer hold
#ifndef NDEBUG
  {
    const Vector2s delta0_t1{
        x0_t1 - q1.segment<2>(3 * teleported_collision.bodyIndex0())};
    assert((delta0_t0 - delta0_t1).lpNorm<Eigen::Infinity>() <= 1.0e-6);
    const Vector2s delta1_t1{
        x1_t1 - q1.segment<2>(3 * teleported_collision.bodyIndex1())};
    assert((delta1_t0 - delta1_t1).lpNorm<Eigen::Infinity>() <= 1.0e-6);
  }
#endif

  // At least one of the bodies must have been teleported
  assert(teleported_collision.portalIndex0() !=
             std::numeric_limits<unsigned>::max() ||
         teleported_collision.portalIndex1() !=
             std::numeric_limits<unsigned>::max());

  // Determine whether the first body was teleported
  bool portal0_is_lees_edwards{false};
#ifndef NDEBUG
  bool first_was_teleported{false};
#endif
  if (teleported_collision.portalIndex0() !=
      std::numeric_limits<unsigned>::max()) {
#ifndef NDEBUG
    first_was_teleported = true;
#endif
    assert(teleported_collision.portalIndex0() <
           m_state.planarPortals().size());
    if (m_state.planarPortals()[teleported_collision.portalIndex0()]
            .isLeesEdwards()) {
      portal0_is_lees_edwards = true;
    }
  }

  // Determine whether the second body was teleported
  bool portal1_is_lees_edwards{false};
#ifndef NDEBUG
  bool second_was_teleported{false};
#endif
  if (teleported_collision.portalIndex1() !=
      std::numeric_limits<unsigned>::max()) {
#ifndef NDEBUG
    second_was_teleported = true;
#endif
    assert(teleported_collision.portalIndex1() <
           m_state.planarPortals().size());
    if (m_state.planarPortals()[teleported_collision.portalIndex1()]
            .isLeesEdwards()) {
      portal1_is_lees_edwards = true;
    }
  }

  // If neither portal is Lees-Edwards
  if (!portal0_is_lees_edwards && !portal1_is_lees_edwards) {
    switch (geo0->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      const CircleGeometry &circle_geo0{static_cast<CircleGeometry &>(*geo0)};
      switch (geo1->type()) {
      case RigidBody2DGeometryType::CIRCLE: {
        const CircleGeometry &circle_geo1{static_cast<CircleGeometry &>(*geo1)};
        if (CircleCircleConstraint::isActive(x0_t1, x1_t1, circle_geo0.r(),
                                             circle_geo1.r())) {
          // Creation of constraints at q0 to preserve angular momentum
          active_set.emplace_back(new TeleportedCircleCircleConstraint{
              teleported_collision.bodyIndex0(),
              teleported_collision.bodyIndex1(), x0_t0, x1_t0, circle_geo0.r(),
              circle_geo1.r(), delta0_t0, delta1_t0, circle_geo0.r(),
              circle_geo1.r()});
        }
        break;
      }
      case RigidBody2DGeometryType::BOX: {
        std::cerr << "CIRCLE-BOX case not handled in "
                     "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      case RigidBody2DGeometryType::ANNULUS: {
        std::cerr << "CIRCLE-ANNULUS case not handled in "
                     "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      }
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      std::cerr << "BOX case not handled in "
                   "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    case RigidBody2DGeometryType::ANNULUS: {
      std::cerr << "ANNULUS case not handled in "
                   "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    }
  }
  // Otherwise, there is a relative velocity contribution from the Lees-Edwards
  // boundary condition
  else {
    Vector2s kinematic_kick;
    assert(portal0_is_lees_edwards != portal1_is_lees_edwards);
    if (portal1_is_lees_edwards) {
      assert(second_was_teleported);
      // N.B. q1 because collision detection was performed with q1
      Array2s min;
      Array2s max;
      geo1->computeAABB(q1.segment<2>(3 * teleported_collision.bodyIndex1()),
                        q1(3 * teleported_collision.bodyIndex1() + 2), min,
                        max);
      kinematic_kick =
          m_state.planarPortals()[teleported_collision.portalIndex1()]
              .getKinematicVelocityOfAABB(min, max);
    } else // portal0_is_lees_edwards
    {
      assert(first_was_teleported);
      // N.B. q1 because collision detection was performed with q1
      Array2s min;
      Array2s max;
      geo0->computeAABB(q1.segment<2>(3 * teleported_collision.bodyIndex0()),
                        q1(3 * teleported_collision.bodyIndex0() + 2), min,
                        max);
      kinematic_kick =
          -m_state.planarPortals()[teleported_collision.portalIndex0()]
               .getKinematicVelocityOfAABB(min, max);
    }
    switch (geo0->type()) {
    case RigidBody2DGeometryType::CIRCLE: {
      const CircleGeometry &circle_geo0{static_cast<CircleGeometry &>(*geo0)};
      switch (geo1->type()) {
      case RigidBody2DGeometryType::CIRCLE: {
        const CircleGeometry &circle_geo1{static_cast<CircleGeometry &>(*geo1)};
        if (CircleCircleConstraint::isActive(x0_t1, x1_t1, circle_geo0.r(),
                                             circle_geo1.r())) {
          // Creation of constraints at q0 to preserve angular momentum
          active_set.emplace_back(new KinematicKickCircleCircleConstraint{
              teleported_collision.bodyIndex0(),
              teleported_collision.bodyIndex1(), x0_t0, x1_t0, circle_geo0.r(),
              circle_geo1.r(), kinematic_kick});
        }
        break;
      }
      case RigidBody2DGeometryType::BOX: {
        std::cerr << "CIRCLE-BOX case not handled in "
                     "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      case RigidBody2DGeometryType::ANNULUS: {
        std::cerr << "CIRCLE-ANNULUS case not handled in "
                     "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      }
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      std::cerr << "BOX case not handled in "
                   "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    case RigidBody2DGeometryType::ANNULUS: {
      std::cerr << "ANNULUS case not handled in "
                   "RigidBody2DSim::dispatchTeleportedNarrowPhaseCollision"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    }
  }
}

void RigidBody2DSim::computeBodyPlaneActiveSetAllPairs(
    const VectorXs &q0, const VectorXs &q1,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  assert(q0.size() % 3 == 0);
  assert(q0.size() == q1.size());

  const unsigned nbodies{static_cast<unsigned>(q0.size() / 3)};

  // Check all body-plane pairs
  for (unsigned plane_idx = 0; plane_idx < m_state.planes().size();
       ++plane_idx) {
    const RigidBody2DStaticPlane &plane{m_state.planes()[plane_idx]};
    for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
      // Skip kinematically scripted bodies
      if (isKinematicallyScripted(bdy_idx)) {
        continue;
      }

      switch (m_state.geometry()[m_state.geometryIndices()(bdy_idx)]->type()) {
      case RigidBody2DGeometryType::CIRCLE: {
        const CircleGeometry &circle_geo{static_cast<CircleGeometry &>(
            *m_state.geometry()[m_state.geometryIndices()(bdy_idx)])};
        if (StaticPlaneCircleConstraint::isActive(q1.segment<2>(3 * bdy_idx),
                                                  circle_geo.r(), plane)) {
          active_set.emplace_back(new StaticPlaneCircleConstraint{
              bdy_idx, plane_idx, circle_geo.r(), plane});
        }
        break;
      }
      case RigidBody2DGeometryType::BOX: {
        // TODO: Make this faster, if needed
        const BoxGeometry &box_geo{static_cast<BoxGeometry &>(
            *m_state.geometry()[m_state.geometryIndices()(bdy_idx)])};

        const Vector2s x{q1.segment<2>(3 * bdy_idx)};
        const Eigen::Rotation2D<scalar> R{q1(3 * bdy_idx + 2)};
        const Array2s r{box_geo.r()};

        // Check each vertex of the box
        for (int i = -1; i < 2; i += 2) {
          for (int j = -1; j < 2; j += 2) {
            const Vector2s body_space_arm{(Array2s{i, j} * r).matrix()};
            const Vector2s transformed_vertex{x + R * body_space_arm};
            const scalar dist{plane.n().dot(transformed_vertex - plane.x())};
            if (dist <= 0.0) {
              active_set.emplace_back(new StaticPlaneBodyConstraint{
                  bdy_idx, body_space_arm, plane, plane_idx});
            }
          }
        }
        break;
      }
      case RigidBody2DGeometryType::ANNULUS: {
        const AnnulusGeometry &annulus_geo{static_cast<AnnulusGeometry &>(
            *m_state.geometry()[m_state.geometryIndices()(bdy_idx)])};
        if (StaticPlaneCircleConstraint::isActive(q1.segment<2>(3 * bdy_idx),
                                                  annulus_geo.r1(), plane)) {
          active_set.emplace_back(new StaticPlaneCircleConstraint{
              bdy_idx, plane_idx, annulus_geo.r1(), plane});
        }
        break;
      }
      }
    }
  }
}

void RigidBody2DSim::computeBodyDrumActiveSetAllPairs(
    const VectorXs &q0, const VectorXs &q1,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  assert(q0.size() % 3 == 0);
  assert(q0.size() == q1.size());

  const unsigned nbodies{static_cast<unsigned>(q0.size() / 3)};

  // Check all body-drum pairs
  for (unsigned drum_idx = 0; drum_idx < m_state.drums().size(); drum_idx++) {
    const RigidBody2DStaticDrum &drum{m_state.drums()[drum_idx]};
    for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
      // Skip kinematically scripted bodies
      if (isKinematicallyScripted(bdy_idx)) {
        continue;
      }

      switch (m_state.geometry()[m_state.geometryIndices()(bdy_idx)]->type()) {
      case RigidBody2DGeometryType::CIRCLE: {
        const CircleGeometry &circle_geo{static_cast<CircleGeometry &>(
            *m_state.geometry()[m_state.geometryIndices()(bdy_idx)])};
        if (StaticDrumCircleConstraint::isActive(q1.segment<2>(3 * bdy_idx),
                                                 circle_geo.r(), drum.x(),
                                                 drum.r())) {
          active_set.emplace_back(new StaticDrumCircleConstraint{
              bdy_idx, drum_idx, circle_geo.r(), drum});
        }
        break;
      }
      case RigidBody2DGeometryType::BOX: {
        std::cerr << "Box vs. static drum not yet coded up!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      case RigidBody2DGeometryType::ANNULUS: {
        std::cerr << "Annulus vs. static drum not yet coded up!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      }
    }
  }
}

void RigidBody2DSim::computeActiveSet(
    const VectorXs &q0, const VectorXs &qp, const VectorXs &v,
    const bool reduce_bandwidth,
    std::vector<std::unique_ptr<Constraint>> &active_set) {
  assert(q0.size() % 3 == 0);
  assert(q0.size() == qp.size());

  active_set.clear();

  // Detect body-body collisions
  if (m_state.planarPortals().empty()) {
    computeBodyBodyActiveSetSpatialGrid(q0, qp, v, active_set);
  } else {
    computeBodyBodyActiveSetSpatialGridWithPortals(q0, qp, v, active_set);
  }

  // Check all body-plane pairs
  computeBodyPlaneActiveSetAllPairs(q0, qp, active_set);

  // Check all body-drum pairs
  computeBodyDrumActiveSetAllPairs(q0, qp, active_set);

  if (reduce_bandwidth && !active_set.empty()) {
    // Reorder the constraints to reduce the bandwidth
    ContactGraphTools::reduceBandwidth(numBodies(), active_set);
  }
}

void RigidBody2DSim::computeImpactBases(
    const VectorXs &q,
    const std::vector<std::unique_ptr<Constraint>> &active_set,
    MatrixXXsc &impact_bases) const {
  const unsigned ncols{static_cast<unsigned>(active_set.size())};
  impact_bases.resize(2, ncols);
  for (unsigned col_num = 0; col_num < ncols; ++col_num) {
    VectorXs current_normal;
    active_set[col_num]->getWorldSpaceContactNormal(q, current_normal);
    assert(fabs(current_normal.norm() - 1.0) <= 1.0e-6);
    impact_bases.col(col_num) = current_normal;
  }
}

void RigidBody2DSim::computeContactBases(
    const VectorXs &q, const VectorXs &v,
    const std::vector<std::unique_ptr<Constraint>> &active_set,
    MatrixXXsc &contact_bases) const {
  const unsigned ncols{static_cast<unsigned>(active_set.size())};
  contact_bases.resize(2, 2 * ncols);
  for (unsigned col_num = 0; col_num < ncols; ++col_num) {
    MatrixXXsc basis;
    active_set[col_num]->computeBasis(q, v, basis);
    assert(basis.rows() == basis.cols());
    assert(basis.rows() == 2);
    assert((basis * basis.transpose() - MatrixXXsc::Identity(2, 2))
               .lpNorm<Eigen::Infinity>() <= 1.0e-6);
    assert(fabs(basis.determinant() - 1.0) <= 1.0e-6);
    contact_bases.block<2, 2>(0, 2 * col_num) = basis;
  }
}

void RigidBody2DSim::clearConstraintCache() { m_constraint_cache.clear(); }

void RigidBody2DSim::cacheConstraint(const Constraint &constraint,
                                     const VectorXs &r) {
  m_constraint_cache.cacheConstraint(constraint, r);
}

void RigidBody2DSim::getCachedConstraintImpulse(const Constraint &constraint,
                                                VectorXs &r) const {
  m_constraint_cache.getCachedConstraint(constraint, r);
}

bool RigidBody2DSim::constraintCacheEmpty() const {
  return m_constraint_cache.empty();
}

void RigidBody2DSim::flow(PythonScripting &call_back, const unsigned iteration,
                          const Rational<std::intmax_t> &dt,
                          UnconstrainedMap &umap) {
  call_back.setState(m_state);
  call_back.startOfStepCallback(iteration, dt);
  call_back.forgetState();

  VectorXs q1{m_state.q().size()};
  VectorXs v1{m_state.v().size()};

  updatePeriodicBoundaryConditionsStartOfStep(iteration, scalar(dt));

  umap.flow(m_state.q(), m_state.v(), *this, iteration, scalar(dt), q1, v1);

  q1.swap(m_state.q());
  v1.swap(m_state.v());

  enforcePeriodicBoundaryConditions(m_state.q(), m_state.v());

  call_back.setState(m_state);
  call_back.endOfStepCallback(iteration, dt);
  call_back.forgetState();
}

void RigidBody2DSim::flow(PythonScripting &call_back, const unsigned iteration,
                          const Rational<std::intmax_t> &dt,
                          UnconstrainedMap &umap, ImpactOperator &iop,
                          const scalar &CoR, ImpactMap &imap,
                          const bool reduce_bandwidth) {
  call_back.setState(m_state);
  call_back.startOfStepCallback(iteration, dt);
  call_back.forgetState();

  VectorXs q1{m_state.q().size()};
  VectorXs v1{m_state.v().size()};

  updatePeriodicBoundaryConditionsStartOfStep(iteration, scalar(dt));

  imap.flow(call_back, *this, *this, umap, iop, iteration, scalar(dt), CoR,
            reduce_bandwidth, m_state.q(), m_state.v(), q1, v1);

  q1.swap(m_state.q());
  v1.swap(m_state.v());

  enforcePeriodicBoundaryConditions(m_state.q(), m_state.v());

  call_back.setState(m_state);
  call_back.endOfStepCallback(iteration, dt);
  call_back.forgetState();
}

void RigidBody2DSim::flow(PythonScripting &call_back, const unsigned iteration,
                          const Rational<std::intmax_t> &dt,
                          UnconstrainedMap &umap, const scalar &CoR,
                          const scalar &mu, FrictionSolver &solver,
                          ImpactFrictionMap &ifmap,
                          const bool reduce_bandwidth) {
  call_back.setState(m_state);
  call_back.startOfStepCallback(iteration, dt);
  call_back.forgetState();

  VectorXs q1{m_state.q().size()};
  VectorXs v1{m_state.v().size()};

  updatePeriodicBoundaryConditionsStartOfStep(iteration, scalar(dt));

  ifmap.flow(call_back, *this, *this, umap, solver, iteration, scalar(dt), CoR,
             mu, reduce_bandwidth, m_state.q(), m_state.v(), q1, v1);

  q1.swap(m_state.q());
  v1.swap(m_state.v());

  enforcePeriodicBoundaryConditions(m_state.q(), m_state.v());

  call_back.setState(m_state);
  call_back.endOfStepCallback(iteration, dt);
  call_back.forgetState();
}

void RigidBody2DSim::flowWithWeight(PythonScripting &call_back,
                                    const unsigned iteration,
                                    const Rational<std::intmax_t> &dt,
                                    UnconstrainedMap &umap, const scalar &CoR,
                                    const scalar &mu, FrictionSolver &solver,
                                    ImpactFrictionMap &ifmap,
                                    const bool reduce_bandwidth) {
  call_back.setState(m_state);
  call_back.startOfStepCallback(iteration, dt);
  call_back.forgetState();

  VectorXs q1{m_state.q().size()};
  VectorXs v1{m_state.v().size()};

  updatePeriodicBoundaryConditionsStartOfStep(iteration, scalar(dt));

  ifmap.flowWithWeights(call_back, *this, *this, umap, solver, iteration,
                        scalar(dt), CoR, mu, reduce_bandwidth, m_state.q(),
                        m_state.v(), m_state.massWeights(), q1, v1);

  q1.swap(m_state.q());
  v1.swap(m_state.v());

  enforcePeriodicBoundaryConditions(m_state.q(), m_state.v());

  call_back.setState(m_state);
  call_back.endOfStepCallback(iteration, dt);
  call_back.forgetState();
}

void RigidBody2DSim::updatePeriodicBoundaryConditionsStartOfStep(
    const unsigned next_iteration, const scalar &dt) {
  for (PlanarPortal &planar_portal : m_state.planarPortals()) {
    planar_portal.updateMovingPortals(dt);
  }
}

void RigidBody2DSim::enforcePeriodicBoundaryConditions(VectorXs &q,
                                                       VectorXs &v) {
  assert(q.size() % 3 == 0);
  assert(q.size() == v.size());

  const unsigned nbodies{static_cast<unsigned>(q.size() / 3)};

  // TODO: Probably faster to invert the loop here, only cache xin once per body
  // For each portal
  for (const PlanarPortal &planar_portal : m_state.planarPortals()) {
    // For each body
    for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
      const Vector2s xin{q.segment<2>(3 * bdy_idx)};
      // TODO: Calling pointInsidePortal and teleportPointInsidePortal is a bit
      // redundant, clean this up! If the body is inside a portal
      if (planar_portal.pointInsidePortal(xin)) {
        // Teleport to the other side of the portal
        Vector2s x_out;
        planar_portal.teleportPointInsidePortal(xin, x_out);
        q.segment<2>(3 * bdy_idx) = x_out;
        // TODO: This check probably isn't needed, additional_vel should be 0
        // for non-LE portals Lees-Edwards Boundary conditions also update the
        // velocity
        if (planar_portal.isLeesEdwards()) {
          const Vector2s additional_vel{
              planar_portal.getKinematicVelocityOfPoint(xin)};
          v.segment<2>(3 * bdy_idx) += additional_vel;
        }
      }
    }
  }
}

void RigidBody2DSim::getAllCircleBodies(Matrix2Xsc &pos, VectorXs &radii,
                                        VectorXu &bdy_idx) const {
  state().getAllCircleBodies(pos, radii, bdy_idx);
}

void RigidBody2DSim::getAllNonFixedCircleBodies(Matrix2Xsc &pos,
                                                VectorXs &radii,
                                                VectorXu &bdy_idx) const {
  state().getAllNonFixedCircleBodies(pos, radii, bdy_idx);
}

VectorXs RigidBody2DSim::getAllRadii() const {
  VectorXs radii{numBodies()};
  for (int bdy_idx = 0; bdy_idx < radii.size(); ++bdy_idx) {
    const std::unique_ptr<RigidBody2DGeometry> &current_geo{
        state().bodyGeometry(bdy_idx)};
    if (current_geo->type() == RigidBody2DGeometryType::CIRCLE) {
      const CircleGeometry &circle_geo{
          static_cast<const CircleGeometry &>(*current_geo)};
      radii(bdy_idx) = circle_geo.r();
    } else {
      radii(bdy_idx) = SCALAR_NAN;
    }
  }
  assert((radii.array() > 0.0).all());
  return radii;
}
//*/

void RigidBody2DSim::computeBodyBodyActiveSetSpatialGrid(
    const VectorXs &q0, const VectorXs &q1, const VectorXs &v,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  assert(q0.size() % 3 == 0);
  assert(q0.size() == q1.size());

  const unsigned nbodies{static_cast<unsigned>(q0.size() / 3)};

  // Candidate bodies that might overlap
  std::set<std::pair<unsigned, unsigned>> possible_overlaps;
  {
    // Compute an AABB for each body
    std::vector<AABB> aabbs;
    aabbs.reserve(nbodies);
    for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
      Array2s min;
      Array2s max;
      m_state.bodyGeometry(bdy_idx)->computeCollisionAABB(
          q0.segment<2>(3 * bdy_idx), q0(3 * bdy_idx + 2),
          q1.segment<2>(3 * bdy_idx), q1(3 * bdy_idx + 2), min, max);
      aabbs.emplace_back(min, max);
    }
    assert(aabbs.size() == nbodies);

    // Determine which bodies possibly overlap
    m_spatial_grid.getPotentialOverlaps(aabbs, possible_overlaps);
  }

  // Create constraints for bodies that actually overlap
  for (const auto &possible_overlap_pair : possible_overlaps) {
    assert(possible_overlap_pair.first < nbodies);
    assert(possible_overlap_pair.second < nbodies);

    // We can run standard narrow phase
    dispatchNarrowPhaseCollision(possible_overlap_pair.first,
                                 possible_overlap_pair.second, q0, q1, v,
                                 active_set);
  }
}

void RigidBody2DSim::computeBodyBodyActiveSetSpatialGridWithPortals(
    const VectorXs &q0, const VectorXs &q1, const VectorXs &v,
    std::vector<std::unique_ptr<Constraint>> &active_set) const {
  assert(q0.size() % 3 == 0);
  assert(q0.size() == q1.size());

  const unsigned nbodies{static_cast<unsigned>(q0.size() / 3)};

  // Candidate bodies that might overlap
  std::set<std::pair<unsigned, unsigned>> possible_overlaps;
  // Map from teleported AABB indices and body and portal indices
  std::map<unsigned, TeleportedBody> teleported_aabb_body_indices;
  {
    // Compute an AABB for each body
    std::vector<AABB> aabbs;
    aabbs.reserve(nbodies);
    for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
      Array2s min;
      Array2s max;
      m_state.bodyGeometry(bdy_idx)->computeAABB(q1.segment<2>(3 * bdy_idx),
                                                 q1(3 * bdy_idx + 2), min, max);
      aabbs.emplace_back(min, max);
    }
    assert(aabbs.size() == nbodies);

    // Compute an AABB for each teleported body
    auto aabb_bdy_map_itr = teleported_aabb_body_indices.cbegin();
    // For each portal
    for (const PlanarPortal &planar_portal : m_state.planarPortals()) {
      // For each body
      for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
        // If the body is inside a portal
        bool intersecting_plane_index;
        if (planar_portal.aabbTouchesPortal(aabbs[bdy_idx].min(),
                                            aabbs[bdy_idx].max(),
                                            intersecting_plane_index)) {
          // Teleport to the other side of the portal
          Vector2s x_out;
          planar_portal.teleportPoint(q1.segment<2>(3 * bdy_idx),
                                      intersecting_plane_index, x_out);
          // Compute an AABB for the teleported particle
          Array2s min;
          Array2s max;
          m_state.bodyGeometry(bdy_idx)->computeAABB(x_out, q1(3 * bdy_idx + 2),
                                                     min, max);
          aabbs.emplace_back(min, max);

          const unsigned prtl_idx{
              Utilities::index(m_state.planarPortals(), planar_portal)};
          aabb_bdy_map_itr = teleported_aabb_body_indices.insert(
              aabb_bdy_map_itr,
              std::make_pair(
                  aabbs.size() - 1,
                  TeleportedBody{bdy_idx, prtl_idx, intersecting_plane_index}));
        }
      }
    }

    // Determine which bodies possibly overlap
    m_spatial_grid.getPotentialOverlaps(aabbs, possible_overlaps);
  }

  std::set<TeleportedCollision> teleported_collisions;

#ifndef NDEBUG
  std::vector<std::pair<unsigned, unsigned>> duplicate_indices;
#endif

  // Create constraints for bodies that actually overlap
  for (const auto &possible_overlap_pair : possible_overlaps) {
    const bool first_teleported{possible_overlap_pair.first >= nbodies};
    const bool second_teleported{possible_overlap_pair.second >= nbodies};

    // If neither body in the current collision was teleported
    if (!first_teleported && !second_teleported) {
      // We can run standard narrow phase
      dispatchNarrowPhaseCollision(possible_overlap_pair.first,
                                   possible_overlap_pair.second, q0, q1, v,
                                   active_set);
    }
    // If at least one of the balls was teleported
    else {
      unsigned bdy_idx_0{possible_overlap_pair.first};
      unsigned bdy_idx_1{possible_overlap_pair.second};
      unsigned prtl_idx_0{std::numeric_limits<unsigned>::max()};
      unsigned prtl_idx_1{std::numeric_limits<unsigned>::max()};
      bool prtl_plane_0{0};
      bool prtl_plane_1{0};

      if (first_teleported) {
        using itr_type = std::map<unsigned, TeleportedBody>::const_iterator;
        const itr_type map_itr{
            teleported_aabb_body_indices.find(possible_overlap_pair.first)};
        assert(map_itr != teleported_aabb_body_indices.cend());
        bdy_idx_0 = map_itr->second.bodyIndex();
        assert(bdy_idx_0 < nbodies);
        prtl_idx_0 = map_itr->second.portalIndex();
        assert(prtl_idx_0 < m_state.planarPortals().size());
        prtl_plane_0 = map_itr->second.planeIndex();
      }
      if (second_teleported) {
        using itr_type = std::map<unsigned, TeleportedBody>::const_iterator;
        const itr_type map_itr{
            teleported_aabb_body_indices.find(possible_overlap_pair.second)};
        assert(map_itr != teleported_aabb_body_indices.cend());
        bdy_idx_1 = map_itr->second.bodyIndex();
        assert(bdy_idx_1 < nbodies);
        prtl_idx_1 = map_itr->second.portalIndex();
        assert(prtl_idx_1 < m_state.planarPortals().size());
        prtl_plane_1 = map_itr->second.planeIndex();
      }

      // Check if the collision will be detected in the unteleported state
      if (first_teleported && second_teleported) {
        if (collisionIsActive(bdy_idx_0, bdy_idx_1,
                              m_state.bodyGeometry(bdy_idx_0),
                              m_state.bodyGeometry(bdy_idx_1), q1)) {
#ifndef NDEBUG
          duplicate_indices.push_back(std::make_pair(bdy_idx_0, bdy_idx_1));
#endif
          continue;
        }
      }

      // Check if the teleported collision happens
      if (isKinematicallyScripted(bdy_idx_0) &&
          isKinematicallyScripted(bdy_idx_1)) {
        continue;
      }
      const TeleportedCollision possible_collision{bdy_idx_0,    bdy_idx_1,
                                                   prtl_idx_0,   prtl_idx_1,
                                                   prtl_plane_0, prtl_plane_1};
      if (teleportedCollisionIsActive(possible_collision,
                                      m_state.bodyGeometry(bdy_idx_0),
                                      m_state.bodyGeometry(bdy_idx_1), q1)) {
        teleported_collisions.insert(possible_collision);
      }
    }
  }
  possible_overlaps.clear();
  teleported_aabb_body_indices.clear();

#ifndef NDEBUG
  // Double check that non-teleport duplicate collisions were actually
  // duplicates
  for (const auto &dup_col : duplicate_indices) {
    bool entry_found{false};
    const unsigned dup_idx_0{std::min(dup_col.first, dup_col.second)};
    const unsigned dup_idx_1{std::max(dup_col.first, dup_col.second)};
    for (const std::unique_ptr<Constraint> &col : active_set) {
      std::pair<int, int> bodies;
      col->getBodyIndices(bodies);
      if (bodies.first < 0 || bodies.second < 0) {
        continue;
      }
      const unsigned idx_0{
          std::min(unsigned(bodies.first), unsigned(bodies.second))};
      const unsigned idx_1{
          std::max(unsigned(bodies.first), unsigned(bodies.second))};
      if (dup_idx_0 == idx_0 && dup_idx_1 == idx_1) {
        entry_found = true;
        break;
      }
    }
    assert(entry_found);
  }
#endif

  // Create constraints for teleported collisions
  for (const TeleportedCollision &teleported_collision :
       teleported_collisions) {
    assert(teleported_collision.bodyIndex0() < nbodies);
    assert(teleported_collision.bodyIndex1() < nbodies);
    assert(teleported_collision.bodyIndex0() !=
           teleported_collision.bodyIndex1());
    dispatchTeleportedNarrowPhaseCollision(
        teleported_collision,
        m_state.bodyGeometry(teleported_collision.bodyIndex0()),
        m_state.bodyGeometry(teleported_collision.bodyIndex1()), q0, q1,
        active_set);
  }
}

// TODO: 0 size plane matrices are not output due to a bug in an older version
// of HDF5
void RigidBody2DSim::writeBinaryState(HDF5File &output_file) const {
#ifdef USE_HDF5
  // Output the configuration
  output_file.writeMatrix("discrete", "q", m_state.q());
  // Output the velocity
  output_file.writeMatrix("discrete", "v", m_state.v());
  // Output the mass
  {
    // Assemble the mass into a single flat vector like q, v, and r
    assert(m_state.M().nonZeros() == m_state.q().size());
    const VectorXs m{Eigen::Map<const VectorXs>{&m_state.M().data().value(0),
                                                m_state.q().size()}};
    output_file.writeMatrix("discrete", "m", m);
  }
  {
    VectorXu fixed{numBodies()};
    for (int body_index = 0; body_index < fixed.size(); ++body_index) {
      fixed(body_index) = isKinematicallyScripted(body_index) ? 1 : 0;
    }
    output_file.writeMatrix("", "discrete/kinematically_scripted", fixed);
  }
  {
    VectorXi unique_ids{numBodies()};
    for (int body_index = 0; body_index < unique_ids.size(); body_index++) {
      unique_ids(body_index) = m_state.getUniqueBodyIndex(body_index);
    }
    output_file.writeMatrix("", "discrete/unique_index", unique_ids);
  }
  // Output the simulated geometry
  RigidBody2DStateOutput::writeGeometryIndices(
      m_state.geometry(), m_state.geometryIndices(), "discrete/geometry",
      output_file);
  RigidBody2DStateOutput::writeGeometry(m_state.geometry(), "discrete/geometry",
                                        output_file);
  // Output the static geometry
  if (!m_state.planes().empty()) {
    RigidBody2DStateOutput::writeStaticPlanes(
        m_state.planes(), "discrete/static_geometry", output_file);
  }
  if (!m_state.drums().empty()) {
    RigidBody2DStateOutput::writeStaticDrums(
        m_state.drums(), "discrete/static_geometry", output_file);
  }
  if (!m_state.planarPortals().empty()) {
    RigidBody2DStateOutput::writePlanarPortals(
        m_state.planarPortals(), "discrete/static_geometry", output_file);
  }
  // For hybridization
  output_file.writeMatrix("discrete", "hybrid_factors",
                          m_state.hybridFactors());
  {
    VectorXu always_fixed{numBodies()};
    for (int body_index = 0; body_index < always_fixed.size(); ++body_index) {
      always_fixed(body_index) = m_state.alwaysFixed(body_index) ? 1 : 0;
    }
    output_file.writeMatrix("", "discrete/always_fixed", always_fixed);
  }
#else
  std::cerr << "Error, RigidBody2DSim::writeBinaryState requires HDF5 support."
            << std::endl;
  std::exit(EXIT_FAILURE);
#endif
}

void RigidBody2DSim::writeBinaryState(const std::string &prefix,
                                      HDF5File &output_file) const {
#ifdef USE_HDF5
  const std::string output_path = prefix.empty() ? "" : prefix + "/";
  // Output the simulated geometry
  RigidBody2DStateOutput::writeGeometryIndices(
      m_state.geometry(), m_state.geometryIndices(), output_path + "geometry",
      output_file);
  RigidBody2DStateOutput::writeGeometry(m_state.geometry(),
                                        output_path + "geometry", output_file);
  // Output the static geometry
  if (!m_state.planes().empty()) {
    RigidBody2DStateOutput::writeStaticPlanes(
        m_state.planes(), output_path + "static_geometry", output_file);
  }
  if (!m_state.drums().empty()) {
    RigidBody2DStateOutput::writeStaticDrums(
        m_state.drums(), output_path + "static_geometry", output_file);
  }
  if (!m_state.planarPortals().empty()) {
    RigidBody2DStateOutput::writePlanarPortals(
        m_state.planarPortals(), output_path + "static_geometry", output_file);
  }
  // Write out the state of each body
  output_file.writeMatrix(output_path + "state", "q", m_state.q());
  output_file.writeMatrix(output_path + "state", "v", m_state.v());
  // Output the mass
  {
    // Assemble the mass into a single flat vector like q, v, and r
    assert(m_state.M().nonZeros() == m_state.q().size());
    const VectorXs m{Eigen::Map<const VectorXs>{&m_state.M().data().value(0),
                                                m_state.q().size()}};
    output_file.writeMatrix(output_path + "state", "m", m);
  }
  {
    VectorXu fixed{m_state.nbodies()};
    for (unsigned body_index = 0; body_index < m_state.nbodies();
         ++body_index) {
      fixed(body_index) = m_state.alwaysFixed(body_index) ? 1 : 0;
    }
    output_file.writeMatrix(output_path + "state", "kinematically_scripted",
                            fixed);
  }
#else
  std::cerr << "Error, RigidBody3DSim::writeBinaryState requires HDF5 support."
            << std::endl;
  std::exit(EXIT_FAILURE);
#endif
}

void RigidBody2DSim::serialize(std::ostream &output_stream) const {
  assert(output_stream.good());
  m_state.serialize(output_stream);
  m_spatial_grid.serialize(output_stream);
  m_constraint_cache.serialize(output_stream);
}

void RigidBody2DSim::deserialize(std::istream &input_stream) {
  assert(input_stream.good());
  m_state.deserialize(input_stream);
  m_spatial_grid.deserialize(input_stream);
  m_constraint_cache.deserialize(input_stream);
}

void RigidBody2DSim::computeContactPoints(std::vector<Vector2s> &points,
                                          std::vector<Vector2s> &normals) {
  points.clear();
  normals.clear();

  std::vector<std::unique_ptr<Constraint>> active_set;
  computeActiveSet(m_state.q(), m_state.q(), m_state.v(), false, active_set);

  for (const std::unique_ptr<Constraint> &con : active_set) {
    VectorXs point;
    con->getWorldSpaceContactPoint(m_state.q(), point);
    points.emplace_back(point);
    VectorXs normal;
    con->getWorldSpaceContactNormal(m_state.q(), normal);
    normals.emplace_back(normal);
  }
}

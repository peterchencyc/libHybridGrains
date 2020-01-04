#include "PenaltyImpactFrictionMap.h"

#include <iostream>

#include "scisim/Constraints/ConstrainedSystem.h"
#include "scisim/Constraints/Constraint.h"
#include "scisim/Math/MathUtilities.h"
#include "scisim/UnconstrainedMaps/FlowableSystem.h"

#include "StaticDrumCircleConstraint.h"
#include "StaticPlaneCircleConstraint.h"

#include "scisim/Timer/TimeUtils.h"

PenaltyImpactFrictionMap::PenaltyImpactFrictionMap(const scalar &kn,
                                                   const scalar &kt,
                                                   const scalar &gamman,
                                                   const scalar &gammat,
                                                   const scalar &mu)
    : m_kn(kn), m_kt(kt), m_gamman(gamman), m_gammat(gammat), m_mu(mu) {
  assert(m_kn >= 0.0);
  assert(m_kt >= 0.0);
  assert(m_gamman >= 0.0);
  assert(m_gammat >= 0.0);
  assert(m_mu >= 0.0);
}

void PenaltyImpactFrictionMap::updateCollisionCache(const bool reduce_bandwidth,
                                                    const scalar &dt,
                                                    ConstrainedSystem &csys,
                                                    const VectorXs &q,
                                                    const VectorXs &v) {
  assert(q.size() % 3 == 0);
  const int nbodies = q.size() / 3;

  // CollisionCache new_cache;

  // assert( m_old_cache.empty() );
  m_old_cache.clear();
  m_old_cache.resizeNumBodies(nbodies);

  std::vector<std::unique_ptr<Constraint>> active_set;
  csys.computeActiveSet(q, q, v, reduce_bandwidth, active_set);
  for (const std::unique_ptr<Constraint> &con : active_set) {
    if (con->name() == "static_plane_circle") {
      const StaticPlaneCircleConstraint &spc_con =
          *static_cast<StaticPlaneCircleConstraint *>(con.get());

      std::pair<int, int> indices;
      con->getBodyIndices(indices);
      indices.second = con->getStaticObjectIndex();
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const Vector2s r = -spc_con.circleRadius() * n;
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      CircleHalfPlaneCollision *result =
          m_collision_cache.circle_half_plane_collisions.find(indices.first,
                                                              indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_half_plane_collisions.insert(
            indices.first, indices.second,
            CircleHalfPlaneCollision(depth, n, delta_s, vrel, r, 1.0));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_half_plane_collisions.insert(
            indices.first, indices.second,
            CircleHalfPlaneCollision(depth, n, delta_s, vrel, r, 1.0));
      }
    } else if (con->name() == "static_drum_circle") {
      const StaticDrumCircleConstraint &spc_con =
          *static_cast<StaticDrumCircleConstraint *>(con.get());

      std::pair<int, int> indices;
      con->getBodyIndices(indices);
      indices.second = con->getStaticObjectIndex();
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const Vector2s r = -spc_con.circleRadius() * n;
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      CircleDrumCollision *result =
          m_collision_cache.circle_drum_collisions.find(indices.first,
                                                        indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_drum_collisions.insert(
            indices.first, indices.second,
            CircleDrumCollision(depth, n, delta_s, vrel, r, 1.0));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_drum_collisions.insert(
            indices.first, indices.second,
            CircleDrumCollision(depth, n, delta_s, vrel, r, 1.0));
      }
    } else if (con->name() == "circle_circle") {
      // const StaticPlaneCircleConstraint& spc_con =
      // *static_cast<StaticPlaneCircleConstraint*>( con.get() );
      VectorXs contact_point;
      con->getWorldSpaceContactPoint(q, contact_point);
      assert(contact_point.size() == 2);

      std::pair<int, int> indices;
      con->getBodyIndices(indices);
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 2);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      CircleCircleCollision *result =
          m_collision_cache.circle_circle_collisions.find(indices.first,
                                                          indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_circle_collisions.insert(
            indices.first, indices.second,
            CircleCircleCollision(depth, n, delta_s, vrel, r.first, r.second,
                                  1.0, contact_point));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_circle_collisions.insert(
            indices.first, indices.second,
            CircleCircleCollision(depth, n, delta_s, vrel, r.first, r.second,
                                  1.0, contact_point)); //, delta_s, vrel, r
      }
    } else if (con->name() == "kinematic_circle_circle") {
      std::pair<int, int> indices;
      con->getBodyIndices(indices);
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 0);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      KinematicCircleCircleCollision *result =
          m_collision_cache.kinematic_circle_circle_collisions.find(
              indices.first, indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_circle_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                           1.0));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_circle_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                           1.0));
      }
    } else if (con->name() == "teleported_circle_circle") {
      std::pair<int, int> indices;
      con->getBodyIndices(indices);
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 2);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      TeleportedCircleCircleCollision *result =
          m_collision_cache.teleported_circle_circle_collisions.find(
              indices.first, indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.teleported_circle_circle_collisions.insert(
            indices.first, indices.second,
            TeleportedCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                            r.second, 1.0));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.teleported_circle_circle_collisions.insert(
            indices.first, indices.second,
            TeleportedCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                            r.second, 1.0));
      }
    } else if (con->name() == "kinematic_box_circle") {
      std::pair<int, int> indices;
      con->getBodyIndices(indices);
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 0);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      KinematicCircleBoxCollision *result =
          m_collision_cache.kinematic_box_circle_collisions.find(
              indices.first, indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_box_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleBoxCollision(depth, n, delta_s, vrel, r.first, 1.0));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_box_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleBoxCollision(depth, n, delta_s, vrel, r.first, 1.0));
      }
    } else {
      std::cerr
          << "Error, unsupported constraint type in penalty collisions 1: "
          << con->name() << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::swap(m_collision_cache, m_old_cache);

  // m_old_cache.clear();
}

void PenaltyImpactFrictionMap::updateCollisionCache(
    const bool reduce_bandwidth, const scalar &dt, ConstrainedSystem &csys,
    const VectorXs &q, const VectorXs &v, const VectorXs &w) {
  assert(q.size() % 3 == 0);
  const int nbodies = q.size() / 3;

  // CollisionCache new_cache;

  // assert( m_old_cache.empty() );
  m_old_cache.clear();
  m_old_cache.resizeNumBodies(nbodies);

  TimingTools timing_tools;
  timing_tools.start();
  std::vector<std::unique_ptr<Constraint>> active_set;
  csys.computeActiveSet(q, q, v, reduce_bandwidth, active_set);
  timing_tools.stop("");
  timing_tools.start();
  for (const std::unique_ptr<Constraint> &con : active_set) {
    if (con->name() == "static_plane_circle") {
      const StaticPlaneCircleConstraint &spc_con =
          *static_cast<StaticPlaneCircleConstraint *>(con.get());

      VectorXs contact_point;
      con->getWorldSpaceContactPoint(q, contact_point);
      assert(contact_point.size() == 2);

      std::pair<int, int> indices;
      con->getBodyIndices(indices);

      const scalar weight = w(indices.first);

      indices.second = con->getStaticObjectIndex();
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const Vector2s r = -spc_con.circleRadius() * n;
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      CircleHalfPlaneCollision *result =
          m_collision_cache.circle_half_plane_collisions.find(indices.first,
                                                              indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_half_plane_collisions.insert(
            indices.first, indices.second,
            CircleHalfPlaneCollision(depth, n, delta_s, vrel, r, weight));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_half_plane_collisions.insert(
            indices.first, indices.second,
            CircleHalfPlaneCollision(depth, n, delta_s, vrel, r, weight));
      }
    } else if (con->name() == "static_drum_circle") {
      const StaticDrumCircleConstraint &spc_con =
          *static_cast<StaticDrumCircleConstraint *>(con.get());

      std::pair<int, int> indices;
      con->getBodyIndices(indices);

      const scalar weight = w(indices.first);

      indices.second = con->getStaticObjectIndex();
      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const Vector2s r = -spc_con.circleRadius() * n;
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      CircleDrumCollision *result =
          m_collision_cache.circle_drum_collisions.find(indices.first,
                                                        indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_drum_collisions.insert(
            indices.first, indices.second,
            CircleDrumCollision(depth, n, delta_s, vrel, r, weight));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_drum_collisions.insert(
            indices.first, indices.second,
            CircleDrumCollision(depth, n, delta_s, vrel, r, weight));
      }
    } else if (con->name() == "circle_circle") {
      // const StaticPlaneCircleConstraint& spc_con =
      // *static_cast<StaticPlaneCircleConstraint*>( con.get() );

      VectorXs contact_point;
      con->getWorldSpaceContactPoint(q, contact_point);
      assert(contact_point.size() == 2);

      std::pair<int, int> indices;
      con->getBodyIndices(indices);

      const scalar weight1 = w(indices.first);
      const scalar weight2 = w(indices.second);
      const scalar weight = (weight1 + weight2) * 0.5;

      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 2);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      CircleCircleCollision *result =
          m_collision_cache.circle_circle_collisions.find(indices.first,
                                                          indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_circle_collisions.insert(
            indices.first, indices.second,
            CircleCircleCollision(depth, n, delta_s, vrel, r.first, r.second,
                                  weight, contact_point));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.circle_circle_collisions.insert(
            indices.first, indices.second,
            CircleCircleCollision(depth, n, delta_s, vrel, r.first, r.second,
                                  weight, contact_point)); //, delta_s, vrel, r
      }
    } else if (con->name() == "kinematic_circle_circle") {
      VectorXs contact_point;
      con->getWorldSpaceContactPoint(q, contact_point);
      assert(contact_point.size() == 2);

      std::pair<int, int> indices;
      con->getBodyIndices(indices);

      const scalar weight = w(indices.first);

      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 0);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      KinematicCircleCircleCollision *result =
          m_collision_cache.kinematic_circle_circle_collisions.find(
              indices.first, indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_circle_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                           weight));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_circle_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                           weight));
      }
    } else if (con->name() == "teleported_circle_circle") {
      std::pair<int, int> indices;
      con->getBodyIndices(indices);

      const scalar weight = w(indices.first);

      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 2);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      TeleportedCircleCircleCollision *result =
          m_collision_cache.teleported_circle_circle_collisions.find(
              indices.first, indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.teleported_circle_circle_collisions.insert(
            indices.first, indices.second,
            TeleportedCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                            r.second, weight));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.teleported_circle_circle_collisions.insert(
            indices.first, indices.second,
            TeleportedCircleCircleCollision(depth, n, delta_s, vrel, r.first,
                                            r.second, weight));
      }
    } else if (con->name() == "kinematic_box_circle") {
      VectorXs contact_point;
      con->getWorldSpaceContactPoint(q, contact_point);
      assert(contact_point.size() == 2);

      std::pair<int, int> indices;
      con->getBodyIndices(indices);

      const scalar weight = w(indices.first);

      // TODO: Check that the indices are ok
      const scalar depth = con->penetrationDepth(q);
      assert(depth <= 0);
      VectorXs n;
      con->getWorldSpaceContactNormal(q, n);
      assert(n.size() == 2);
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      // Compute the lever arm
      const std::pair<VectorXs, VectorXs> r = con->computeLeverArms(q);
      assert(r.first.size() == 2);
      assert(r.second.size() == 0);
      // Compute the relative velocity at the contact point
      const Vector2s vrel = con->computeRelativeVelocity(q, v);
      // Check if this collision existed in the last time-step
      KinematicCircleBoxCollision *result =
          m_collision_cache.kinematic_box_circle_collisions.find(
              indices.first, indices.second);
      if (result == nullptr) {
        Vector2s delta_s = dt * vrel;
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_box_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleBoxCollision(depth, n, delta_s, vrel, r.first,
                                        weight));
      } else {
        // Get the old delta s
        Vector2s delta_s = result->delta_s;
        // Add the current displacement to the total delta s
        delta_s += dt * vrel;
        // Project out the component along n
        delta_s -= n.dot(delta_s) * n;
        m_old_cache.kinematic_box_circle_collisions.insert(
            indices.first, indices.second,
            KinematicCircleBoxCollision(depth, n, delta_s, vrel, r.first,
                                        weight));
      }
    } else {
      std::cerr
          << "Error, unsupported constraint type in penalty collisions 0: "
          << con->name() << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  timing_tools.stop("");
  std::swap(m_collision_cache, m_old_cache);

  // m_old_cache.clear();
}

static scalar secondRootOfQuadratic(const scalar &a, const scalar &b,
                                    const scalar &c, const scalar &dscr_sqrt) {
  scalar root;
  if (b > 0.0) {
    assert((-b - dscr_sqrt) != 0.0);
    root = (2.0 * c) / (-b - dscr_sqrt);
    assert(root == root);
  } else {
    assert(a != 0.0);
    root = (-b + dscr_sqrt) / (2.0 * a);
    assert(root == root);
  }
  return root;
}

void PenaltyImpactFrictionMap::updateAnchorsAndAccumulateCollisionForces(
    const VectorXs &v, VectorXs &F) {
  assert(v.size() % 3 == 0);
  const int nbodies = v.size() / 3;

  // Circle vs plane collisions
  for (int bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    for (std::pair<int, CircleHalfPlaneCollision> &circle_vs_plane :
         m_collision_cache.circle_half_plane_collisions.collisions()[bdy_idx]) {
      const scalar w{circle_vs_plane.second.weight};

      // Cache out some useful state for this collision
      const Vector2s &n{circle_vs_plane.second.n};
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      const Vector2s &deltas{circle_vs_plane.second.delta_s};
      assert(std::fabs(n.dot(deltas)) <= 1.0e-6);
      const scalar &pen_depth{circle_vs_plane.second.pen_depth};
      assert(pen_depth <= 0.0);

      // Break the relative velocity into normal and tangential components
      const Vector2s &vel{circle_vs_plane.second.vrel};
      const Vector2s vt{vel - n.dot(vel) * n};
      const Vector2s vn{vel - vt};

      // Compute the total normal force
      const Vector2s normal_force{-m_kn * w * pen_depth * n -
                                  0.5 * m_gamman * w * vn};

      // Tangential friction
      Vector2s friction_force{-m_kt * w * deltas - 0.5 * m_gammat * w * vt};
      // Project to respect Coulomb friction
      {
        const scalar mu_fn = m_mu * normal_force.norm();
        const scalar ft = friction_force.norm();
        // If modifying delta_s so the force lies within the cone is impossible,
        // zero out delta_s and set the force to the maximum possible value
        // opposing the slip direction
        if (0.5 * m_gammat * w * vt.norm() > mu_fn) {
          circle_vs_plane.second.delta_s.setZero();
          friction_force = -0.5 * m_gammat * w * vt;
          friction_force.normalize();
          friction_force *= mu_fn;
        }
        // Otherwise, shrink delta_s so the force lies within the cone
        else if (ft > mu_fn) {
          // Elastic and with damping case
          const scalar a{m_kt * w * m_kt * w * deltas.squaredNorm()};
          assert(a >= 0.0);
          const scalar b{m_kt * w * m_gammat * w * deltas.dot(vt)};
          const scalar c{0.25 * m_gammat * w * m_gammat * w * vt.dot(vt) -
                         m_mu * m_mu * normal_force.squaredNorm()};
          const scalar dscr{b * b - 4 * a * c};
          assert(dscr >= 0.0);
          const scalar dscr_sqrt{std::sqrt(dscr)};
          // const scalar root0 = firstRootOfQuadratic( a, b, c, dscr_sqrt );
          const scalar root1{secondRootOfQuadratic(a, b, c, dscr_sqrt)};
          assert(root1 >= 0.0);
          assert(root1 <= 1.0);
          assert(std::fabs(a * root1 * root1 + b * root1 + c) <= 1.0e-5);

          const Vector2s new_delta_s{root1 * deltas};
          const Vector2s new_friction_force{-m_kt * w * new_delta_s -
                                            0.5 * m_gammat * w * vt};
          assert(std::fabs(new_friction_force.norm() - mu_fn) <= 1.0e-5);
          friction_force = new_friction_force;
          circle_vs_plane.second.delta_s = new_delta_s;
        }
      }
      assert(friction_force.norm() <= m_mu * normal_force.norm() + 1.0e-5);

      const Vector2s total_force = normal_force + friction_force;
      F.segment<2>(3 * bdy_idx) += total_force;
      F(3 * bdy_idx + 2) +=
          MathUtilities::cross(circle_vs_plane.second.r, total_force);

      circle_vs_plane.second.normal_force = normal_force;
      circle_vs_plane.second.friction_force = friction_force;
    }
  }

  // Circle vs drum collisions
  for (int bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    for (std::pair<int, CircleDrumCollision> &circle_vs_drum :
         m_collision_cache.circle_drum_collisions.collisions()[bdy_idx]) {
      const scalar w{circle_vs_drum.second.weight};

      // Cache out some useful state for this collision
      const Vector2s &n{circle_vs_drum.second.n};
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      const Vector2s &deltas{circle_vs_drum.second.delta_s};
      assert(std::fabs(n.dot(deltas)) <= 1.0e-6);
      const scalar &pen_depth{circle_vs_drum.second.pen_depth};
      assert(pen_depth <= 0.0);

      // Break the relative velocity into normal and tangential components
      const Vector2s &vel{circle_vs_drum.second.vrel};
      const Vector2s vt{vel - n.dot(vel) * n};
      const Vector2s vn{vel - vt};

      // Compute the total normal force
      const Vector2s normal_force{-m_kn * w * pen_depth * n -
                                  0.5 * m_gamman * w * vn};

      // Tangential friction
      Vector2s friction_force{-m_kt * w * deltas - 0.5 * m_gammat * w * vt};
      // Project to respect Coulomb friction
      {
        const scalar mu_fn = m_mu * normal_force.norm();
        const scalar ft = friction_force.norm();
        // If modifying delta_s so the force lies within the cone is impossible,
        // zero out delta_s and set the force to the maximum possible value
        // opposing the slip direction
        if (0.5 * m_gammat * w * vt.norm() > mu_fn) {
          circle_vs_drum.second.delta_s.setZero();
          friction_force = -0.5 * m_gammat * w * vt;
          friction_force.normalize();
          friction_force *= mu_fn;
        }
        // Otherwise, shrink delta_s so the force lies within the cone
        else if (ft > mu_fn) {
          // Elastic and with damping case
          const scalar a{m_kt * w * m_kt * w * deltas.squaredNorm()};
          assert(a >= 0.0);
          const scalar b{m_kt * w * m_gammat * w * deltas.dot(vt)};
          const scalar c{0.25 * m_gammat * w * m_gammat * w * vt.dot(vt) -
                         m_mu * m_mu * normal_force.squaredNorm()};
          const scalar dscr{b * b - 4 * a * c};
          assert(dscr >= 0.0);
          const scalar dscr_sqrt{std::sqrt(dscr)};
          // const scalar root0 = firstRootOfQuadratic( a, b, c, dscr_sqrt );
          const scalar root1{secondRootOfQuadratic(a, b, c, dscr_sqrt)};
          assert(root1 >= 0.0);
          assert(root1 <= 1.0);
          assert(std::fabs(a * root1 * root1 + b * root1 + c) <= 1.0e-5);

          const Vector2s new_delta_s{root1 * deltas};
          const Vector2s new_friction_force{-m_kt * w * new_delta_s -
                                            0.5 * m_gammat * w * vt};
          assert(std::fabs(new_friction_force.norm() - mu_fn) <= 1.0e-5);
          friction_force = new_friction_force;
          circle_vs_drum.second.delta_s = new_delta_s;
        }
      }
      assert(friction_force.norm() <= m_mu * normal_force.norm() + 1.0e-5);

      const Vector2s total_force = normal_force + friction_force;
      F.segment<2>(3 * bdy_idx) += total_force;
      F(3 * bdy_idx + 2) +=
          MathUtilities::cross(circle_vs_drum.second.r, total_force);

      circle_vs_drum.second.normal_force = normal_force;
      circle_vs_drum.second.friction_force = friction_force;
    }
  }

  // Circle vs circle collisions
  for (int bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    for (std::pair<int, CircleCircleCollision> &circle_vs_circle :
         m_collision_cache.circle_circle_collisions.collisions()[bdy_idx]) {
      const scalar w{circle_vs_circle.second.weight};

      // Cache out some useful state for this collision
      const Vector2s &n{circle_vs_circle.second.n};
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      const Vector2s &deltas{circle_vs_circle.second.delta_s};
      assert(std::fabs(n.dot(deltas)) <= 1.0e-6);
      const scalar &pen_depth{circle_vs_circle.second.pen_depth};
      assert(pen_depth <= 0.0);

      // Break the relative velocity into normal and tangential components
      const Vector2s &vel{circle_vs_circle.second.vrel};
      const Vector2s vt{vel - n.dot(vel) * n};
      const Vector2s vn{vel - vt};

      // Compute the total normal force
      const Vector2s normal_force{-m_kn * w * pen_depth * n -
                                  0.5 * m_gamman * w * vn};

      // Tangential friction
      Vector2s friction_force{-m_kt * w * deltas - 0.5 * m_gammat * w * vt};
      // Project to respect Coulomb friction
      {
        const scalar mu_fn = m_mu * normal_force.norm();
        const scalar ft = friction_force.norm();
        // If modifying delta_s so the force lies within the cone is impossible,
        // zero out delta_s and set the force to the maximum possible value
        // opposing the slip direction
        if (0.5 * m_gammat * w * vt.norm() > mu_fn) {
          circle_vs_circle.second.delta_s.setZero();
          friction_force = -0.5 * m_gammat * w * vt;
          friction_force.normalize();
          friction_force *= mu_fn;
        }
        // Otherwise, shrink delta_s so the force lies within the cone
        else if (ft > mu_fn) {
          // Elastic and with damping case
          const scalar a{m_kt * w * m_kt * w * deltas.squaredNorm()};
          assert(a >= 0.0);
          const scalar b{m_kt * w * m_gammat * w * deltas.dot(vt)};
          const scalar c{0.25 * m_gammat * w * m_gammat * w * vt.dot(vt) -
                         m_mu * m_mu * normal_force.squaredNorm()};
          const scalar dscr{b * b - 4 * a * c};
          assert(dscr >= 0.0);
          const scalar dscr_sqrt{std::sqrt(dscr)};
          // const scalar root0 = firstRootOfQuadratic( a, b, c, dscr_sqrt );
          const scalar root1{secondRootOfQuadratic(a, b, c, dscr_sqrt)};
          assert(root1 >= 0.0);
          assert(root1 <= 1.0);
          assert(std::fabs(a * root1 * root1 + b * root1 + c) <= 1.0e-5);

          const Vector2s new_delta_s{root1 * deltas};
          const Vector2s new_friction_force{-m_kt * w * new_delta_s -
                                            0.5 * m_gammat * w * vt};
          assert(std::fabs(new_friction_force.norm() - mu_fn) <= 1.0e-5);
          friction_force = new_friction_force;
          circle_vs_circle.second.delta_s = new_delta_s;
        }
      }
      assert(friction_force.norm() <= m_mu * normal_force.norm() + 1.0e-5);

      const Vector2s total_force = normal_force + friction_force;
      F.segment<2>(3 * bdy_idx) += total_force;
      F(3 * bdy_idx + 2) +=
          MathUtilities::cross(circle_vs_circle.second.r0, total_force);
      F.segment<2>(3 * circle_vs_circle.first) -= total_force;
      F(3 * circle_vs_circle.first + 2) -=
          MathUtilities::cross(circle_vs_circle.second.r1, total_force);

      circle_vs_circle.second.normal_force = normal_force;
      circle_vs_circle.second.friction_force = friction_force;
    }
  }

  // Teleported circle vs circle collisions
  for (int bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    for (std::pair<int, TeleportedCircleCircleCollision>
             &teleported_circle_vs_circle :
         m_collision_cache.teleported_circle_circle_collisions
             .collisions()[bdy_idx]) {
      const scalar w{teleported_circle_vs_circle.second.weight};

      // Cache out some useful state for this collision
      const Vector2s &n{teleported_circle_vs_circle.second.n};
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      const Vector2s &deltas{teleported_circle_vs_circle.second.delta_s};
      assert(std::fabs(n.dot(deltas)) <= 1.0e-6);
      const scalar &pen_depth{teleported_circle_vs_circle.second.pen_depth};
      assert(pen_depth <= 0.0);

      // Break the relative velocity into normal and tangential components
      const Vector2s &vel{teleported_circle_vs_circle.second.vrel};
      const Vector2s vt{vel - n.dot(vel) * n};
      const Vector2s vn{vel - vt};

      // Compute the total normal force
      const Vector2s normal_force{-m_kn * w * pen_depth * n -
                                  0.5 * m_gamman * w * vn};

      // Tangential friction
      Vector2s friction_force{-m_kt * w * deltas - 0.5 * m_gammat * w * vt};
      // Project to respect Coulomb friction
      {
        const scalar mu_fn = m_mu * normal_force.norm();
        const scalar ft = friction_force.norm();
        // If modifying delta_s so the force lies within the cone is impossible,
        // zero out delta_s and set the force to the maximum possible value
        // opposing the slip direction
        if (0.5 * m_gammat * w * vt.norm() > mu_fn) {
          teleported_circle_vs_circle.second.delta_s.setZero();
          friction_force = -0.5 * m_gammat * w * vt;
          friction_force.normalize();
          friction_force *= mu_fn;
        }
        // Otherwise, shrink delta_s so the force lies within the cone
        else if (ft > mu_fn) {
          // Elastic and with damping case
          const scalar a{m_kt * w * m_kt * w * deltas.squaredNorm()};
          assert(a >= 0.0);
          const scalar b{m_kt * w * m_gammat * w * deltas.dot(vt)};
          const scalar c{0.25 * m_gammat * w * m_gammat * w * vt.dot(vt) -
                         m_mu * m_mu * normal_force.squaredNorm()};
          const scalar dscr{b * b - 4 * a * c};
          assert(dscr >= 0.0);
          const scalar dscr_sqrt{std::sqrt(dscr)};
          // const scalar root0 = firstRootOfQuadratic( a, b, c, dscr_sqrt );
          const scalar root1{secondRootOfQuadratic(a, b, c, dscr_sqrt)};
          assert(root1 >= 0.0);
          assert(root1 <= 1.0);
          assert(std::fabs(a * root1 * root1 + b * root1 + c) <= 1.0e-5);

          const Vector2s new_delta_s{root1 * deltas};
          const Vector2s new_friction_force{-m_kt * w * new_delta_s -
                                            0.5 * m_gammat * w * vt};
          assert(std::fabs(new_friction_force.norm() - mu_fn) <= 1.0e-5);
          friction_force = new_friction_force;
          teleported_circle_vs_circle.second.delta_s = new_delta_s;
        }
      }
      assert(friction_force.norm() <= m_mu * normal_force.norm() + 1.0e-5);

      const Vector2s total_force = normal_force + friction_force;
      F.segment<2>(3 * bdy_idx) += total_force;
      F(3 * bdy_idx + 2) += MathUtilities::cross(
          teleported_circle_vs_circle.second.r0, total_force);
      F.segment<2>(3 * teleported_circle_vs_circle.first) -= total_force;
      F(3 * teleported_circle_vs_circle.first + 2) -= MathUtilities::cross(
          teleported_circle_vs_circle.second.r1, total_force);

      teleported_circle_vs_circle.second.normal_force = normal_force;
      teleported_circle_vs_circle.second.friction_force = friction_force;
    }
  }

  // Kinematic circle vs circle collisions
  for (int bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    for (std::pair<int, KinematicCircleCircleCollision>
             &kinematic_circle_vs_circle :
         m_collision_cache.kinematic_circle_circle_collisions
             .collisions()[bdy_idx]) {
      const scalar w{kinematic_circle_vs_circle.second.weight};

      // Cache out some useful state for this collision
      const Vector2s &n{kinematic_circle_vs_circle.second.n};
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      const Vector2s &deltas{kinematic_circle_vs_circle.second.delta_s};
      assert(std::fabs(n.dot(deltas)) <= 1.0e-6);
      const scalar &pen_depth{kinematic_circle_vs_circle.second.pen_depth};
      assert(pen_depth <= 0.0);

      // Break the relative velocity into normal and tangential components
      const Vector2s &vel{kinematic_circle_vs_circle.second.vrel};
      const Vector2s vt{vel - n.dot(vel) * n};
      const Vector2s vn{vel - vt};

      // Compute the total normal force
      const Vector2s normal_force{-m_kn * w * pen_depth * n -
                                  0.5 * m_gamman * w * vn};

      // Tangential friction
      Vector2s friction_force{-m_kt * w * deltas - 0.5 * m_gammat * w * vt};
      // Project to respect Coulomb friction
      {
        const scalar mu_fn = m_mu * normal_force.norm();
        const scalar ft = friction_force.norm();
        // If modifying delta_s so the force lies within the cone is impossible,
        // zero out delta_s and set the force to the maximum possible value
        // opposing the slip direction
        if (0.5 * m_gammat * w * vt.norm() > mu_fn) {
          kinematic_circle_vs_circle.second.delta_s.setZero();
          friction_force = -0.5 * m_gammat * w * vt;
          friction_force.normalize();
          friction_force *= mu_fn;
        }
        // Otherwise, shrink delta_s so the force lies within the cone
        else if (ft > mu_fn) {
          // Elastic and with damping case
          const scalar a{m_kt * w * m_kt * w * deltas.squaredNorm()};
          assert(a >= 0.0);
          const scalar b{m_kt * w * m_gammat * w * deltas.dot(vt)};
          const scalar c{0.25 * m_gammat * w * m_gammat * w * vt.dot(vt) -
                         m_mu * m_mu * normal_force.squaredNorm()};
          const scalar dscr{b * b - 4 * a * c};
          assert(dscr >= 0.0);
          const scalar dscr_sqrt{std::sqrt(dscr)};
          // const scalar root0 = firstRootOfQuadratic( a, b, c, dscr_sqrt );
          const scalar root1{secondRootOfQuadratic(a, b, c, dscr_sqrt)};
          assert(root1 >= 0.0);
          assert(root1 <= 1.0);
          assert(std::fabs(a * root1 * root1 + b * root1 + c) <= 1.0e-5);

          const Vector2s new_delta_s{root1 * deltas};
          const Vector2s new_friction_force{-m_kt * w * new_delta_s -
                                            0.5 * m_gammat * w * vt};
          assert(std::fabs(new_friction_force.norm() - mu_fn) <= 1.0e-5);
          friction_force = new_friction_force;
          kinematic_circle_vs_circle.second.delta_s = new_delta_s;
        }
      }
      assert(friction_force.norm() <= m_mu * normal_force.norm() + 1.0e-5);

      const Vector2s total_force = normal_force + friction_force;
      F.segment<2>(3 * bdy_idx) += total_force;
      F(3 * bdy_idx + 2) += MathUtilities::cross(
          kinematic_circle_vs_circle.second.r0, total_force);

      kinematic_circle_vs_circle.second.normal_force = normal_force;
      kinematic_circle_vs_circle.second.friction_force = friction_force;
    }
  }

  // Kinematic box vs circle collisions
  for (int bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    for (std::pair<int, KinematicCircleBoxCollision> &kinematic_box_vs_circle :
         m_collision_cache.kinematic_box_circle_collisions
             .collisions()[bdy_idx]) {
      const scalar w{kinematic_box_vs_circle.second.weight};

      // Cache out some useful state for this collision
      const Vector2s &n{kinematic_box_vs_circle.second.n};
      assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
      const Vector2s &deltas{kinematic_box_vs_circle.second.delta_s};
      assert(std::fabs(n.dot(deltas)) <= 1.0e-6);
      const scalar &pen_depth{kinematic_box_vs_circle.second.pen_depth};
      assert(pen_depth <= 0.0);

      // Break the relative velocity into normal and tangential components
      const Vector2s &vel{kinematic_box_vs_circle.second.vrel};
      const Vector2s vt{vel - n.dot(vel) * n};
      const Vector2s vn{vel - vt};

      // Compute the total normal force
      const Vector2s normal_force{-m_kn * w * pen_depth * n -
                                  0.5 * m_gamman * w * vn};

      // Tangential friction
      Vector2s friction_force{-m_kt * w * deltas - 0.5 * m_gammat * w * vt};
      // Project to respect Coulomb friction
      {
        const scalar mu_fn = m_mu * normal_force.norm();
        const scalar ft = friction_force.norm();
        // If modifying delta_s so the force lies within the cone is impossible,
        // zero out delta_s and set the force to the maximum possible value
        // opposing the slip direction
        if (0.5 * m_gammat * w * vt.norm() > mu_fn) {
          kinematic_box_vs_circle.second.delta_s.setZero();
          friction_force = -0.5 * m_gammat * w * vt;
          friction_force.normalize();
          friction_force *= mu_fn;
        }
        // Otherwise, shrink delta_s so the force lies within the cone
        else if (ft > mu_fn) {
          // Elastic and with damping case
          const scalar a{m_kt * w * m_kt * w * deltas.squaredNorm()};
          assert(a >= 0.0);
          const scalar b{m_kt * w * m_gammat * w * deltas.dot(vt)};
          const scalar c{0.25 * m_gammat * w * m_gammat * w * vt.dot(vt) -
                         m_mu * m_mu * normal_force.squaredNorm()};
          const scalar dscr{b * b - 4 * a * c};
          assert(dscr >= 0.0);
          const scalar dscr_sqrt{std::sqrt(dscr)};
          // const scalar root0 = firstRootOfQuadratic( a, b, c, dscr_sqrt );
          const scalar root1{secondRootOfQuadratic(a, b, c, dscr_sqrt)};
          assert(root1 >= 0.0);
          assert(root1 <= 1.0);
          assert(std::fabs(a * root1 * root1 + b * root1 + c) <= 1.0e-5);

          const Vector2s new_delta_s{root1 * deltas};
          const Vector2s new_friction_force{-m_kt * w * new_delta_s -
                                            0.5 * m_gammat * w * vt};
          assert(std::fabs(new_friction_force.norm() - mu_fn) <= 1.0e-5);
          friction_force = new_friction_force;
          kinematic_box_vs_circle.second.delta_s = new_delta_s;
        }
      }
      assert(friction_force.norm() <= m_mu * normal_force.norm() + 1.0e-5);

      const Vector2s total_force = normal_force + friction_force;
      F.segment<2>(3 * bdy_idx) += total_force;
      F(3 * bdy_idx + 2) +=
          MathUtilities::cross(kinematic_box_vs_circle.second.r0, total_force);

      kinematic_box_vs_circle.second.normal_force = normal_force;
      kinematic_box_vs_circle.second.friction_force = friction_force;
    }
  }
}

void PenaltyImpactFrictionMap::flow(
    ScriptingCallback & /*call_back*/, FlowableSystem &fsys,
    ConstrainedSystem &csys, UnconstrainedMap & /*umap*/,
    FrictionSolver & /*friction_solver*/, const unsigned /*iteration*/,
    const scalar &dt, const scalar & /*CoR_default*/,
    const scalar & /*mu_default*/, const bool reduce_bandwidth,
    const VectorXs &q0, const VectorXs &v0, VectorXs &q1, VectorXs &v1) {
  assert(dt > 0.0);
  assert(q0.size() == fsys.nqdofs());
  assert(q1.size() == fsys.nqdofs());
  assert(v0.size() == fsys.nvdofs());
  assert(v1.size() == fsys.nvdofs());
  assert(fsys.nqdofs() % 3 == 0);

  const int nbodies = fsys.nqdofs() / 3;

  updateCollisionCache(reduce_bandwidth, dt, csys, q0, v0);

  VectorXs F(fsys.nqdofs());

  fsys.computeForce(q0, v0, dt, F);

  updateAnchorsAndAccumulateCollisionForces(v0, F);

  for (int bdyIdx = 0; bdyIdx < nbodies; bdyIdx++) {
    if (fsys.isKinematicallyScripted(bdyIdx)) {
      F.segment<3>(3 * bdyIdx).setZero();
    }
  }

  v1 = v0 + dt * fsys.Minv() * F;
  q1 = q0 + dt * v1;
}

void PenaltyImpactFrictionMap::flowWithWeights(
    ScriptingCallback &call_back, FlowableSystem &fsys, ConstrainedSystem &csys,
    UnconstrainedMap &umap, FrictionSolver &friction_solver,
    const unsigned iteration, const scalar &dt, const scalar &CoR,
    const scalar &mu_default, const bool reduce_bandwidth, const VectorXs &q0,
    const VectorXs &v0, const VectorXs &w, VectorXs &q1, VectorXs &v1) {
  assert(dt > 0.0);
  assert(q0.size() == fsys.nqdofs());
  assert(q1.size() == fsys.nqdofs());
  assert(v0.size() == fsys.nvdofs());
  assert(v1.size() == fsys.nvdofs());
  assert(fsys.nqdofs() % 3 == 0);

  const int nbodies = fsys.nqdofs() / 3;
  TimingTools timing_tools;
  timing_tools.start();
  updateCollisionCache(reduce_bandwidth, dt, csys, q0, v0, w);
  timing_tools.stop("");
  VectorXs F(fsys.nqdofs());

  fsys.computeForce(q0, v0, dt, F);

  updateAnchorsAndAccumulateCollisionForces(v0, F);

  for (int bdyIdx = 0; bdyIdx < nbodies; bdyIdx++) {
    if (fsys.isKinematicallyScripted(bdyIdx)) {
      F.segment<3>(3 * bdyIdx).setZero();
    }
  }

  v1 = v0 + dt * fsys.Minv() * F;
  q1 = q0 + dt * v1;
}

bool PenaltyImpactFrictionMap::stabilize() const { return false; }

void PenaltyImpactFrictionMap::resetCachedData() { m_collision_cache.clear(); }

void PenaltyImpactFrictionMap::serialize(std::ostream &output_stream) const {
  std::cerr << "PenaltyImpactFrictionMap::serialize not coded up." << std::endl;
  std::exit(EXIT_FAILURE);
}

std::string PenaltyImpactFrictionMap::name() const {
  return "penalty_impact_friction_map";
}

void PenaltyImpactFrictionMap::exportForcesNextStep(HDF5File &output_file) {
  std::cerr << "PenaltyImpactFrictionMap::exportForcesNextStep not coded up."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

void PenaltyImpactFrictionMap::saveImpulsesNextStep(
    CollisionImpulses &impulses) {
  std::cerr << "PenaltyImpactFrictionMap::saveImpulsesNextStep not coded up."
            << std::endl;
  std::exit(EXIT_FAILURE);
}

std::unique_ptr<ImpactFrictionMap> PenaltyImpactFrictionMap::clone() const {
  return std::unique_ptr<ImpactFrictionMap>{
      new PenaltyImpactFrictionMap{m_kn, m_kt, m_gamman, m_gammat, m_mu}};
}

void PenaltyImpactFrictionMap::enlargeCache(const int new_size) {
#ifndef NDEBUG
  const int original_collision_count_0 = m_collision_cache.numTotalCollisions();
  const int original_collision_count_1 = m_old_cache.numTotalCollisions();
#endif

  m_collision_cache.resizeNumBodies(new_size);
  m_old_cache.resizeNumBodies(new_size);

#ifndef NDEBUG
  const int new_collision_count_0 = m_collision_cache.numTotalCollisions();
  const int new_collision_count_1 = m_old_cache.numTotalCollisions();
#endif
  assert(new_collision_count_0 == original_collision_count_0);
  assert(new_collision_count_1 == original_collision_count_1);
}

void PenaltyImpactFrictionMap::deleteCacheEntries(const VectorXs &q) {
#ifndef NDEBUG
  const int num_collisions_initial_0 = m_collision_cache.numTotalCollisions();
  const int num_collisions_initial_1 = m_old_cache.numTotalCollisions();
#endif

  assert(q.size() % 3 == 0);
  const int nbodies_initial = q.size() / 3;

  if (m_collision_cache.numBodies() != m_old_cache.numBodies()) {
    std::cerr << "Error, inconsistent body count in "
                 "PenaltyImpactFrictionMap::deleteCacheEntries. This is a bug."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::vector<int> old_to_new_idx_map;
  old_to_new_idx_map.resize(nbodies_initial);

  int copy_to = 0;
  int copy_from = 0;
  for (; copy_from < nbodies_initial; copy_from++, copy_to++) {
    while (copy_from < nbodies_initial && std::isnan(q(3 * copy_from))) {
      old_to_new_idx_map[copy_from] = -1;
      copy_from++;
    }
    if (copy_from == nbodies_initial) {
      break;
    }

    old_to_new_idx_map[copy_from] = copy_to;

    if (!m_old_cache.empty()) {
      m_collision_cache.circle_half_plane_collisions.collisions(copy_to) =
          std::move(m_collision_cache.circle_half_plane_collisions.collisions(
              copy_from));
      m_collision_cache.circle_drum_collisions.collisions(copy_to) = std::move(
          m_collision_cache.circle_drum_collisions.collisions(copy_from));
      m_collision_cache.circle_circle_collisions.collisions(copy_to) =
          std::move(
              m_collision_cache.circle_circle_collisions.collisions(copy_from));
      m_collision_cache.kinematic_circle_circle_collisions.collisions(copy_to) =
          std::move(
              m_collision_cache.kinematic_circle_circle_collisions.collisions(
                  copy_from));
      m_collision_cache.kinematic_box_circle_collisions.collisions(copy_to) =
          std::move(
              m_collision_cache.kinematic_box_circle_collisions.collisions(
                  copy_from));

      m_old_cache.circle_half_plane_collisions.collisions(copy_to) = std::move(
          m_old_cache.circle_half_plane_collisions.collisions(copy_from));
      m_old_cache.circle_drum_collisions.collisions(copy_to) =
          std::move(m_old_cache.circle_drum_collisions.collisions(copy_from));
      m_old_cache.circle_circle_collisions.collisions(copy_to) =
          std::move(m_old_cache.circle_circle_collisions.collisions(copy_from));
      m_old_cache.kinematic_circle_circle_collisions.collisions(copy_to) =
          std::move(m_old_cache.kinematic_circle_circle_collisions.collisions(
              copy_from));
      m_old_cache.kinematic_box_circle_collisions.collisions(copy_to) =
          std::move(m_old_cache.kinematic_box_circle_collisions.collisions(
              copy_from));
    }
  }

  if (!m_old_cache.empty()) {
    m_collision_cache.resizeNumBodies(copy_to);
    m_old_cache.resizeNumBodies(copy_to);

    m_collision_cache.relabelEntries(old_to_new_idx_map);
    m_old_cache.relabelEntries(old_to_new_idx_map);
  }

#ifndef NDEBUG
  const int num_collisions_final_0 = m_collision_cache.numTotalCollisions();
  const int num_collisions_final_1 = m_old_cache.numTotalCollisions();
#endif
  assert(num_collisions_final_0 <= num_collisions_initial_0);
  assert(num_collisions_final_1 <= num_collisions_initial_1);
}

CollisionCache &PenaltyImpactFrictionMap::getCollisionCache() {
  return m_collision_cache;
}

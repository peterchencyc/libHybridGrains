#include "StaticDrumCircleConstraint.h"

#include "scisim/Math/MathUtilities.h"

#include "RigidBody2DStaticDrum.h"

bool StaticDrumCircleConstraint::isActive(const Vector2s &x_circle,
                                          const scalar &r_circle,
                                          const Vector2s &x_drum,
                                          const scalar &r_drum) {
  return (x_drum - x_circle).squaredNorm() >=
         (r_drum - r_circle) * (r_drum - r_circle);
}

StaticDrumCircleConstraint::StaticDrumCircleConstraint(
    const unsigned body_idx, const unsigned drum_idx, const scalar &circle_r,
    const RigidBody2DStaticDrum &drum)
    : m_circle_idx(body_idx), m_r(circle_r), m_drum(drum),
      m_drum_idx(drum_idx) {
  assert(m_r >= 0.0);
}

scalar StaticDrumCircleConstraint::evalNdotV(const VectorXs &q,
                                             const VectorXs &v) const {
  const Vector2s n = (m_drum.x() - q.segment<2>(3 * m_circle_idx)).normalized();
  assert(fabs(n.norm() - 1.0) <= 1.0e-6);
  return n.dot(v.segment<2>(3 * m_circle_idx) -
               computeDrumCollisionPointVelocity(q));
}

int StaticDrumCircleConstraint::impactStencilSize() const { return 2; }

void StaticDrumCircleConstraint::getSimulatedBodyIndices(
    std::pair<int, int> &bodies) const {
  bodies.first = m_circle_idx;
  bodies.second = -1;
}

void StaticDrumCircleConstraint::getBodyIndices(
    std::pair<int, int> &bodies) const {
  bodies.first = m_circle_idx;
  bodies.second = -1;
}

void StaticDrumCircleConstraint::evalH(const VectorXs &q,
                                       const MatrixXXsc &basis, MatrixXXsc &H0,
                                       MatrixXXsc &H1) const {
  assert(H0.rows() == 2);
  assert(H0.cols() == 3);
  assert(H1.rows() == 2);
  assert(H1.cols() == 3);
  assert((basis * basis.transpose() - MatrixXXsc::Identity(2, 2))
             .lpNorm<Eigen::Infinity>() <= 1.0e-6);
  assert(fabs(basis.determinant() - 1.0) <= 1.0e-6);

  // Grab the contact normal
  const Vector2s n{basis.col(0)};
  // Grab the tangent basis
  const Vector2s t{basis.col(1)};

  // Compute the displacement from the center of mass to the point of contact
  assert(m_r >= 0.0);
  const Vector2s r_world{-m_r * n};

  H0.block<1, 2>(0, 0) = n;
  H0(0, 2) = 0.0;

  H0.block<1, 2>(1, 0) = t;
  H0(1, 2) = MathUtilities::cross(r_world, t);
}

bool StaticDrumCircleConstraint::conservesTranslationalMomentum() const {
  return false;
}

bool StaticDrumCircleConstraint::conservesAngularMomentumUnderImpact() const {
  return false;
}

bool StaticDrumCircleConstraint::
    conservesAngularMomentumUnderImpactAndFriction() const {
  return false;
}

std::string StaticDrumCircleConstraint::name() const {
  return "static_drum_circle";
}

Vector2s StaticDrumCircleConstraint::computeDrumCollisionPointVelocity(
    const VectorXs &q) const {
  VectorXs contact_point;
  getWorldSpaceContactPoint(q, contact_point);
  assert(contact_point.size() == 2);
  const Vector2s collision_arm{contact_point - m_drum.x()};
  const Vector2s t0{-collision_arm.y(), collision_arm.x()};
  return m_drum.v() + m_drum.omega() * t0;
}

void StaticDrumCircleConstraint::computeContactBasis(const VectorXs &q,
                                                     const VectorXs &v,
                                                     MatrixXXsc &basis) const {
  const Vector2s n = (m_drum.x() - q.segment<2>(3 * m_circle_idx)).normalized();
  assert(fabs(n.norm() - 1.0) <= 1.0e-6);
  const Vector2s t{-n.y(), n.x()};
  assert(fabs(t.norm() - 1.0) <= 1.0e-6);
  assert(fabs(n.dot(t)) <= 1.0e-6);

  basis.resize(2, 2);
  basis.col(0) = n;
  basis.col(1) = t;
}

VectorXs
StaticDrumCircleConstraint::computeRelativeVelocity(const VectorXs &q,
                                                    const VectorXs &v) const {
  const Vector2s n = (m_drum.x() - q.segment<2>(3 * m_circle_idx)).normalized();

  // Point of contact relative to first body's center of mass
  const Vector2s r0{-m_r * n};
  // Rotate 90 degrees counter clockwise for computing the torque
  const Vector2s t0{-r0.y(), r0.x()};

  // v_point + omega_point x r_point - v_plane_collision_point
  return v.segment<2>(3 * m_circle_idx) + v(3 * m_circle_idx + 2) * t0 -
         computeDrumCollisionPointVelocity(q);
}

void StaticDrumCircleConstraint::setBodyIndex0(const unsigned idx) {
  m_circle_idx = idx;
}

scalar
StaticDrumCircleConstraint::computePenetrationDepth(const VectorXs &q) const {
  assert(3 * m_circle_idx + 2 < q.size());
  return std::min(0.0, -((m_drum.x() - q.segment<2>(3 * m_circle_idx)).norm() -
                         (m_drum.r() - m_r)));
}

VectorXs StaticDrumCircleConstraint::computeKinematicRelativeVelocity(
    const VectorXs &q, const VectorXs &v) const {
  return computeDrumCollisionPointVelocity(q);
}

void StaticDrumCircleConstraint::getWorldSpaceContactPoint(
    const VectorXs &q, VectorXs &contact_point) const {
  const Vector2s n = (m_drum.x() - q.segment<2>(3 * m_circle_idx)).normalized();
  contact_point = q.segment<2>(3 * m_circle_idx) - m_r * n;
}

void StaticDrumCircleConstraint::getWorldSpaceContactNormal(
    const VectorXs &q, VectorXs &contact_normal) const {
  contact_normal = (m_drum.x() - q.segment<2>(3 * m_circle_idx)).normalized();
}

unsigned StaticDrumCircleConstraint::getStaticObjectIndex() const {
  return m_drum_idx;
}

const scalar &StaticDrumCircleConstraint::circleRadius() const { return m_r; }

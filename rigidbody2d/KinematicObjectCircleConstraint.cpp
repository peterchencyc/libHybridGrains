#include "KinematicObjectCircleConstraint.h"

#include "scisim/Math/MathUtilities.h"

#include "CircleBoxTools.h"

KinematicObjectCircleConstraint::KinematicObjectCircleConstraint(
    const unsigned sim_bdy_idx, const scalar &sim_bdy_r, const Vector2s &n,
    const unsigned knmtc_bdy_idx, const Vector2s &x, const Vector2s &v,
    const scalar &omega)
    : m_sim_idx(sim_bdy_idx), m_r(sim_bdy_r), m_n(n),
      m_kinematic_index(knmtc_bdy_idx), m_kinematic_x(x), m_v(v),
      m_omega(omega) {
  assert(m_r > 0.0);
  assert(fabs(m_n.norm() - 1.0) <= 1.0e-6);
}

scalar KinematicObjectCircleConstraint::evalNdotV(const VectorXs &q,
                                                  const VectorXs &v) const {
  assert(v.size() % 3 == 0);
  assert(3 * m_sim_idx + 1 < v.size());
  return m_n.dot(v.segment<2>(3 * m_sim_idx) -
                 computeKinematicCollisionPointVelocity(q));
}

void KinematicObjectCircleConstraint::evalgradg(
    const VectorXs &q, const int col, SparseMatrixsc &G,
    const FlowableSystem &fsys) const {
  assert(col >= 0);
  assert(col < G.cols());

  // MUST BE ADDED GOING DOWN THE COLUMN. DO NOT TOUCH ANOTHER COLUMN.
  assert(3 * m_sim_idx + 1 < unsigned(G.rows()));
  G.insert(3 * m_sim_idx + 0, col) = m_n.x();
  G.insert(3 * m_sim_idx + 1, col) = m_n.y();
}

int KinematicObjectCircleConstraint::impactStencilSize() const { return 2; }

void KinematicObjectCircleConstraint::getSimulatedBodyIndices(
    std::pair<int, int> &bodies) const {
  bodies.first = m_sim_idx;
  bodies.second = -1;
}

void KinematicObjectCircleConstraint::getBodyIndices(
    std::pair<int, int> &bodies) const {
  bodies.first = m_sim_idx;
  bodies.second = m_kinematic_index;
}

void KinematicObjectCircleConstraint::evalKinematicNormalRelVel(
    const VectorXs &q, const int strt_idx, VectorXs &gdotN) const {
  assert(strt_idx >= 0);
  assert(strt_idx < gdotN.size());
  gdotN(strt_idx) = -m_n.dot(computeKinematicCollisionPointVelocity(q));
}

void KinematicObjectCircleConstraint::evalH(const VectorXs &q,
                                            const MatrixXXsc &basis,
                                            MatrixXXsc &H0,
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

  // Generate arms
  const Vector2s r{-m_r * n};

  // Format for H:
  //   n^T  r x n
  //   t^T  r x t

  H0.block<1, 2>(0, 0) = n;
  assert(fabs(MathUtilities::cross(r, n)) <= 1.0e-6);
  H0(0, 2) = 0.0;

  H0.block<1, 2>(1, 0) = t;
  H0(1, 2) = MathUtilities::cross(r, t);
}

void KinematicObjectCircleConstraint::computeContactBasis(
    const VectorXs &q, const VectorXs &v, MatrixXXsc &basis) const {
  assert(fabs(m_n.norm() - 1.0) <= 1.0e-6);
  const Vector2s t{-m_n.y(), m_n.x()};
  assert(fabs(t.norm() - 1.0) <= 1.0e-6);
  assert(fabs(m_n.dot(t)) <= 1.0e-6);
  basis.resize(2, 2);
  basis.col(0) = m_n;
  basis.col(1) = t;
}

bool KinematicObjectCircleConstraint::conservesTranslationalMomentum() const {
  return false;
}

bool KinematicObjectCircleConstraint::conservesAngularMomentumUnderImpact()
    const {
  return false;
}

bool KinematicObjectCircleConstraint::
    conservesAngularMomentumUnderImpactAndFriction() const {
  return false;
}

std::string KinematicObjectCircleConstraint::name() const {
  return "kinematic_object_circle";
}

Vector2s
KinematicObjectCircleConstraint::computeKinematicCollisionPointVelocity(
    const VectorXs &q) const {
  VectorXs contact_point;
  getWorldSpaceContactPoint(q, contact_point);
  assert(contact_point.size() == 2);
  const Vector2s collision_arm{contact_point - m_kinematic_x};
  const Vector2s t0{-collision_arm.y(), collision_arm.x()};
  return m_v + m_omega * t0;
}

VectorXs KinematicObjectCircleConstraint::computeRelativeVelocity(
    const VectorXs &q, const VectorXs &v) const {
  assert(v.size() % 3 == 0);
  assert(3 * m_sim_idx + 2 < v.size());

  // Point of contact relative to the simulated body's center of mass
  const Vector2s r{-m_r * m_n};
  assert(fabs(MathUtilities::cross(m_n, r)) <= 1.0e-6);

  // Rotate 90 degrees counter clockwise for computing the torque
  const Vector2s t{-r.y(), r.x()};

  // v + omega x r - kinematic_vel
  return v.segment<2>(3 * m_sim_idx) + v(3 * m_sim_idx + 2) * t -
         computeKinematicCollisionPointVelocity(q);
}

void KinematicObjectCircleConstraint::setBodyIndex0(const unsigned idx) {
  m_sim_idx = idx;
}

VectorXs KinematicObjectCircleConstraint::computeKinematicRelativeVelocity(
    const VectorXs &q, const VectorXs &v) const {
  return computeKinematicCollisionPointVelocity(q);
}

void KinematicObjectCircleConstraint::getWorldSpaceContactPoint(
    const VectorXs &q, VectorXs &contact_point) const {
  contact_point = q.segment<2>(3 * m_sim_idx) - m_r * m_n;
}

void KinematicObjectCircleConstraint::getWorldSpaceContactNormal(
    const VectorXs &q, VectorXs &contact_normal) const {
  contact_normal = m_n;
}

KinematicCircleCircleConstraint::KinematicCircleCircleConstraint(
    const unsigned sim_bdy_idx, const scalar &sim_bdy_r, const Vector2s &n,
    const unsigned knmtc_bdy_idx, const Vector2s &x, const Vector2s &v,
    const scalar &omega, const scalar &r_kinematic)
    : KinematicObjectCircleConstraint(sim_bdy_idx, sim_bdy_r, n, knmtc_bdy_idx,
                                      x, v, omega),
      m_r_kinematic(r_kinematic) {
  assert(m_r_kinematic >= 0.0);
}

scalar
KinematicCircleCircleConstraint::evaluateGapFunction(const VectorXs &q) const {
  return (q.segment<2>(3 * m_sim_idx) - q.segment<2>(3 * m_kinematic_index))
             .norm() -
         m_r - m_r_kinematic;
}

std::string KinematicCircleCircleConstraint::name() const {
  return "kinematic_circle_circle";
}

scalar KinematicCircleCircleConstraint::computePenetrationDepth(
    const VectorXs &q) const {
  return std::min(
      0.0, (q.segment<2>(3 * m_sim_idx) - q.segment<2>(3 * m_kinematic_index))
                   .norm() -
               m_r - m_r_kinematic);
}

std::pair<VectorXs, VectorXs>
KinematicCircleCircleConstraint::computeLeverArms(const VectorXs &q) const {
  return std::make_pair<VectorXs, VectorXs>(-m_r * m_n, VectorXs());
}

KinematicBoxCircleConstraint::KinematicBoxCircleConstraint(
    const unsigned sim_bdy_idx, const scalar &sim_bdy_r, const Vector2s &n,
    const unsigned knmtc_bdy_idx, const Vector2s &x, const Vector2s &v,
    const scalar &omega, const Vector2s &r_box)
    : KinematicObjectCircleConstraint(sim_bdy_idx, sim_bdy_r, n, knmtc_bdy_idx,
                                      x, v, omega),
      m_r_box(r_box) {
  assert((r_box.array() >= 0.0).all());
}

std::string KinematicBoxCircleConstraint::name() const {
  return "kinematic_box_circle";
}

scalar
KinematicBoxCircleConstraint::computePenetrationDepth(const VectorXs &q) const {
  Vector2s n;
  Vector2s p;
  scalar pen_depth;
  CircleBoxTools::isActive(
      q.segment<2>(3 * m_sim_idx), m_r, q.segment<2>(3 * m_kinematic_index),
      q(3 * m_kinematic_index + 2), m_r_box, n, p, pen_depth);
  return pen_depth;
}

std::pair<VectorXs, VectorXs>
KinematicBoxCircleConstraint::computeLeverArms(const VectorXs &q) const {
  return std::make_pair<VectorXs, VectorXs>(-m_r * m_n, VectorXs());
}

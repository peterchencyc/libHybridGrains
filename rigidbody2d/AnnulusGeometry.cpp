#include "AnnulusGeometry.h"

#include "scisim/Utilities.h"

AnnulusGeometry::AnnulusGeometry(const scalar &r0, const scalar &r1)
    : m_r0(r0), m_r1(r1) {
  assert(m_r0 >= 0.0);
  assert(m_r1 > 0.0);
  assert(m_r0 <= m_r1);
}

AnnulusGeometry::AnnulusGeometry(std::istream &input_stream)
    : m_r0(Utilities::deserialize<scalar>(input_stream)),
      m_r1(Utilities::deserialize<scalar>(input_stream)) {
  assert(m_r0 >= 0.0);
  assert(m_r1 > 0.0);
  assert(m_r0 <= m_r1);
}

RigidBody2DGeometryType AnnulusGeometry::type() const {
  return RigidBody2DGeometryType::ANNULUS;
}

std::unique_ptr<RigidBody2DGeometry> AnnulusGeometry::clone() const {
  return std::unique_ptr<RigidBody2DGeometry>{new AnnulusGeometry{m_r0, m_r1}};
}

void AnnulusGeometry::computeCollisionAABB(const Vector2s &x0,
                                           const scalar &theta0,
                                           const Vector2s &x1,
                                           const scalar &theta1, Array2s &min,
                                           Array2s &max) const {
  min = x0.array().min(x1.array()) - m_r1;
  max = x0.array().max(x1.array()) + m_r1;
  assert((min <= max).all());
}

scalar AnnulusGeometry::computeArea() const {
  return MathDefines::PI<scalar>() * (m_r1 * m_r1 - m_r0 * m_r0);
}

void AnnulusGeometry::AABB(const Vector2s &x, const scalar &theta, Array2s &min,
                           Array2s &max) const {
  min = x.array() - m_r1;
  max = x.array() + m_r1;
}

void AnnulusGeometry::massAndInertia(const scalar &density, scalar &m,
                                     scalar &I) const {
  m = density * MathDefines::PI<scalar>() * (m_r1 * m_r1 - m_r0 * m_r0);
  I = 0.5 * m * (m_r1 * m_r1 + m_r0 * m_r0);
}

void AnnulusGeometry::serializeState(std::ostream &output_stream) const {
  Utilities::serializeBuiltInType(RigidBody2DGeometryType::ANNULUS,
                                  output_stream);
  Utilities::serializeBuiltInType(m_r0, output_stream);
  Utilities::serializeBuiltInType(m_r1, output_stream);
}

const scalar &AnnulusGeometry::r0() const { return m_r0; }

const scalar &AnnulusGeometry::r1() const { return m_r1; }

#include "CircleGeometry.h"

#include "scisim/Utilities.h"

CircleGeometry::CircleGeometry(const scalar &r) : m_r(r) { assert(m_r > 0.0); }

CircleGeometry::CircleGeometry(std::istream &input_stream)
    : m_r(Utilities::deserialize<scalar>(input_stream)) {
  assert(m_r > 0.0);
}

RigidBody2DGeometryType CircleGeometry::type() const {
  return RigidBody2DGeometryType::CIRCLE;
}

std::unique_ptr<RigidBody2DGeometry> CircleGeometry::clone() const {
  return std::unique_ptr<RigidBody2DGeometry>{new CircleGeometry{m_r}};
}

const scalar &CircleGeometry::r() const { return m_r; }

void CircleGeometry::computeCollisionAABB(const Vector2s &x0,
                                          const scalar &theta0,
                                          const Vector2s &x1,
                                          const scalar &theta1, Array2s &min,
                                          Array2s &max) const {
  min = x0.array().min(x1.array()) - m_r;
  max = x0.array().max(x1.array()) + m_r;
  assert((min <= max).all());
}

scalar CircleGeometry::computeArea() const {
  return MathDefines::PI<scalar>() * m_r * m_r;
}

void CircleGeometry::AABB(const Vector2s &x, const scalar &theta, Array2s &min,
                          Array2s &max) const {
  // Circles are invariant to rotations about the center of mass, so theta is
  // unused here
  min = x.array() - m_r;
  max = x.array() + m_r;
}

void CircleGeometry::massAndInertia(const scalar &density, scalar &m,
                                    scalar &I) const {
  m = density * MathDefines::PI<scalar>() * m_r * m_r;
  I = 0.5 * m * m_r * m_r;
}

void CircleGeometry::serializeState(std::ostream &output_stream) const {
  Utilities::serializeBuiltInType(RigidBody2DGeometryType::CIRCLE,
                                  output_stream);
  Utilities::serializeBuiltInType(m_r, output_stream);
}

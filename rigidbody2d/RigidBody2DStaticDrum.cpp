#include "RigidBody2DStaticDrum.h"

RigidBody2DStaticDrum::RigidBody2DStaticDrum(const Vector2s &x, const scalar &r)
    : m_x(x), m_r(r), m_v(Vector2s::Zero()), m_omega(0.0), m_theta(0.0) {
  assert(m_r > 0.0);
}

const Vector2s &RigidBody2DStaticDrum::x() const { return m_x; }

const scalar &RigidBody2DStaticDrum::r() const { return m_r; }

const Vector2s &RigidBody2DStaticDrum::v() const { return m_v; }

const scalar &RigidBody2DStaticDrum::omega() const { return m_omega; }

const scalar &RigidBody2DStaticDrum::theta() const { return m_theta; }

void RigidBody2DStaticDrum::setX(const Vector2s &x) { m_x = x; }

void RigidBody2DStaticDrum::setV(const Vector2s &v) { m_v = v; }

void RigidBody2DStaticDrum::setTheta(const scalar &theta) { m_theta = theta; }

void RigidBody2DStaticDrum::setOmega(const scalar &omega) { m_omega = omega; }

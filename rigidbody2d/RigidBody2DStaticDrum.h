#ifndef RIGID_BODY_2D_STATIC_DRUM_H
#define RIGID_BODY_2D_STATIC_DRUM_H

#include "scisim/Math/MathDefines.h"

class RigidBody2DStaticDrum final {

public:
  RigidBody2DStaticDrum(const Vector2s &x, const scalar &r);

  const Vector2s &x() const;

  const scalar &r() const;

  const Vector2s &v() const;

  const scalar &omega() const;

  const scalar &theta() const;

  void setX(const Vector2s &x);

  void setV(const Vector2s &v);

  void setTheta(const scalar &theta);

  void setOmega(const scalar &omega);

private:
  Vector2s m_x;
  scalar m_r;
  Vector2s m_v;
  scalar m_omega;
  scalar m_theta;
};

#endif

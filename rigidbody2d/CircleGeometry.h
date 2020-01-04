#ifndef CIRCLE_GEOMETRY_H
#define CIRCLE_GEOMETRY_H

#include "RigidBody2DGeometry.h"

class CircleGeometry final : public RigidBody2DGeometry {

public:
  explicit CircleGeometry(const scalar &r);
  explicit CircleGeometry(std::istream &input_stream);

  virtual ~CircleGeometry() override = default;

  virtual RigidBody2DGeometryType type() const override;

  virtual std::unique_ptr<RigidBody2DGeometry> clone() const override;

  virtual void computeCollisionAABB(const Vector2s &x0, const scalar &theta0,
                                    const Vector2s &x1, const scalar &theta1,
                                    Array2s &min, Array2s &max) const override;

  virtual scalar computeArea() const override;

  const scalar &r() const;

private:
  virtual void AABB(const Vector2s &x, const scalar &theta, Array2s &min,
                    Array2s &max) const override;

  virtual void massAndInertia(const scalar &density, scalar &m,
                              scalar &I) const override;

  virtual void serializeState(std::ostream &output_stream) const override;

  const scalar m_r;
};

#endif

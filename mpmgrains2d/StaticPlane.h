#ifndef MPM_STATIC_PLANE_2D_H
#define MPM_STATIC_PLANE_2D_H

#include "scisim/Math/MathDefines.h"

struct MPMStaticPlane final {
  enum class BoundaryBehavior { SLIDING, STICKING };

  MPMStaticPlane(const Vector2s &point, const Vector2s &normal,
                 const BoundaryBehavior boundary_treatment,
                 const scalar &lower_bound);
  MPMStaticPlane(std::istream &input_stream);

  void serialize(std::ostream &output_stream) const;

  double distanceToPoint(const Vector2s &p) const;

  Vector2s x;
  Vector2s n;
  BoundaryBehavior boundary_behavior;
  scalar lower_bound;
};

#endif

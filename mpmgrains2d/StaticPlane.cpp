#include "StaticPlane.h"

#include "scisim/Math/MathUtilities.h"
#include "scisim/Utilities.h"

MPMStaticPlane::MPMStaticPlane(const Vector2s &point, const Vector2s &normal,
                               const BoundaryBehavior boundary_treatment,
                               const scalar &lower_bound_in)
    : x(point), n(normal.normalized()), boundary_behavior(boundary_treatment),
      lower_bound(lower_bound_in) {
  assert(fabs(n.norm() - 1.0) <= 1.0e-9);
}

MPMStaticPlane::MPMStaticPlane(std::istream &input_stream)
    : x(MathUtilities::deserialize<Vector2s>(input_stream)),
      n(MathUtilities::deserialize<Vector2s>(input_stream)),
      boundary_behavior(Utilities::deserialize<BoundaryBehavior>(input_stream)),
      lower_bound(Utilities::deserialize<scalar>(input_stream)) {
  assert(fabs(n.norm() - 1.0) <= 1.0e-9);
}

void MPMStaticPlane::serialize(std::ostream &output_stream) const {
  MathUtilities::serialize(x, output_stream);
  MathUtilities::serialize(n, output_stream);
  Utilities::serializeBuiltInType(boundary_behavior, output_stream);
  Utilities::serializeBuiltInType(lower_bound, output_stream);
}

double MPMStaticPlane::distanceToPoint(const Vector2s &p) const {
  return n.dot(p - x);
}

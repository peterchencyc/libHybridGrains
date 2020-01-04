#include "uGIMPLinearBasisFunctions.h"

#include "MaterialPoints.h"
#include "PhysicsGrid.h"

static scalar LinearIntegral(const scalar &xp, const scalar &hl,
                             const scalar &xi, const scalar &w) {
  assert(w > 2.0 * hl);
  assert(w > 0.0);
  assert(hl > 0.0);

  const scalar diff = fabs(xp - xi);
  if (diff >= w + hl) {
    return 0.0;
  } else if (diff >= w - hl) {
    return (w + hl - diff) * (w + hl - diff) / (2.0 * w);
  } else if (diff >= hl) {
    return 2.0 * hl * (1.0 - diff / w);
  } else {
    return 2.0 * hl - (hl * hl + diff * diff) / w;
  }
}

static scalar LinearIntegralGrad(const scalar &xp, const scalar &hl,
                                 const scalar &xi, const scalar &w) {
  assert(w > 2.0 * hl);
  assert(w > 0.0);
  assert(hl > 0.0);

  const scalar diff = fabs(xp - xi);
  const scalar sgn = (xp - xi) >= 0.0 ? 1.0 : -1.0;

  if (diff >= w + hl) {
    return 0.0;
  } else if (diff >= w - hl) {
    return -sgn * (w + hl - diff) / w;
  } else if (diff >= hl) {
    return -sgn * 2.0 * hl / w;
  } else {
    return 2.0 * (xi - xp) / w;
  }
}

std::unique_ptr<BasisFunctions> uGIMPLinearBasisFunctions::clone() const {
  return std::make_unique<uGIMPLinearBasisFunctions>();
}

BasisFunctionType uGIMPLinearBasisFunctions::type() const {
  return BasisFunctionType::uGIMPLinear;
}

std::pair<Array2u, Array2u>
uGIMPLinearBasisFunctions::computeStencil(const unsigned int pidx,
                                          const MaterialPoints &points,
                                          const PhysicsGrid &grid) const {
  const Vector2s xp = points.q.col(pidx);
  const scalar hl = points.hl(pidx);
  return computeStencil(xp, hl, grid.min, grid.cell_width);
}

std::pair<Array2u, Array2u>
uGIMPLinearBasisFunctions::computeStencil(const Vector2s &xp,
                                          const scalar &half_width,
                                          const PhysicsGrid &grid) const {
  return computeStencil(xp, half_width, grid.min, grid.cell_width);
}

std::pair<Array2u, Array2u> uGIMPLinearBasisFunctions::computeStencil(
    const Vector2s &xp, const scalar &half_width, const Vector2s &grid_min,
    const scalar &grid_cell_width) const {
  std::pair<Array2u, Array2u> stencil;

  const Vector2s xp_m = xp.array() - half_width;
  const Vector2s xp_M = xp.array() + half_width;

  stencil.first = ((xp_m.array() - grid_min.array()) / grid_cell_width)
                      .unaryExpr([](const scalar &s) {
                        using std::floor;
                        return floor(s);
                      })
                      .cast<unsigned>();
  stencil.second = ((xp_M.array() - grid_min.array()) / grid_cell_width)
                       .unaryExpr([](const scalar &s) {
                         using std::ceil;
                         return ceil(s);
                       })
                       .cast<unsigned>();
  ;
  return stencil;
}

scalar uGIMPLinearBasisFunctions::weight(const unsigned int pidx,
                                         const MaterialPoints &points,
                                         const Vector2u &nidx,
                                         const PhysicsGrid &grid) const {
  const Vector2s xp = points.q.col(pidx);
  const scalar hl = points.hl(pidx);
  return weight(xp, hl, nidx, grid.min, grid.cell_width);
}

scalar uGIMPLinearBasisFunctions::weight(const Vector2s &xp,
                                         const scalar &half_width,
                                         const Vector2u &nidx,
                                         const PhysicsGrid &grid) const {
  return weight(xp, half_width, nidx, grid.min, grid.cell_width);
}

scalar uGIMPLinearBasisFunctions::weight(const Vector2s &xp,
                                         const scalar &half_width,
                                         const Vector2u &nidx,
                                         const Vector2s &grid_min,
                                         const scalar &grid_cell_width) const {
  const Vector2s xi = grid_cell_width * nidx.cast<scalar>() + grid_min;
  const scalar w = grid_cell_width;
  const scalar Vp = 4.0 * half_width * half_width;

  const scalar weight{LinearIntegral(xp.x(), half_width, xi.x(), w) *
                      LinearIntegral(xp.y(), half_width, xi.y(), w) / Vp};

  assert(weight >= 0.0);
  assert(weight <= 1.0);
  return weight;
}

Vector2s uGIMPLinearBasisFunctions::weightGrad(const unsigned int pidx,
                                               const MaterialPoints &points,
                                               const Vector2u &nidx,
                                               const PhysicsGrid &grid) const {
  const Vector2s xp = points.q.col(pidx);
  const scalar hl = points.hl(pidx);
  return weightGrad(xp, hl, nidx, grid.min, grid.cell_width);
}

Vector2s uGIMPLinearBasisFunctions::weightGrad(const Vector2s &xp,
                                               const scalar &half_width,
                                               const Vector2u &nidx,
                                               const PhysicsGrid &grid) const {
  return weightGrad(xp, half_width, nidx, grid.min, grid.cell_width);
}

Vector2s uGIMPLinearBasisFunctions::weightGrad(
    const Vector2s &xp, const scalar &half_width, const Vector2u &nidx,
    const Vector2s &grid_min, const scalar &grid_cell_width) const {
  const Vector2s xi = grid_cell_width * nidx.cast<scalar>() + grid_min.matrix();
  const scalar w = grid_cell_width;
  const scalar Vp = 4.0 * half_width * half_width;

  const scalar wx{LinearIntegral(xp.x(), half_width, xi.x(), w)};
  const scalar wy{LinearIntegral(xp.y(), half_width, xi.y(), w)};

  Vector2s grad;
  grad.x() = wy * LinearIntegralGrad(xp.x(), half_width, xi.x(), w) / Vp;
  grad.y() = wx * LinearIntegralGrad(xp.y(), half_width, xi.y(), w) / Vp;

  return grad;
}

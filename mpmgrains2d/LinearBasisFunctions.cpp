#include "LinearBasisFunctions.h"

#include "MaterialPoints.h"
#include "PhysicsGrid.h"

static scalar linearN(const scalar &x) {
  const scalar absx{fabs(x)};

  if (absx < 1.0) {
    return 1.0 - absx;
  } else {
    return 0.0;
  }
}

// TODO: Also pass this 1.0 / grid.cell_width to avoid the divide in weightGrad
static scalar lineardN(const scalar &x) {
  if (fabs(x) <= 1.0) {
    return (x >= 0.0) ? -1.0 : 1.0;
  } else {
    return 0.0;
  }
}

std::unique_ptr<BasisFunctions> LinearBasisFunctions::clone() const {
  return std::make_unique<LinearBasisFunctions>();
}

BasisFunctionType LinearBasisFunctions::type() const {
  return BasisFunctionType::Linear;
}

std::pair<Array2u, Array2u>
LinearBasisFunctions::computeStencil(const unsigned int pidx,
                                     const MaterialPoints &points,
                                     const PhysicsGrid &grid) const {
  constexpr scalar half_width{0.0}; // NB: This is unused
  return computeStencil(points.q.col(pidx), half_width, grid.min,
                        grid.cell_width);
}

std::pair<Array2u, Array2u>
LinearBasisFunctions::computeStencil(const Vector2s &xp,
                                     const scalar &half_width,
                                     const PhysicsGrid &grid) const {
  return computeStencil(xp, half_width, grid.min, grid.cell_width);
}

std::pair<Array2u, Array2u> LinearBasisFunctions::computeStencil(
    const Vector2s &xp, const scalar &half_width, const Vector2s &grid_min,
    const scalar &grid_cell_width) const {
  std::pair<Array2u, Array2u> stencil;
  using std::floor;
  stencil.first = ((xp.array() - grid_min.array()) / grid_cell_width)
                      .unaryExpr([](const scalar &s) {
                        using std::floor;
                        return floor(s);
                      })
                      .cast<unsigned>();
  stencil.second = stencil.first + 1;
  return stencil;
}

scalar LinearBasisFunctions::weight(const unsigned int pidx,
                                    const MaterialPoints &points,
                                    const Vector2u &nidx,
                                    const PhysicsGrid &grid) const {
  constexpr scalar half_width{0.0}; // NB: This is unused
  return weight(points.q.col(pidx), half_width, nidx, grid.min,
                grid.cell_width);
}

scalar LinearBasisFunctions::weight(const Vector2s &xp,
                                    const scalar &half_width,
                                    const Vector2u &nidx,
                                    const PhysicsGrid &grid) const {
  return weight(xp, half_width, nidx, grid.min, grid.cell_width);
}

scalar LinearBasisFunctions::weight(const Vector2s &xp,
                                    const scalar &half_width,
                                    const Vector2u &nidx,
                                    const Vector2s &grid_min,
                                    const scalar &grid_cell_width) const {
  // TODO: Cache 1 / cell_width^2? Can replace two divides with one?
  const Vector2s delta_over_h{
      (xp - grid_cell_width * nidx.cast<scalar>() - grid_min) /
      grid_cell_width};
  // TODO: If any of these things is less than 1 or greater than 1, return 0
  const scalar weight{linearN(delta_over_h.x()) * linearN(delta_over_h.y())};
  // return (1.0 - fabs(x)) * (1.0 - fabs(y)) * (1.0 - fabs(z))
  // OR OR OR, if violated, final value is negative. Just clamp negative values
  // to 0 and return?
  assert(weight >= 0.0);
  assert(weight <= 1.0);
  return weight;
}

Vector2s LinearBasisFunctions::weightGrad(const unsigned int pidx,
                                          const MaterialPoints &points,
                                          const Vector2u &nidx,
                                          const PhysicsGrid &grid) const {
  constexpr scalar half_width{0.0}; // NB: This is unused
  return weightGrad(points.q.col(pidx), half_width, nidx, grid.min,
                    grid.cell_width);
}

Vector2s LinearBasisFunctions::weightGrad(const Vector2s &xp,
                                          const scalar &half_width,
                                          const Vector2u &nidx,
                                          const PhysicsGrid &grid) const {
  return weightGrad(xp, half_width, nidx, grid.min, grid.cell_width);
}

Vector2s LinearBasisFunctions::weightGrad(const Vector2s &xp,
                                          const scalar &half_width,
                                          const Vector2u &nidx,
                                          const Vector2s &grid_min,
                                          const scalar &grid_cell_width) const {
  // TODO: cache 1 / cell_width in PhysicsGrid
  const Vector2s delta_over_h{
      (xp - grid_cell_width * nidx.cast<scalar>() - grid_min) /
      grid_cell_width};

  const scalar wx{linearN(delta_over_h.x())};
  const scalar wy{linearN(delta_over_h.y())};

  Vector2s grad;
  grad.x() = wy * lineardN(delta_over_h.x());
  grad.y() = wx * lineardN(delta_over_h.y());
  // TODO: What breaks if I nuke this?
  // this division comes from the derivative
  grad /= grid_cell_width;

  return grad;
}

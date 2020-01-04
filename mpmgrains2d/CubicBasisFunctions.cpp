#include "CubicBasisFunctions.h"

#include "MaterialPoints.h"
#include "PhysicsGrid.h"

static scalar cubicN(const scalar &x) {
  const scalar absx{fabs(x)};

  if (absx < 1.0) {
    assert(fabs((0.5 * absx * absx * absx - absx * absx + 2.0 / 3.0) -
                ((0.5 * absx - 1.0) * absx * absx + 2.0 / 3.0)) <= 1.0e-6);
    return (0.5 * absx - 1.0) * absx * absx + 2.0 / 3.0;
  } else if (absx < 2.0) {
    assert(fabs((((-absx / 6.0 + 1.0) * absx - 2.0) * absx + 4.0 / 3.0) -
                (-absx * absx * absx / 6.0 + absx * absx - 2.0 * absx +
                 4.0 / 3.0)) <= 1.0e-6);
    return ((-absx / 6.0 + 1.0) * absx - 2.0) * absx + 4.0 / 3.0;
  } else {
    return 0.0;
  }
}

static scalar cubicdN(const scalar &x) {
  const scalar absx{fabs(x)};
  const scalar sgnx{(x >= 0.0) ? 1.0 : -1.0};

  if (absx < 1.0) {
    // return sgnx * ( 1.5 * absx * absx - 2.0 * absx );
    assert(fabs((sgnx * (1.5 * absx * absx - 2.0 * absx)) -
                (sgnx * (1.5 * absx - 2.0) * absx)) <= 1.0e-6);
    return sgnx * (1.5 * absx - 2.0) * absx;
  } else if (absx < 2.0) {
    // return sgnx * ( - 0.5 * absx * absx + 2.0 * absx - 2.0 );
    assert(fabs((sgnx * (-0.5 * absx * absx + 2.0 * absx - 2.0)) -
                (sgnx * ((-0.5 * absx + 2.0) * absx - 2.0))) <= 1.0e-6);
    return sgnx * ((-0.5 * absx + 2.0) * absx - 2.0);
  } else {
    return 0.0;
  }
}

std::unique_ptr<BasisFunctions> ThirdOrderBasisFunctions::clone() const {
  return std::make_unique<ThirdOrderBasisFunctions>();
}

BasisFunctionType ThirdOrderBasisFunctions::type() const {
  return BasisFunctionType::ThirdOrder;
}

std::pair<Array2u, Array2u>
ThirdOrderBasisFunctions::computeStencil(const unsigned int pidx,
                                         const MaterialPoints &points,
                                         const PhysicsGrid &grid) const {
  constexpr scalar half_width{0.0}; // NB: This is unused
  return computeStencil(points.q.col(pidx), half_width, grid.min,
                        grid.cell_width);
}

std::pair<Array2u, Array2u>
ThirdOrderBasisFunctions::computeStencil(const Vector2s &xp,
                                         const scalar &half_width,
                                         const PhysicsGrid &grid) const {
  return computeStencil(xp, half_width, grid.min, grid.cell_width);
}

std::pair<Array2u, Array2u> ThirdOrderBasisFunctions::computeStencil(
    const Vector2s &xp, const scalar &half_width, const Vector2s &grid_min,
    const scalar &grid_cell_width) const {
  std::pair<Array2u, Array2u> stencil;
  // TODO: Assert that this will not underflow
  stencil.first = ((xp.array() - grid_min.array()) / grid_cell_width)
                      .unaryExpr([](const scalar &s) {
                        using std::floor;
                        return floor(s);
                      })
                      .cast<unsigned>() -
                  1;
  stencil.second = stencil.first + 3;
  assert((stencil.second == stencil.first + 3).all());
  return stencil;
}

scalar ThirdOrderBasisFunctions::weight(const unsigned int pidx,
                                        const MaterialPoints &points,
                                        const Vector2u &nidx,
                                        const PhysicsGrid &grid) const {
  constexpr scalar half_width{0.0}; // NB: This is unused
  return weight(points.q.col(pidx), half_width, nidx, grid.min,
                grid.cell_width);
}

scalar ThirdOrderBasisFunctions::weight(const Vector2s &xp,
                                        const scalar &half_width,
                                        const Vector2u &nidx,
                                        const PhysicsGrid &grid) const {
  return weight(xp, half_width, nidx, grid.min, grid.cell_width);
}

scalar ThirdOrderBasisFunctions::weight(const Vector2s &xp,
                                        const scalar &half_width,
                                        const Vector2u &nidx,
                                        const Vector2s &grid_min,
                                        const scalar &grid_cell_width) const {
  const Vector2s delta_over_h{
      (xp - grid_cell_width * nidx.cast<scalar>() - grid_min) /
      grid_cell_width};
  const scalar weight{cubicN(delta_over_h.x()) * cubicN(delta_over_h.y())};
  assert(weight >= -3.0e-16);
  assert(weight <= 1.0);
  return weight;
}

Vector2s ThirdOrderBasisFunctions::weightGrad(const unsigned int pidx,
                                              const MaterialPoints &points,
                                              const Vector2u &nidx,
                                              const PhysicsGrid &grid) const {
  constexpr scalar half_width{0.0}; // NB: This is unused
  return weightGrad(points.q.col(pidx), half_width, nidx, grid.min,
                    grid.cell_width);
}

Vector2s ThirdOrderBasisFunctions::weightGrad(const Vector2s &xp,
                                              const scalar &half_width,
                                              const Vector2u &nidx,
                                              const PhysicsGrid &grid) const {
  return weightGrad(xp, half_width, nidx, grid.min, grid.cell_width);
}

Vector2s ThirdOrderBasisFunctions::weightGrad(
    const Vector2s &xp, const scalar &half_width, const Vector2u &nidx,
    const Vector2s &grid_min, const scalar &grid_cell_width) const {
  const Vector2s delta_over_h{
      (xp - grid_cell_width * nidx.cast<scalar>() - grid_min) /
      grid_cell_width};

  const scalar wx{cubicN(delta_over_h.x())};
  const scalar wy{cubicN(delta_over_h.y())};

  Vector2s grad;
  grad.x() = wy * cubicdN(delta_over_h.x());
  grad.y() = wx * cubicdN(delta_over_h.y());
  // this division comes from the derivative
  grad /= grid_cell_width;

  return grad;
}

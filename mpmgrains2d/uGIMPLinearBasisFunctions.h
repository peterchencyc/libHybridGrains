#ifndef UGIMP_BASIS_FUNCTIONS_H
#define UGIMP_BASIS_FUNCTIONS_H

#include "BasisFunctions.h"

struct uGIMPLinearBasisFunctions final : public BasisFunctions {

  uGIMPLinearBasisFunctions() = default;

  virtual ~uGIMPLinearBasisFunctions() override = default;

  virtual std::unique_ptr<BasisFunctions> clone() const override;

  virtual BasisFunctionType type() const override;

  virtual std::pair<Array2u, Array2u>
  computeStencil(const unsigned int pidx, const MaterialPoints &points,
                 const PhysicsGrid &grid) const override;
  virtual std::pair<Array2u, Array2u>
  computeStencil(const Vector2s &xp, const scalar &half_width,
                 const PhysicsGrid &grid) const override;
  virtual std::pair<Array2u, Array2u>
  computeStencil(const Vector2s &xp, const scalar &half_width,
                 const Vector2s &grid_min,
                 const scalar &grid_cell_width) const override;

  virtual scalar weight(const unsigned int pidx, const MaterialPoints &points,
                        const Vector2u &nidx,
                        const PhysicsGrid &grid) const override;
  virtual scalar weight(const Vector2s &xp, const scalar &half_width,
                        const Vector2u &nidx,
                        const PhysicsGrid &grid) const override;
  virtual scalar weight(const Vector2s &xp, const scalar &half_width,
                        const Vector2u &nidx, const Vector2s &grid_min,
                        const scalar &grid_cell_width) const override;

  virtual Vector2s weightGrad(const unsigned int pidx,
                              const MaterialPoints &points,
                              const Vector2u &nidx,
                              const PhysicsGrid &grid) const override;
  virtual Vector2s weightGrad(const Vector2s &xp, const scalar &half_width,
                              const Vector2u &nidx,
                              const PhysicsGrid &grid) const override;
  virtual Vector2s weightGrad(const Vector2s &xp, const scalar &half_width,
                              const Vector2u &nidx, const Vector2s &grid_min,
                              const scalar &grid_cell_width) const override;
};

#endif

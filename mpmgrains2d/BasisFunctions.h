#ifndef BASIS_FUNCTIONS_2D_H
#define BASIS_FUNCTIONS_2D_H

#include "scisim/Math/MathDefines.h"
#include <iostream>
#include <memory>

struct PhysicsGrid;
struct MaterialPoints;

enum class BasisFunctionCategory { Standard, uGIMP };

enum class BasisFunctionType { Linear, ThirdOrder, uGIMPLinear };

struct BasisFunctions {

  BasisFunctions() = default;

  BasisFunctions(const BasisFunctions &) = delete;
  BasisFunctions(BasisFunctions &&) = delete;
  BasisFunctions &operator=(const BasisFunctions &) = delete;
  BasisFunctions &operator=(BasisFunctions &&) = delete;

  virtual ~BasisFunctions() = 0;

  virtual std::unique_ptr<BasisFunctions> clone() const = 0;

  virtual BasisFunctionType type() const = 0;

  virtual std::pair<Array2u, Array2u>
  computeStencil(const unsigned int pidx, const MaterialPoints &points,
                 const PhysicsGrid &grid) const = 0;
  virtual std::pair<Array2u, Array2u>
  computeStencil(const Vector2s &xp, const scalar &half_width,
                 const PhysicsGrid &grid) const = 0;
  virtual std::pair<Array2u, Array2u>
  computeStencil(const Vector2s &xp, const scalar &half_width,
                 const Vector2s &grid_min,
                 const scalar &grid_cell_width) const = 0;

  virtual scalar weight(const unsigned int pidx, const MaterialPoints &points,
                        const Vector2u &nidx,
                        const PhysicsGrid &grid) const = 0;
  virtual scalar weight(const Vector2s &xp, const scalar &half_width,
                        const Vector2u &nidx,
                        const PhysicsGrid &grid) const = 0;
  virtual scalar weight(const Vector2s &xp, const scalar &half_width,
                        const Vector2u &nidx, const Vector2s &grid_min,
                        const scalar &grid_cell_width) const = 0;

  virtual Vector2s weightGrad(const unsigned int pidx,
                              const MaterialPoints &points,
                              const Vector2u &nidx,
                              const PhysicsGrid &grid) const = 0;
  virtual Vector2s weightGrad(const Vector2s &xp, const scalar &half_width,
                              const Vector2u &nidx,
                              const PhysicsGrid &grid) const = 0;
  virtual Vector2s weightGrad(const Vector2s &xp, const scalar &half_width,
                              const Vector2u &nidx, const Vector2s &grid_min,
                              const scalar &grid_cell_width) const = 0;
};

namespace BasisFunctionsTools {
void serialize(const std::unique_ptr<BasisFunctions> &basis_functions,
               std::ostream &output_stream);
std::unique_ptr<BasisFunctions> deserialize(std::istream &input_stream);
} // namespace BasisFunctionsTools

#endif

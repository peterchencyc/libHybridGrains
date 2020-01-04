#ifndef GRID_H
#define GRID_H

#include "mpmgrains2d/BasisFunctions.h"
#include "scisim/Math/MathDefines.h"
#include <vector>

struct MaterialPoints;

class Grid final {

public:
  Grid(const Vector2s &min, const Vector2s &max, const Vector2i &cell_counts,
       const unsigned val_length);

  std::vector<bool> &cellOccupied();
  const std::vector<bool> &cellOccupied() const;

  const VectorXs &m() const;
  scalar &m(const int node_idx);
  const scalar &m(const int node_idx) const;
  scalar totalMass() const;

  MatrixXXsc &vals();

  // Interpolated value on the grid
  VectorXs interpolateValue(const std::unique_ptr<BasisFunctions> &bf,
                            const MaterialPoints &points,
                            const Vector2s &x) const;

  void rasterizeMass(const MaterialPoints &points, const int npoints,
                     const std::unique_ptr<BasisFunctions> &basis_funcs);
  void rasterizeMassWeightedValuesVector(
      const MaterialPoints &points, const int npoints,
      const std::unique_ptr<BasisFunctions> &basis_funcs,
      const VectorXs &values);
  void rasterizeMassWeightedValuesMatrix(
      const MaterialPoints &points, const int npoints,
      const std::unique_ptr<BasisFunctions> &basis_funcs,
      const MatrixXXsc &values);
  void rasterizeMassWeightedValues(
      const MaterialPoints &points, const int npoints,
      const std::unique_ptr<BasisFunctions> &basis_funcs,
      const std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>
          &values);
  void normalizationByMass();

  void normalizationByMass(const VectorXs &default_values_for_zero_mass);

  // Interpolated gradient on the grid
  MatrixXXsc interpolateGradient(const std::unique_ptr<BasisFunctions> &bf,
                                 const MaterialPoints &points,
                                 const Vector2s &x) const;

  const Vector2s &start() const;

  const scalar &cellWidth() const;

  Vector2i containingCell(const Vector2s &x) const;

  int nodeIndex(const int i_idx, const int j_idx) const;
  int cellIndex(const int i_idx, const int j_idx) const;

  // Cell index that contains the point x
  int cellIndex(const Vector2s &x) const;

  unsigned numNodes() const;

  const Vector2i &cellCount() const;

#ifndef NDEBUG
  bool nodeIndicesValid(const int xidx, const int yidx) const;
#endif

private:
  // Number of cells along each dimension
  const Vector2i m_N;

  // Region bounded by this grid
  const Vector2s m_min;
  const Vector2s m_max;

  // Width of grid cells
  const scalar m_h;

  // Velocity stored at each grid point
  MatrixXXsc m_vals;

  // Mass stored at each grid point
  VectorXs m_m;

  std::vector<bool> m_cell_occupied;

  const int m_val_length;
};

#endif

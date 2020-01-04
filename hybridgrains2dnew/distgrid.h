#ifndef DIST_GRID_H
#define DIST_GRID_H

#include "scisim/Math/MathDefines.h"
#include <Eigen/Dense>

struct DistGrid final {
  void resize(const Eigen::Vector2d &in_min_coords, const Vector2u &in_nCells,
              const double &in_h);

  Vector2u nCells() const;
  int num_points() const;
  int num_cells() const;

  Vector2u N;
  Eigen::Vector2d min_coords;
  Eigen::Vector2d max_coords;
  double h;

  Eigen::VectorXd minDist;
  Eigen::VectorXi closestPoint;
};

#endif

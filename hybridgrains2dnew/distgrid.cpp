#include "distgrid.h"

void DistGrid::resize(const Eigen::Vector2d &in_min_coords,
                      const Vector2u &in_nCells, const double &in_h) {
  min_coords = in_min_coords;
  N = in_nCells.array();
  h = in_h;
  max_coords = min_coords.array() + h * (in_nCells.array()).cast<double>();

  minDist.setConstant(num_points(), 0.0);
  closestPoint.setConstant(num_points(), 0.0);
}

Vector2u DistGrid::nCells() const { return N.array(); }

int DistGrid::num_points() const { return (N.x() + 1) * (N.y() + 1); }

int DistGrid::num_cells() const { return (N.x()) * (N.y()); }

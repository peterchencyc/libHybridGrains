#ifndef LEVEL_SET_UTIL_H
#define LEVEL_SET_UTIL_H

#include <cmath>
#include <map>

#include "distgrid.h"
#include <Eigen/Dense>

struct CrossingPoints {
  unsigned nPoints;
  unsigned nAlloc;
  Eigen::Matrix2Xd x;
};

struct SIDPair {
  int id1, id2;
};

inline bool operator<(const SIDPair &lhs, const SIDPair &rhs) {
  return (lhs.id1 < rhs.id1) || ((lhs.id1 == rhs.id1) && (lhs.id2 < rhs.id2));
}

inline bool operator>(const SIDPair &lhs, const SIDPair &rhs) {
  return (lhs.id1 > rhs.id1) || ((lhs.id1 == rhs.id1) && (lhs.id2 > rhs.id2));
}

inline double point2LineSegDist(const Eigen::Vector2d &x,
                                const Eigen::Vector2d &lineSegEnd1,
                                const Eigen::Vector2d &lineSegEnd2) {
  const double numer = (lineSegEnd1 - x).dot(lineSegEnd1 - lineSegEnd2);
  const double denom =
      (lineSegEnd1 - lineSegEnd2).dot(lineSegEnd1 - lineSegEnd2);
  const double s = numer / denom;

  if (s <= 0.0)
    return (x - lineSegEnd1).norm();
  else if (s >= 1.0)
    return (x - lineSegEnd2).norm();
  else {
    const Eigen::Vector2d h = s * lineSegEnd2 + (1.0 - s) * lineSegEnd1;
    return (x - h).norm();
  }
}

inline double point2CircleDist(const Eigen::Vector2d &x,
                               const Eigen::Vector2d &center, const double &r) {
  return fabs(r - (center - x).norm());
}

void initCrossingPoints(CrossingPoints &io_XPoints);
void reallocCrossingPoints(CrossingPoints &io_XPoints, unsigned nPoints);

// void genCrossPoints2(const DistGrid& grid, CrossingPoints& io_XPoints, const
// MaterialPoints& points);
void redistanceGrid(DistGrid &grid, const CrossingPoints &in_XPoints,
                    const double &particleSize);

#endif

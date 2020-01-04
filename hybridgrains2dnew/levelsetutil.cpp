#include "levelsetutil.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <queue>

void initCrossingPoints(CrossingPoints &io_XPoints) {
  io_XPoints.nPoints = 0;
  io_XPoints.nAlloc = 0;
  io_XPoints.x.resize(Eigen::NoChange, 0);
}

void reallocCrossingPoints(CrossingPoints &io_XPoints, unsigned nPoints) {
  if (io_XPoints.nAlloc >= nPoints)
    return;

  if (io_XPoints.nAlloc < 4096)
    io_XPoints.nAlloc = std::max(io_XPoints.nAlloc * 2, nPoints);
  else
    io_XPoints.nAlloc = std::max(io_XPoints.nAlloc + 1024, nPoints);

  io_XPoints.x.conservativeResize(Eigen::NoChange, io_XPoints.nAlloc);
}

void redistanceGrid(DistGrid &grid, const CrossingPoints &in_XPoints,
                    const double &particleSize) {
  double *absDist = (double *)malloc(sizeof(double) * grid.num_points());

  for (int i = 0; i < grid.num_points(); i++) {
    grid.closestPoint(i) = -1;
    absDist[i] = 1.0e10;
  }

  for (unsigned p = 0; p < in_XPoints.nPoints; p++) {
    double r = 3.0 * particleSize;
    Eigen::Vector2d min_x(in_XPoints.x(0, p) - r, in_XPoints.x(1, p) - r);
    Eigen::Vector2d max_x(in_XPoints.x(0, p) + r, in_XPoints.x(1, p) + r);

    const Eigen::Vector2i min_grid_idx(
        std::max<int>(
            0, std::min<int>(
                   grid.N.x(),
                   int(floor((grid.N.x()) * (min_x(0) - grid.min_coords(0)) /
                             (grid.max_coords(0) - grid.min_coords(0)))))),
        std::max<int>(
            0, std::min<int>(
                   grid.N.y(),
                   int(floor((grid.N.y()) * (min_x(1) - grid.min_coords(1)) /
                             (grid.max_coords(1) - grid.min_coords(1)))))));

    const Eigen::Vector2i max_grid_idx(
        std::max<int>(
            0, std::min<int>(
                   grid.N.x(),
                   int(ceil((grid.N.x()) * (max_x(0) - grid.min_coords(0)) /
                            (grid.max_coords(0) - grid.min_coords(0)))))),
        std::max<int>(
            0, std::min<int>(
                   grid.N.y(),
                   int(ceil((grid.N.y()) * (max_x(1) - grid.min_coords(1)) /
                            (grid.max_coords(1) - grid.min_coords(1)))))));

    for (int j = min_grid_idx.y(); j <= max_grid_idx.y(); j++) {
      for (int i = min_grid_idx.x(); i <= max_grid_idx.x(); i++) {
        const unsigned flat_idx = j * (grid.N.x() + 1) + i;

        const Eigen::Vector2d gp((grid.max_coords.x() - grid.min_coords.x()) *
                                         double(i) / double(grid.N.x()) +
                                     grid.min_coords.x(),
                                 (grid.max_coords.y() - grid.min_coords.y()) *
                                         double(j) / double(grid.N.y()) +
                                     grid.min_coords.y());

        const double dist = (gp - in_XPoints.x.col(p)).norm();

        if (dist < absDist[flat_idx]) {
          absDist[flat_idx] = dist;
          grid.closestPoint(flat_idx) = p;
        }
      }
    }
  }

  // build initial queue of cells with known closest points
  std::queue<unsigned> gridPointQueue;
  for (int i = 0; i < grid.num_points(); i++) {
    if (grid.closestPoint(i) >= 0)
      gridPointQueue.push(i);
  }

  // Repeat until everyone has pushed to all their neighbours
  // and no new changes have occurred.
  while (!gridPointQueue.empty()) {
    unsigned gpID = gridPointQueue.front();
    gridPointQueue.pop();

    int cp = grid.closestPoint(gpID);

    int neighborIDs[4];
    unsigned gci = gpID % (grid.N.x() + 1);
    unsigned gcj = gpID / (grid.N.x() + 1);

    neighborIDs[0] = (gci > 0) ? int(gcj * (grid.N.x() + 1) + gci - 1) : -1;
    neighborIDs[1] = (gcj > 0) ? int((gcj - 1) * (grid.N.x() + 1) + gci) : -1;
    neighborIDs[2] = (int(gci) < int(grid.N.x()))
                         ? int(gcj * (grid.N.x() + 1) + gci + 1)
                         : -1;
    neighborIDs[3] = (int(gcj) < int(grid.N.y()))
                         ? int((gcj + 1) * (grid.N.x() + 1) + gci)
                         : -1;

    // Check neighbours;
    for (int l = 0; l < 4; l++) {
      if (neighborIDs[l] >= 0) {
        int _i = neighborIDs[l] % (grid.N.x() + 1);
        int _j = neighborIDs[l] / (grid.N.x() + 1);

        const Eigen::Vector2d gp((grid.max_coords.x() - grid.min_coords.x()) *
                                         double(_i) / double(grid.N.x()) +
                                     grid.min_coords.x(),
                                 (grid.max_coords.y() - grid.min_coords.y()) *
                                         double(_j) / double(grid.N.y()) +
                                     grid.min_coords.y());

        const double distance = (gp - in_XPoints.x.col(cp)).norm();
        if ((distance < absDist[neighborIDs[l]]) ||
            (grid.closestPoint(neighborIDs[l]) < 0)) {
          absDist[neighborIDs[l]] = distance;
          grid.closestPoint(neighborIDs[l]) = cp;
          gridPointQueue.push(neighborIDs[l]);
        }
      }
    }
  }

  // Take existing signs for starters.
  for (int i = 0; i < grid.num_points(); i++) {
    // grid.minDist(i) = (grid.minDist(i) >= 0) ? absDist[i] : -absDist[i];
    grid.minDist(i) = absDist[i];
  }

  free(absDist);
}

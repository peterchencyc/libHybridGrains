// SpatialGridDetector.cpp
//
// Breannan Smith
// Last updated: 09/14/2015

#include "SpatialGrid.h"

#include "scisim/Utilities.h"

AABB::AABB(const Array2s &min, const Array2s &max) : m_min(min), m_max(max) {
  assert((m_min < m_max).all());
}

bool AABB::overlaps(const AABB &other) const {
  // Temporary sanity check: internal code shouldn't compare an AABB to itself
  assert(&other != this);
  // If no separating axis exists, the AABBs overlap
  return (m_max >= other.m_min).all() && (other.m_max >= m_min).all();
}

SpatialGrid::SpatialGrid() : m_h(0.0) {}

SpatialGrid::SpatialGrid(const scalar &h) : m_h(h) { assert(m_h > 0.0); }

void SpatialGrid::setCellWidth(const scalar &h) {
  assert(h > 0.0);
  m_h = h;
}

static void rasterizeAABBs(const std::vector<AABB> &aabbs,
                           const Array2s &min_coord, const scalar &h,
                           const Array2u &dimensions,
                           std::map<unsigned, std::vector<unsigned>> &voxels) {
  // For each bounding box
  for (std::vector<AABB>::size_type aabb_idx = 0; aabb_idx < aabbs.size();
       ++aabb_idx) {
    // Compute the cells the AABB overlaps with. Slightly enlarge the boxes to
    // account for FPA errors.
    const Array2u index_lower{
        ((aabbs[aabb_idx].min() - min_coord - 1.0e-6) / h).cast<unsigned>()};
    const Array2u index_upper{
        ((aabbs[aabb_idx].max() - min_coord + 1.0e-6) / h).cast<unsigned>()};
    assert((index_lower <= index_upper).all());

    for (unsigned x_idx = index_lower.x(); x_idx <= index_upper.x(); ++x_idx) {
      for (unsigned y_idx = index_lower.y(); y_idx <= index_upper.y();
           ++y_idx) {
        voxels[x_idx + dimensions.x() * y_idx].emplace_back(aabb_idx);
      }
    }
  }
}

static void initializeSpatialGrid(const std::vector<AABB> &aabbs,
                                  const scalar &h_scale, Array2s &min_coord,
                                  Array2u &dimensions, scalar &h) {
  // Compute a bounding box for all AABBs
  min_coord.setConstant(SCALAR_INFINITY);
  Array2s max_coord{-SCALAR_INFINITY, -SCALAR_INFINITY};
  for (const AABB &aabb : aabbs) {
    assert((aabb.min() < aabb.max()).all());
    min_coord = min_coord.min(aabb.min());
    max_coord = max_coord.max(aabb.max());
    assert((min_coord < max_coord).all());
  }
  // Inflate the AABB to account for FPA quantization errors
  min_coord -= 2.0e-6;
  max_coord += 2.0e-6;

  // Compute the grid cell width
  {
    Array2s delta{0.0, 0.0};
    for (const AABB &aabb : aabbs) {
      delta += aabb.max() - aabb.min();
    }
    h = delta.maxCoeff() / scalar(aabbs.size());
  }

  h *= h_scale;

  // Compute the number of cells in the grid
  dimensions = ((max_coord - min_coord) / h)
                   .unaryExpr([](const scalar &s) {
                     using std::ceil;
                     return ceil(s);
                   })
                   .cast<unsigned>();
}

void SpatialGrid::getPotentialOverlaps(
    const std::vector<AABB> &aabbs,
    std::set<std::pair<unsigned, unsigned>> &overlaps) const {
  Array2s min_coord;
  Array2u dimensions;
  scalar h;
  initializeSpatialGrid(aabbs, m_h, min_coord, dimensions, h);

  std::map<unsigned, std::vector<unsigned>> voxels;
  rasterizeAABBs(aabbs, min_coord, h, dimensions, voxels);

  // For each voxel
  for (const auto &voxel : voxels) {
    // Visit each pair of AABBs in this voxel
    for (std::vector<unsigned>::size_type idx0 = 0;
         idx0 + 1 < voxel.second.size(); ++idx0) {
      for (std::vector<unsigned>::size_type idx1 = idx0 + 1;
           idx1 < voxel.second.size(); ++idx1) {
        // If the AABBs overlap
        if (aabbs[voxel.second[idx0]].overlaps(aabbs[voxel.second[idx1]])) {
          assert(voxel.second[idx1] < aabbs.size());
          assert(voxel.second[idx0] < voxel.second[idx1]);
          overlaps.insert(
              std::make_pair(voxel.second[idx0], voxel.second[idx1]));
        }
      }
    }
  }
}

void SpatialGrid::getPotentialOverlaps(const AABB &trial_aabb,
                                       const std::vector<AABB> &aabbs,
                                       std::vector<unsigned> &overlaps) const {
  // Compare the trial AABB to each input AABB
  for (unsigned aabb_idx = 0; aabb_idx < aabbs.size(); ++aabb_idx) {
    if (aabbs[aabb_idx].overlaps(trial_aabb)) {
      overlaps.emplace_back(aabb_idx);
    }
  }
}

void SpatialGrid::serialize(std::ostream &output_stream) const {
  Utilities::serializeBuiltInType(m_h, output_stream);
}

void SpatialGrid::deserialize(std::istream &input_stream) {
  m_h = Utilities::deserialize<scalar>(input_stream);
}

void SpatialGrid::setCellScale(const scalar &h) { m_h = h; }

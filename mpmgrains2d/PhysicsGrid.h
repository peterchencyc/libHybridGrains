#ifndef PHYSICS_GRID_2D_H
#define PHYSICS_GRID_2D_H

#include "scisim/Math/MathDefines.h"

#include <memory>
#include <vector>

struct MaterialPoints;
struct BasisFunctions;
struct MPMStaticPlane;

struct PhysicsGrid final {

  void setDimensions(const Array2s &grid_min, const Array2s &grid_max,
                     const scalar &cell_w);

  void clear();
  void clearRasterizedData();

  // for implicit coupling
  void clearCouplingForce();

  void serialize(std::ostream &strm) const;
  void deserialize(std::istream &input_stream);

  Vector2s computeTotalMomentum() const;
  scalar computeTotalAngularMomentum() const;

  void rasterizePointMasses(const MaterialPoints &points,
                            const std::unique_ptr<BasisFunctions> &basis_funcs);
  void
  rasterizePointMomentum(const MaterialPoints &points,
                         const std::unique_ptr<BasisFunctions> &basis_funcs);

  void computeForces(const MaterialPoints &points,
                     const std::unique_ptr<BasisFunctions> &basis_funcs,
                     const Vector2s &near_earth_gravity);
  void updateMomentum(const scalar &dt);
  void lumpedMassVelocityAndAccelerationUpdate(const scalar &dt);

  // Total number of grid *points*
  unsigned numGridPoints() const;

  // Total number of grid *cells*
  unsigned numGridCells() const;

  // Given a flat node index, computes the triplet index
  Vector2u computeNodeTripletIndex(const unsigned flat_idx) const;

  // Computes the location of a grid point
  // Vector2s computeGridPointLocation( const unsigned xidx, const unsigned yidx
  // ) const;
  Vector2s computeGridPointLocation(const unsigned flat_idx) const;

  // Given the x, y, and z index of a grid node, returns the index in linear
  // storage
  unsigned flatNodeIndex(const unsigned xidx, const unsigned yidx) const;

  // #ifndef NDEBUG
  bool nodeIndicesValid(const unsigned xidx, const unsigned yidx) const;
  // #endif

  void
  resolvePlaneCollisions(const std::vector<MPMStaticPlane> &plane_obstacles);

  // Lower left corner of the grid
  Array2s min;
  // Upper right corner of the grid
  Array2s max;
  // Width a grid cell
  scalar cell_width;
  // Number of grid cells along each dimension
  Array2u cell_count;

  // Data rasterized from material points
  VectorXs rasterized_mass;
  Matrix2Xsc rasterized_momentum;
  // Post-force momentum
  Matrix2Xsc rasterized_momentum_new;

  // Force computed on the grid
  Matrix2Xsc force;

  // Integration results
  Matrix2Xsc velocity;
  Matrix2Xsc acceleration;

  // for implicit coupling
  Matrix2Xsc coupling_force;

  VectorXs homog_factor;

  VectorXi hybridized;
};

#endif

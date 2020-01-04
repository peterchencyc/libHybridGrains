#ifndef INITIAL_SIMULATION_STATE_2D_H
#define INITIAL_SIMULATION_STATE_2D_H

#include <Eigen/StdVector>

#include "mpmgrains2d/SimulationState.h"
#include "mpmgrains2d/StaticPlane.h"
#include "scisim/Math/MathUtilities.h"
#include "scisim/Math/Rational.h"
#include "scisim/Utilities.h"

struct RectangleRegion final {
  RectangleRegion(const Array2s &region_min, const Array2s &region_max,
                  const scalar &density, const unsigned num_points_per_dim,
                  const Vector2s &vel, const scalar &angular_vel);

  bool containsPoint(const Vector2s &p) const;

  scalar volume() const;

  Vector2s computeCenter() const;

  Array2s min;
  Array2s max;
  scalar rho;
  unsigned npoints_per_dim;
  Vector2s v;
  scalar omega;
};

struct GridSettings final {
  Array2s min;
  Array2s max;
  scalar cell_width;
};

struct InitialSimulationState final {
  SimulationState generateSimulationState() const;

  // Initial material point settings
  std::vector<RectangleRegion> rectangle_regions;

  // Background grid settings
  GridSettings grid_settings;

  // External obstacles
  std::vector<MPMStaticPlane> static_planes;

  // Material settings
  scalar shear_modulus;
  scalar bulk_modulus;
  scalar Drucker_Prager_alpha;

  // Integrator settings
  std::string dt_string;
  Rational<std::intmax_t> dt;
  scalar end_time;
  // alpha == 0 => pic, alpha == 1 => flip
  scalar alpha;
  // Order of basis functions
  unsigned basis_function_order;
  // Category of basis functions;
  BasisFunctionCategory basis_function_category;

  // near_earth_gravity
  Vector2s near_earth_gravity;
};

#endif

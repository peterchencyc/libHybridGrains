#include "InitialSimulationState.h"

#include "CubicBasisFunctions.h"
#include "LinearBasisFunctions.h"
#include "uGIMPLinearBasisFunctions.h"

#include <iostream>

RectangleRegion::RectangleRegion(const Array2s &region_min,
                                 const Array2s &region_max,
                                 const scalar &density,
                                 const unsigned num_points_per_dim,
                                 const Vector2s &vel, const scalar &angular_vel)
    : min(region_min), max(region_max), rho(density),
      npoints_per_dim(num_points_per_dim), v(vel), omega(angular_vel) {}

bool RectangleRegion::containsPoint(const Vector2s &p) const {
  if ((p.array() < min).any()) {
    return false;
  }
  if ((p.array() > max).any()) {
    return false;
  }
  return true;
}

scalar RectangleRegion::volume() const { return (max - min).prod(); }

Vector2s RectangleRegion::computeCenter() const {
  return 0.5 * (min + max).matrix();
}

static bool pointInsidePlane(const std::vector<MPMStaticPlane> &static_planes,
                             const Vector2s &p) {
  for (const MPMStaticPlane &plane : static_planes) {
    if (p(1) < plane.lower_bound) {
      return false;
    }
    if (plane.distanceToPoint(p) < 0.0) {
      return true;
    }
  }
  return false;
}

static std::vector<InitialMaterialPoint,
                   Eigen::aligned_allocator<InitialMaterialPoint>>
rasterizeRectangle(const RectangleRegion &region,
                   const PhysicsGrid &physics_grid,
                   const std::vector<MPMStaticPlane> &static_planes) {
  // TODO: Move these checks into the parser
  // Ensure that the rectangle is inside the simulated region
  if ((region.min < physics_grid.min).any()) {
    std::cerr << "Error, initial rectangle region is below the simulation "
                 "grid. Exiting."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((region.max > physics_grid.max).any()) {
    std::cerr << "Error, initial rectangle region is above the simulation "
                 "grid. Exiting."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::vector<InitialMaterialPoint,
              Eigen::aligned_allocator<InitialMaterialPoint>>
      material_points;

  // Which grid cells overlap the region to rasterize
  using std::floor;
  const Array2u grid_bounds_start =
      ((region.min - physics_grid.min) / physics_grid.cell_width)
          .unaryExpr([](const scalar &s) { return floor(s); })
          .cast<unsigned>();
  using std::ceil;
  const Array2u grid_bounds_end =
      ((region.max - physics_grid.min) / physics_grid.cell_width)
          .unaryExpr([](const scalar &s) { return ceil(s); })
          .cast<unsigned>();

  unsigned long num_inserted_points{0};

  const Vector2s rectangle_center{region.computeCenter()};

  // For each physics grid cell in in the rectangle region's bounding box
  for (unsigned gridy = grid_bounds_start.y(); gridy < grid_bounds_end.y();
       gridy++) {
    for (unsigned gridx = grid_bounds_start.x(); gridx < grid_bounds_end.x();
         gridx++) {
      // Generate evenly spaced samples
      for (unsigned celly = 0; celly < region.npoints_per_dim; celly++) {
        for (unsigned cellx = 0; cellx < region.npoints_per_dim; cellx++) {
          // Compute the world space position of this sample
          const Vector2s position{
              physics_grid.min +
              physics_grid.cell_width *
                  (Array2s(gridx, gridy) + (Array2s(cellx, celly) + 0.5) /
                                               scalar(region.npoints_per_dim))};

          // Skip points that are outside the actual rectangle
          if (!region.containsPoint(position)) {
            continue;
          }

          // Skip points that intersect a boundary
          if (pointInsidePlane(static_planes, position)) {
            continue;
          }

          const Matrix22sc be_bar{Matrix22sc::Identity()};
          const scalar J{1.0};

          const Vector3s rot_axis{0.0, 0.0, 1.0};
          Vector3s arm;
          arm << position - rectangle_center, 0.0;

          const Vector2s velocity{
              region.v + region.omega * (rot_axis.cross(arm)).segment<2>(0)};

          material_points.push_back(
              InitialMaterialPoint(position, be_bar, J, velocity));

          num_inserted_points++;
        }
      }
    }
  }

  // Set the masses
  const scalar hl =
      0.5 * physics_grid.cell_width / scalar(region.npoints_per_dim);
  const scalar point_volume = hl * hl * 4.0;
  const scalar point_mass = region.rho * point_volume;

  for (InitialMaterialPoint &material_point : material_points) {
    material_point.m = point_mass;
    material_point.volume = point_volume;
    material_point.hl = hl;
    material_point.always_fixed = false;
    material_point.kinematically_scripted = false;
  }

  return material_points;
}

SimulationState InitialSimulationState::generateSimulationState() const {
  SimulationState new_state;

  new_state.static_planes = static_planes;
  new_state.physics_grid.setDimensions(grid_settings.min, grid_settings.max,
                                       grid_settings.cell_width);

  new_state.shear_modulus = shear_modulus;
  new_state.bulk_modulus = bulk_modulus;

  new_state.Drucker_Prager_alpha = Drucker_Prager_alpha;

  new_state.near_earth_gravity = near_earth_gravity;

  new_state.alpha = alpha;

  if (basis_function_category == BasisFunctionCategory::Standard) {
    if (basis_function_order == 1) {
      new_state.basis_functions = std::make_unique<LinearBasisFunctions>();
    } else if (basis_function_order == 3) {
      new_state.basis_functions = std::make_unique<ThirdOrderBasisFunctions>();
    } else {
      std::cerr << "Error in InitialSimulationState::generateSimulationState, "
                   "standard basis functions are not first or third order."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else if (basis_function_category == BasisFunctionCategory::uGIMP) {
    if (basis_function_order == 1) {
      new_state.basis_functions = std::make_unique<uGIMPLinearBasisFunctions>();
    } else {
      std::cerr << "Error in InitialSimulationState::generateSimulationState, "
                   "uGIMP basis functions are not first order."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    std::cerr << "Error in InitialSimulationState::generateSimulationState, "
                 "unsupported basis function category."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::vector<InitialMaterialPoint,
              Eigen::aligned_allocator<InitialMaterialPoint>>
      new_points;
  for (const RectangleRegion &rectangle : rectangle_regions) {
    const std::vector<InitialMaterialPoint,
                      Eigen::aligned_allocator<InitialMaterialPoint>>
        rectangle_points = rasterizeRectangle(rectangle, new_state.physics_grid,
                                              static_planes);
    new_points.insert(new_points.end(), rectangle_points.begin(),
                      rectangle_points.end());
  }
  new_state.material_points.setPoints(new_points);

  if (rectangle_regions.size() == 0) {
    std::cerr << "Provide at least one rectangle region." << std::endl;
    std::exit(EXIT_FAILURE);
  } else {
    // Need a better to consistently specify these variables...

    if (rectangle_regions.size() > 0) {
      new_state.initial_particle_size =
          new_state.physics_grid.cell_width /
          scalar(rectangle_regions[0].npoints_per_dim);
      new_state.material_density = rectangle_regions[0].rho;
      std::cout << "set initial_particle_size from cuboid data: "
                << new_state.initial_particle_size << std::endl;
      std::cout << "set material_density from cuboid data: "
                << new_state.material_density << std::endl;
    }

    std::cout << "new_state.initial_particle_size: "
              << new_state.initial_particle_size << std::endl;

    new_state.initial_particle_volume =
        new_state.initial_particle_size * new_state.initial_particle_size;
  }

  return new_state;
}

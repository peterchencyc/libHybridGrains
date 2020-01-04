#ifndef SIMULATION_STATE_2D_H
#define SIMULATION_STATE_2D_H

#include <Eigen/StdVector>

#include "MaterialPoints.h"
#include "PhysicsGrid.h"
#include "StaticPlane.h"

#include <iosfwd>

#ifdef USE_HDF5
class HDF5File;
#endif

struct SimulationState final {

  SimulationState() = default;

  SimulationState(const SimulationState &other);
  SimulationState &operator=(const SimulationState &other);

  SimulationState(SimulationState &&) = default;
  SimulationState &operator=(SimulationState &&) = default;

  void clear();

  void serialize(std::ostream &strm) const;
  void deserialize(std::istream &input_stream);

#ifdef USE_HDF5
  void writeBinaryState(const std::string &prefix, HDF5File &output_file) const;
#endif

#ifndef NDEBUG
  bool isConsistent() const;
#endif

  void computeBoundingBox(Vector4s &bbox) const;

  MaterialPoints material_points;
  PhysicsGrid physics_grid;

  std::vector<MPMStaticPlane> static_planes;

  // Hyper-elastic material parameters
  scalar shear_modulus;
  scalar bulk_modulus;

  // Plasticity
  scalar Drucker_Prager_alpha;

  // near_earth_gravity
  Vector2s near_earth_gravity;

  // PIC/FLIP interpolation constant
  scalar alpha;

  // Shape functions
  std::unique_ptr<BasisFunctions> basis_functions;

  // For Hybrid resampling
  // Probably need a better and more consistent way for handling these
  // variables.
  scalar initial_particle_volume;
  scalar material_density;
  scalar initial_particle_size;
};

#endif

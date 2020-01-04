#include "SimulationState.h"

#include "scisim/Math/MathUtilities.h"
#include "scisim/Utilities.h"

#ifdef USE_HDF5
#include "scisim/HDF5File.h"
#endif

SimulationState::SimulationState(const SimulationState &other)
    : material_points(other.material_points), physics_grid(other.physics_grid),
      static_planes(other.static_planes), shear_modulus(other.shear_modulus),
      bulk_modulus(other.bulk_modulus),
      Drucker_Prager_alpha(other.Drucker_Prager_alpha),
      near_earth_gravity(other.near_earth_gravity), alpha(other.alpha),
      basis_functions(other.basis_functions->clone()),
      initial_particle_volume(other.initial_particle_volume),
      material_density(other.material_density),
      initial_particle_size(other.initial_particle_size) {}

SimulationState &SimulationState::operator=(const SimulationState &other) {
  SimulationState copy{other};
  using std::swap;
  swap(*this, copy);
  return *this;
}

#ifndef NDEBUG
template <typename Derived>
bool entriesAreValid(const Eigen::DenseBase<Derived> &eigen_variable) {
  for (int row = 0; row < eigen_variable.rows(); row++) {
    for (int col = 0; col < eigen_variable.cols(); col++) {
      using std::isnan;
      if (isnan(eigen_variable(row, col))) {
        return false;
      }
      using std::isinf;
      if (isinf(eigen_variable(row, col))) {
        return false;
      }
    }
  }
  return true;
}

template <typename T> bool entriesAreValid(const std::vector<T> &vec) {
  using std::all_of;
  return all_of(vec.begin(), vec.end(),
                [](const T &entry) { return entriesAreValid(entry); });
}

bool SimulationState::isConsistent() const {
  // Check that all material points are within the grid bounds
  for (unsigned pnt_idx = 0; pnt_idx < material_points.npoints; pnt_idx++) {
    const Array2s point_location{material_points.q.col(pnt_idx)};
    if ((point_location < physics_grid.min).any()) {
      return false;
    }
    if ((point_location > physics_grid.max).any()) {
      return false;
    }
  }

  // Check that the material points do not contain nans or infs
  if (!entriesAreValid(material_points.q)) {
    return false;
  }
  if (!entriesAreValid(material_points.v)) {
    return false;
  }
  if (!entriesAreValid(material_points.m)) {
    return false;
  }
  if (!entriesAreValid(material_points.volume)) {
    return false;
  }
  if (!entriesAreValid(material_points.J)) {
    return false;
  }
  // Material point volumes should be positive
  if ((material_points.volume.array() <= 0.0).any()) {
    return false;
  }

  // Volume preserving left Cauchy-Green should have a determinant of 1...
  using std::all_of;
  if (!all_of(material_points.be_bar.begin(), material_points.be_bar.end(),
              [](const Matrix22sc &b) {
                return fabs(b.determinant() - 1.0) <= 1.0e-6;
              })) {
    return false;
  }
  // ...and should be symmetric
  if (!all_of(material_points.be_bar.begin(), material_points.be_bar.end(),
              [](const Matrix22sc &b) {
                return (b - b.transpose()).lpNorm<Eigen::Infinity>() <= 1.0e-6;
              })) {
    return false;
  }

  // The stress should be symmetric
  if (!all_of(material_points.sigma.begin(), material_points.sigma.end(),
              [](const Matrix22sc &s) {
                return (s - s.transpose()).lpNorm<Eigen::Infinity>() <= 1.0e-6;
              })) {
    return false;
  }

  return true;
}
#endif

void SimulationState::clear() {
  material_points.clear();
  physics_grid.clear();
  static_planes.clear();
}

void SimulationState::serialize(std::ostream &strm) const {
  material_points.serialize(strm);
  physics_grid.serialize(strm);
  Utilities::serializeVectorCustomType(static_planes, strm);
  Utilities::serializeBuiltInType(shear_modulus, strm);
  Utilities::serializeBuiltInType(bulk_modulus, strm);
  Utilities::serializeBuiltInType(Drucker_Prager_alpha, strm);
  MathUtilities::serialize(near_earth_gravity, strm);
  Utilities::serializeBuiltInType(alpha, strm);
  BasisFunctionsTools::serialize(basis_functions, strm);
}

void SimulationState::deserialize(std::istream &input_stream) {
  material_points.deserialize(input_stream);
  physics_grid.deserialize(input_stream);
  Utilities::deserializeVectorCustomType<MPMStaticPlane>(static_planes,
                                                         input_stream);
  shear_modulus = Utilities::deserialize<scalar>(input_stream);
  bulk_modulus = Utilities::deserialize<scalar>(input_stream);
  Drucker_Prager_alpha = Utilities::deserialize<scalar>(input_stream);
  near_earth_gravity = MathUtilities::deserialize<Vector2s>(input_stream);
  alpha = Utilities::deserialize<scalar>(input_stream);
  basis_functions = BasisFunctionsTools::deserialize(input_stream);
}

#ifdef USE_HDF5
void SimulationState::writeBinaryState(const std::string &prefix,
                                       HDF5File &output_file) const {
  // Save out the number of material points
  output_file.writeScalar(prefix, "npoints", material_points.npoints);
  // Output the configuration
  output_file.writeMatrix(prefix, "q", material_points.q);
  // Output the velocity
  output_file.writeMatrix(prefix, "v", material_points.v);
  // Output the mass
  output_file.writeMatrix(prefix, "m", material_points.m);

  // Output the size
  output_file.writeMatrix(prefix, "r", material_points.hl);

  // Output the determinant
  output_file.writeMatrix(prefix, "J", material_points.J);

  // Output the 'just homogenized flag'
  output_file.writeMatrix(prefix, "jhf", material_points.just_homogenized);

  // Output the grid point locations, masses, and homogenization factors
  {
    VectorXs grid_points(2 * physics_grid.numGridPoints());
    for (int node_idx = 0; node_idx < int(physics_grid.numGridPoints());
         node_idx++) {
      grid_points.segment<2>(2 * node_idx) =
          physics_grid.computeGridPointLocation(node_idx);
    }
    output_file.writeMatrix(prefix, "grid_node_positions", grid_points);
    output_file.writeMatrix(prefix, "grid_node_masses",
                            physics_grid.rasterized_mass);
    output_file.writeMatrix(prefix, "grid_node_hybrid_factors",
                            physics_grid.homog_factor);
  }

  // Output the pressure
  VectorXs p;
  p.resize(material_points.q.cols());
  for (unsigned pnt_idx = 0; pnt_idx < material_points.npoints; pnt_idx++) {
    const scalar pressure = -material_points.sigma[pnt_idx].trace() / 2.0;
    p(pnt_idx) = pressure;
  }
  output_file.writeMatrix(prefix, "p", p);
}
#endif

void SimulationState::computeBoundingBox(Vector4s &bbox) const {
  bbox(0) = physics_grid.min.x();
  bbox(1) = physics_grid.max.x();
  bbox(2) = physics_grid.min.y();
  bbox(3) = physics_grid.max.y();
}

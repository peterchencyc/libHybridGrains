#include "ResamplingTools.h"
#include "mpmgrains2d/ConstitutiveModel.h"
#include "scisim/Timer/TimeUtils.h"
#include "uniformgrid.h"

void ResamplingTools::init(const scalar &dem_r_mean, const scalar &dem_r_std) {
  m_Vel_field = nullptr;
  m_J_field = nullptr;
  m_Bebar_field = nullptr;

  m_dem_r_mean = dem_r_mean;
  m_dem_r_std = dem_r_std;
  m_dist = std::normal_distribution<double>(dem_r_mean, dem_r_std);
}

double ResamplingTools::constantRadiusSampler(const double &r) { return r; }

double ResamplingTools::randomRadiusSampler() {
  double r = -1.0e33;
  while ((r < m_dem_r_mean - 2.0 * m_dem_r_std) ||
         (r > m_dem_r_mean + 2.0 * m_dem_r_std)) {
    r = m_dist(m_generator);
  }
  return r;
}

void ResamplingTools::deleteMPMparticlesUsingZoneIndicators(
    ZoneTools &zone_tools, SimulationState &continuum_state,
    std::vector<int> &mpm_old_to_new_id_map) {
  std::vector<unsigned> bodies_to_del;

  zone_tools.determineMPMParticlesToDelete(continuum_state, bodies_to_del);
  continuum_state.material_points.deleteMaterialPoint(bodies_to_del,
                                                      mpm_old_to_new_id_map);
}

void ResamplingTools::deleteDEMparticlesUsingZoneIndicators(
    ZoneTools &zone_tools, RigidBody2DState &discrete_state,
    std::vector<int> &dem_old_to_new_id_map, ImpactFrictionMap *ifmap) {
  std::vector<unsigned> bodies_to_del;

  zone_tools.determineDEMGrainsToDelete(discrete_state, bodies_to_del);
  discrete_state.removeBodies(bodies_to_del, dem_old_to_new_id_map, ifmap);
}

void ResamplingTools::deleteParticlesUsingZoneIndicators(
    ZoneTools &zone_tools, RigidBody2DState &discrete_state,
    SimulationState &continuum_state, ImpactFrictionMap *ifmap) {
  std::vector<int> mpm_old_to_new_id_map;
  std::vector<int> dem_old_to_new_id_map;

  deleteMPMparticlesUsingZoneIndicators(zone_tools, continuum_state,
                                        mpm_old_to_new_id_map);
  deleteDEMparticlesUsingZoneIndicators(zone_tools, discrete_state,
                                        dem_old_to_new_id_map, ifmap);
}

void ResamplingTools::avoidAVoidMPMparticlesInsideZonesWithTentativeQuantities(
    SimulationState &continuum_state,
    const std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
    int &numPointsBeforeSampling, int &numPointsAfterSampling) {
  // sample particles
  const double r_mpm = 0.5 * continuum_state.initial_particle_size;
  numPointsBeforeSampling = continuum_state.material_points.npoints;

  Matrix2Xsc pos;
  pos.resize(2, numPointsBeforeSampling);
  VectorXs r;
  r.resize(numPointsBeforeSampling);

  for (int pnt_idx = 0; pnt_idx < numPointsBeforeSampling; ++pnt_idx) {
    pos.col(pnt_idx) = continuum_state.material_points.q.col(pnt_idx);
    r(pnt_idx) = r_mpm;
  }

  m_avoid_a_void.setTargetRegions(zones);
  m_avoid_a_void.setCheckConflictAgainstRegionBoundary(false);
  m_avoid_a_void.setOldBodies(pos, r, r_mpm);
  m_avoid_a_void.avoidAVoid(
      std::bind(&ResamplingTools::constantRadiusSampler, this, r_mpm));
  m_avoid_a_void.getBodies(pos, r);

  numPointsAfterSampling = pos.cols();

  // assign properties
  for (int pnt_idx = numPointsBeforeSampling; pnt_idx < numPointsAfterSampling;
       ++pnt_idx) {
    const scalar volume = 0.0; // continuum_state.initial_particle_volume;
    const scalar hl = continuum_state.initial_particle_size * 0.5;
    const scalar mass = 0.0; // continuum_state.initial_particle_volume *
                             // continuum_state.material_density;

    const Vector2s x = pos.col(pnt_idx);
    const Vector2s x0 = x;
    const Vector2s v = Vector2s{0.0, 0.0};

    scalar J = 1.0;

    const Matrix22sc be_bar = Matrix22sc::Identity();

    continuum_state.material_points.addMaterialPoint(x0, x, v, mass, volume, hl,
                                                     be_bar, J, false);
  }
}

void ResamplingTools::avoidAVoidDEMparticlesInsideZonesWithTentativeQuantities(
    RigidBody2DState &discrete_state,
    const std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
    const scalar &in_rho_dem, ImpactFrictionMap *ifmap,
    int &numGrainsBeforeSampling, int &numGrainsAfterSampling) {
  // sample particles
  Matrix2Xsc pos;
  VectorXs r;
  VectorXu bdy_idx;
  discrete_state.getAllCircleBodies(pos, r, bdy_idx);

  int internal_numGrainsBeforeSampling = 0;
  int internal_numGrainsAfterSampling = 0;

  internal_numGrainsBeforeSampling = pos.cols();
  numGrainsBeforeSampling = discrete_state.numBodies();

  scalar r_mean = m_dem_r_mean;

  if (m_dem_r_mean < 0.0) {
    std::cout << "r_mean not set, using default sampler for radius."
              << std::endl;
    const double rep_r = r(0);
    m_dist = std::normal_distribution<double>(rep_r, 0.0);
    r_mean = rep_r;
  }

  // TimingTools timing_tools;

  // std::cout << "DEM setTargetRegions..." << std::endl;
  // timing_tools.start();
  m_avoid_a_void.setTargetRegions(zones);
  // timing_tools.stop("");

  // std::cout << "DEM setCheckConflictAgainstRegionBoundary..." << std::endl;
  // timing_tools.start();
  m_avoid_a_void.setCheckConflictAgainstRegionBoundary(false);
  // timing_tools.stop("");

  // std::cout << "DEM setOldBodies..." << std::endl;
  // timing_tools.start();
  m_avoid_a_void.setOldBodies(pos, r, r_mean);
  // timing_tools.stop("");

  // std::cout << "DEM avoidAVoid..." << std::endl;
  // timing_tools.start();
  m_avoid_a_void.avoidAVoid(
      std::bind(&ResamplingTools::randomRadiusSampler, this));
  // timing_tools.stop("");

  // std::cout << "DEM getBodies..." << std::endl;
  // timing_tools.start();
  m_avoid_a_void.getBodies(pos, r);
  // timing_tools.stop("");

  // std::cout << "DEM assign properties..." << std::endl;
  // timing_tools.start();
  internal_numGrainsAfterSampling = pos.cols();

  std::vector<Vector2s> x;
  std::vector<Vector2s> x0;
  std::vector<scalar> theta;
  std::vector<Vector2s> v;
  std::vector<scalar> omega;
  std::vector<scalar> rho;
  std::vector<unsigned> geo_idx;
  std::vector<bool> fixed;

  // assign properties
  for (int pnt_idx = internal_numGrainsBeforeSampling;
       pnt_idx < internal_numGrainsAfterSampling; ++pnt_idx) {
    const Vector2s x_dem{pos.col(pnt_idx)};
    const Vector2s x0_dem{pos.col(pnt_idx)};
    Vector2s v_dem = Vector2s{0.0, 0.0};
    const int geo_idx_dem = discrete_state.addCircleGeometry(r(pnt_idx));

    x.emplace_back(x_dem);
    x0.emplace_back(x0_dem);
    theta.emplace_back(0.0);
    v.emplace_back(v_dem);
    omega.emplace_back(0.0);
    rho.emplace_back(in_rho_dem);
    geo_idx.emplace_back(geo_idx_dem);
    fixed.emplace_back(false);
  }

  discrete_state.addBodies(x, x0, theta, v, omega, rho, geo_idx, fixed, ifmap);
  numGrainsAfterSampling = discrete_state.numBodies();

  // timing_tools.stop("");
}

static void gatherMassAndVolume(SimulationState &continuum_state,
                                const CUniformGridH2D &io_Ugrid,
                                const int in_InsertedParticleId) {
  const scalar sigma = 0.5 * continuum_state.initial_particle_size;

  // const int the_TestID = (*io_TestID);
  const scalar r = 2.0 * continuum_state.initial_particle_size;
  const Vector2s x =
      continuum_state.material_points.q.col(in_InsertedParticleId);

  const Vector2s min_coords(x(0) - r, x(1) - r);
  const Vector2s max_coords(x(0) + r, x(1) + r);

  Vector2u min_grid_idx;
  Vector2u max_grid_idx;
  io_Ugrid.getGridIDRange(min_coords, max_coords, min_grid_idx, max_grid_idx);

  std::vector<int> the_NeighborParticles;

  // scan the particles
  scalar totMass = 0.0;
  scalar totVolume = 0.0;
  scalar weight_sum = 0.0;
  Vector2s q0_sum = Vector2s(0.0, 0.0);

  for (unsigned j = min_grid_idx.y(); j <= max_grid_idx.y(); j++) {
    for (unsigned i = min_grid_idx.x(); i <= max_grid_idx.x(); i++) {
      const Vector2u grid_idx(i, j);
      const int *ids;
      int nElems;
      io_Ugrid.getIDs(grid_idx, &nElems, &ids);

      for (int l = 0; l < nElems; l++) {
        // not needed if r is const
        // if(points.mailbox(ids[l]) == the_TestID) continue;
        // points.mailbox(ids[l]) = the_TestID;

        if (ids[l] >= in_InsertedParticleId)
          continue;

        const Vector2s _x = continuum_state.material_points.q.col(ids[l]);
        const scalar dist2 = (x - _x).squaredNorm();
        if (dist2 < r * r) {
          scalar weight = exp(-(dist2 * dist2) / (2.0 * sigma * sigma));
          weight_sum += weight;
          Vector2s q0 = continuum_state.material_points.q0.col(ids[l]);
          q0_sum += weight * q0;
          scalar _mass = continuum_state.material_points.m0(ids[l]);
          scalar _volume = continuum_state.material_points.volume0(ids[l]);
          the_NeighborParticles.push_back(ids[l]);

          totMass += _mass;
          totVolume += _volume;
        }
      }
    }
  }

  int count = the_NeighborParticles.size();
  if (count == 0) {
    std::cout << "WARNING: material point " << in_InsertedParticleId
              << " inserted in the inner continuum region does not have a "
                 "neighboring material point to gather mass..."
              << std::endl;
    continuum_state.material_points.q0.col(in_InsertedParticleId) = x;
  }

  continuum_state.material_points.q0.col(in_InsertedParticleId) =
      q0_sum / weight_sum;

  const scalar targetMass = totMass / (count + 1.0);
  continuum_state.material_points.m0(in_InsertedParticleId) = targetMass;
  const scalar targetVolume = totVolume / (count + 1.0);
  continuum_state.material_points.volume0(in_InsertedParticleId) = targetVolume;

  scalar factor = scalar(count) / scalar(count + 1.0);

  // redistribute mass and volume
  for (unsigned i = 0; i < the_NeighborParticles.size(); i++) {
    continuum_state.material_points.m0(the_NeighborParticles[i]) *= factor;
    continuum_state.material_points.volume0(the_NeighborParticles[i]) *= factor;
  }
}

void ResamplingTools::computeFieldsForInterpolation(
    SimulationState &continuum_state, const int numPointsBeforeSampling) {
  std::cout << "preparing the interpolation field..." << std::endl;
  if (m_J_field == nullptr) {
    m_J_field = new Grid(
        continuum_state.physics_grid.min, continuum_state.physics_grid.max,
        continuum_state.physics_grid.cell_count.cast<int>(), 1);
  }
  if (m_Vel_field == nullptr) {
    m_Vel_field = new Grid(
        continuum_state.physics_grid.min, continuum_state.physics_grid.max,
        continuum_state.physics_grid.cell_count.cast<int>(), 2);
  }
  if (m_Bebar_field == nullptr) {
    m_Bebar_field = new Grid(
        continuum_state.physics_grid.min, continuum_state.physics_grid.max,
        continuum_state.physics_grid.cell_count.cast<int>(), 4);
  }

  m_J_field->rasterizeMass(continuum_state.material_points,
                           numPointsBeforeSampling,
                           continuum_state.basis_functions);
  m_J_field->rasterizeMassWeightedValuesVector(
      continuum_state.material_points, numPointsBeforeSampling,
      continuum_state.basis_functions, continuum_state.material_points.J);
  VectorXs defaultJ;
  defaultJ.setZero(1);
  defaultJ(0) = 1.0;
  m_J_field->normalizationByMass(defaultJ);

  m_Vel_field->rasterizeMass(continuum_state.material_points,
                             numPointsBeforeSampling,
                             continuum_state.basis_functions);
  m_Vel_field->rasterizeMassWeightedValuesMatrix(
      continuum_state.material_points, numPointsBeforeSampling,
      continuum_state.basis_functions, continuum_state.material_points.v);
  VectorXs defaultV;
  defaultV.setZero(2);
  defaultV(0) = 0.0;
  defaultV(1) = 0.0;
  m_Vel_field->normalizationByMass(defaultV);

  // std::cout << "continuum_state.material_points: " <<
  // continuum_state.material_points.npoints << std::endl; std::cout <<
  // "checking vel_field... " << m_Vel_field->vals().norm() << std::endl;

  m_Bebar_field->rasterizeMass(continuum_state.material_points,
                               numPointsBeforeSampling,
                               continuum_state.basis_functions);
  m_Bebar_field->rasterizeMassWeightedValues(
      continuum_state.material_points, numPointsBeforeSampling,
      continuum_state.basis_functions, continuum_state.material_points.be_bar);
  VectorXs default_be_bar;
  default_be_bar.setZero(4);
  default_be_bar(0) = 1.0;
  default_be_bar(1) = 0.0;
  default_be_bar(2) = 0.0;
  default_be_bar(3) = 1.0;
  m_Bebar_field->normalizationByMass(default_be_bar);
}

void ResamplingTools::reassignJVelAndBebarForInternalMaterialPoint(
    SimulationState &continuum_state, const int pnt_idx,
    bool newly_inserted_continuum_stress_free) {
  const Vector2s x = continuum_state.material_points.q.col(pnt_idx);

  VectorXs _J = m_J_field->interpolateValue(continuum_state.basis_functions,
                                            continuum_state.material_points, x);
  scalar J = _J(0);

  VectorXs _be_bar = m_Bebar_field->interpolateValue(
      continuum_state.basis_functions, continuum_state.material_points, x);
  Eigen::Map<Matrix22sc> be_bar(_be_bar.data());
  scalar det_be_bar = be_bar.determinant();
  be_bar /= sqrt(det_be_bar);

  VectorXs _vel = m_Vel_field->interpolateValue(
      continuum_state.basis_functions, continuum_state.material_points, x);
  Vector2s vel = _vel.segment<2>(0);

  continuum_state.material_points.v.col(pnt_idx) = vel;

  if (newly_inserted_continuum_stress_free) {
    continuum_state.material_points.J(pnt_idx) = 1.0;
    continuum_state.material_points.be_bar[pnt_idx] = Matrix22sc::Identity();
  } else {
    continuum_state.material_points.J(pnt_idx) = J;
    continuum_state.material_points.be_bar[pnt_idx] = be_bar;
  }
}

// requirement:
// 1. the min_coord, res, and cell_count of the uniform grids should be
// identical to those of the mpm grid.
// 2. the mpm points and grains are already registered in the uniform grids.
// void ResamplingTools::determineMPMQuantitiesForInnerContinuumZones(
// SimulationState& continuum_state, const std::vector<Vector3u>& zone_indices,
// const int numPointsBeforeSampling, const int numPointsAfterSampling, const
// CUniformGrid& ugrd )
void ResamplingTools::determineMPMQuantitiesForInnerContinuumZones(
    SimulationState &continuum_state, const int numPointsBeforeSampling,
    const int numPointsAfterSampling, const CUniformGridH2D &ugrd,
    bool newly_inserted_continuum_stress_free) {
  computeFieldsForInterpolation(continuum_state, numPointsBeforeSampling);

  for (unsigned i = numPointsBeforeSampling;
       i < unsigned(numPointsAfterSampling); i++) {
    gatherMassAndVolume(continuum_state, ugrd, i);
  }

  for (unsigned pnt_idx = numPointsBeforeSampling;
       pnt_idx < unsigned(numPointsAfterSampling); ++pnt_idx) {
    reassignJVelAndBebarForInternalMaterialPoint(
        continuum_state, pnt_idx, newly_inserted_continuum_stress_free);
  }
}

void ResamplingTools::
    reassignMassAndVolumesForHybridAndNewlyGeneratedContinuumZone(
        SimulationState &continuum_state,
        const std::vector<Vector2u> &zone_indices,
        const CUniformGridH2D &ugrd_continuum, const scalar cell_volume) {
  const scalar mass_per_cell = continuum_state.material_density * cell_volume;
  for (int i = 0; i < int(zone_indices.size()); i++) {
    int nPts;
    const int *IDs;
    ugrd_continuum.getIDs(zone_indices[i], &nPts, &IDs);

    if (nPts < 0) {
      std::cout << "WARNING: a hybrid cell has no material point" << std::endl;
      continue;
    }

    for (int k = 0; k < nPts; k++) {
      continuum_state.material_points.m0[IDs[k]] = mass_per_cell / nPts;
      continuum_state.material_points.volume0[IDs[k]] = cell_volume / nPts;
    }
  }
}

void ResamplingTools::reassignJVelAndBebarForHybridMaterialPoint(
    SimulationState &continuum_state, const int pnt_idx, const scalar sigma,
    const scalar max_dist, const int numPointsBeforeSampling,
    const int numPointsAfterSampling, const RigidBody2DState &discrete_state,
    const int numGrainsBeforeSampling, const CUniformGridH2D &ugrd_continuum,
    const CUniformGridH2D &ugrd_discrete,
    bool newly_inserted_continuum_stress_free, bool homogenize_stress,
    PenaltyImpactFrictionMap *pifmap, bool homogenize_velocity,
    const scalar kappa, const scalar mu) {
  const Vector2s x = continuum_state.material_points.q.col(pnt_idx);
  const Vector2s min_c = x.array() - max_dist;
  const Vector2s max_c = x.array() + max_dist;

  Vector2s v_sum = Vector2s(0.0, 0.0);
  int nNeighbor = 0;
  scalar J_sum = 0.0;
  Matrix22sc be_bar_sum = Matrix22sc::Zero();
  scalar weight_sum = 0.0;
  Vector2s q0_sum = Vector2s(0.0, 0.0);

  Vector2u min_grid_idx, max_grid_idx;
  ugrd_continuum.getGridIDRange(min_c, max_c, min_grid_idx, max_grid_idx);
  for (unsigned b = min_grid_idx.y(); b <= max_grid_idx.y(); b++) {
    for (unsigned a = min_grid_idx.x(); a <= max_grid_idx.x(); a++) {
      Vector2u cell_idx = Vector2u{a, b};
      int nData;
      const int *IDs;
      ugrd_continuum.getIDs(cell_idx, &nData, &IDs);

      for (int k = 0; k < nData; k++) {
        if (IDs[k] >= numPointsBeforeSampling)
          continue;

        Vector2s p = continuum_state.material_points.q.col(IDs[k]);
        scalar dist = (p - x).norm();
        if (dist <= max_dist) {
          scalar weight = exp(-(dist * dist) / (2.0 * sigma * sigma));
          Vector2s v = continuum_state.material_points.v.col(IDs[k]);
          Vector2s q0 = continuum_state.material_points.q0.col(IDs[k]);
          scalar J = continuum_state.material_points.J(IDs[k]);
          v_sum += v * weight;
          q0_sum += q0 * weight;
          J_sum += J * weight;
          be_bar_sum += continuum_state.material_points.be_bar[IDs[k]] * weight;
          weight_sum += weight;
          nNeighbor++;
        }
      }
    }
  }

  Matrix22sc stress = Matrix22sc::Zero();

  if (homogenize_stress) {
    const scalar hl = continuum_state.physics_grid.cell_width * 0.25;
    // interpolate homogenized stress
    const std::pair<Array2u, Array2u> stencil{
        continuum_state.basis_functions->computeStencil(
            x, hl, continuum_state.physics_grid)};

    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(continuum_state.physics_grid.nodeIndicesValid(x_idx, y_idx));
        const scalar weight{continuum_state.basis_functions->weight(
            x, hl, {x_idx, y_idx}, continuum_state.physics_grid)};
        assert(weight >= 0.0);

        const unsigned int grid_index{
            continuum_state.physics_grid.flatNodeIndex(x_idx, y_idx)};

        // WARNING: This will only work for linear basis functions. cell_area
        // should account for the stencil size.
        stress += weight * m_homogenized_stress[grid_index];
      }
    }
  }

  Vector2s velocity = Vector2s::Zero();

  if (homogenize_velocity) {
    const scalar hl = continuum_state.physics_grid.cell_width * 0.25;
    // interpolate homogenized stress
    const std::pair<Array2u, Array2u> stencil{
        continuum_state.basis_functions->computeStencil(
            x, hl, continuum_state.physics_grid)};

    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(continuum_state.physics_grid.nodeIndicesValid(x_idx, y_idx));
        const scalar weight{continuum_state.basis_functions->weight(
            x, hl, {x_idx, y_idx}, continuum_state.physics_grid)};
        assert(weight >= 0.0);

        const unsigned int grid_index{
            continuum_state.physics_grid.flatNodeIndex(x_idx, y_idx)};

        // WARNING: This will only work for linear basis functions. cell_area
        // should account for the stencil size.
        velocity += weight * m_homogenized_velocity[grid_index];
      }
    }
  }

  scalar J;
  Matrix22sc be_bar;
  Vector2s q0;

  if (nNeighbor == 0) {
    J = 1.0;
    be_bar = Matrix22sc::Identity();
  } else {
    J = J_sum / weight_sum;
    be_bar = be_bar_sum / weight_sum;
    scalar det_be_bar = be_bar.determinant();
    be_bar = be_bar / sqrt(det_be_bar);
  }

  ugrd_discrete.getGridIDRange(min_c, max_c, min_grid_idx, max_grid_idx);

  for (unsigned b = min_grid_idx.y(); b <= max_grid_idx.y(); b++) {
    for (unsigned a = min_grid_idx.x(); a <= max_grid_idx.x(); a++) {
      Vector2u cell_idx = Vector2u{a, b};
      int nData;
      const int *IDs;
      ugrd_discrete.getIDs(cell_idx, &nData, &IDs);

      for (int k = 0; k < nData; k++) {
        if (IDs[k] >= numGrainsBeforeSampling)
          continue;

        Vector2s p = discrete_state.q().segment<2>(3 * IDs[k]);
        scalar dist = (p - x).norm();
        if (dist <= max_dist) {
          scalar weight = exp(-(dist * dist) / (2.0 * sigma * sigma));
          Vector2s v = discrete_state.v().segment<2>(3 * IDs[k]);
          Vector2s q0l = discrete_state.q0().segment<2>(3 * IDs[k]);
          v_sum += v * weight;
          q0_sum += q0l * weight;
          weight_sum += weight;
        }
      }
    }
  }

  Vector2s v;
  if (weight_sum == 0.0) {
    v = Vector2s::Zero();
    std::cout << "WARNING: material point " << pnt_idx
              << " has no material point or grain in its neighbor. This might "
                 "be a bug. "
              << std::endl;
    // potentially wrong q0
    q0 = x;
  } else {
    v = v_sum / weight_sum;
    q0 = q0_sum / weight_sum;
  }

  continuum_state.material_points.q0.col(pnt_idx) = q0;

  if (homogenize_velocity) {
    continuum_state.material_points.v.col(pnt_idx) = velocity;
  } else {
    continuum_state.material_points.v.col(pnt_idx) = v;
  }

  if (homogenize_stress) {
    ConstitutiveModel::computeStrainFromCauchyStress(stress, kappa, mu, J,
                                                     be_bar);
    // std::cout << "[" << stress(0, 0) << ", " << stress(0, 1) << ", " <<
    // stress(1, 0) << ", " << stress(1, 1) << "]=>(" << J << ", [" << be_bar(0,
    // 0) << ", " << be_bar(0, 1) << ", " << be_bar(1, 0) << ", " << be_bar(1,
    // 1)
    // << "]), ";
    continuum_state.material_points.J(pnt_idx) = J;
    continuum_state.material_points.be_bar[pnt_idx] = be_bar;
  } else {
    if (newly_inserted_continuum_stress_free)
      continuum_state.material_points.J(pnt_idx) = 1.0;
    else
      continuum_state.material_points.J(pnt_idx) = J;
    continuum_state.material_points.be_bar[pnt_idx] = be_bar;
  }
}

// requirement:
// 1. the min_coord, res, and cell_count of the uniform grids should be
// identical to those of the mpm grid.
// 2. the mpm points and grains are already registered (as points) in the
// uniform grids.
void ResamplingTools::determineMPMQuantitiesForHybridZones(
    SimulationState &continuum_state, const std::vector<Vector2u> &zone_indices,
    const int numPointsBeforeSampling, const int numPointsAfterSampling,
    const RigidBody2DState &discrete_state, const int numGrainsBeforeSampling,
    const CUniformGridH2D &ugrd_continuum, const CUniformGridH2D &ugrd_discrete,
    const scalar &cell_volume, bool newly_inserted_continuum_stress_free,
    bool homogenize_stress, PenaltyImpactFrictionMap *pifmap,
    bool homogenize_velocity, const scalar kappa, const scalar mu) {
  reassignMassAndVolumesForHybridAndNewlyGeneratedContinuumZone(
      continuum_state, zone_indices, ugrd_continuum, cell_volume);

  const scalar sigma = 0.5 * continuum_state.initial_particle_size;
  const scalar max_dist = 3.0 * continuum_state.initial_particle_size;

  // reassign J, v, be_bar
  // assign properties
  for (unsigned pnt_idx = numPointsBeforeSampling;
       pnt_idx < unsigned(numPointsAfterSampling); ++pnt_idx) {
    reassignJVelAndBebarForHybridMaterialPoint(
        continuum_state, pnt_idx, sigma, max_dist, numPointsBeforeSampling,
        numPointsAfterSampling, discrete_state, numGrainsBeforeSampling,
        ugrd_continuum, ugrd_discrete, newly_inserted_continuum_stress_free,
        homogenize_stress, pifmap, homogenize_velocity, kappa, mu);
  }
}

// requirement:
// 1. the min_coord, res, and cell_count of the uniform grids should be
// identical to those of the mpm grid.
// 2. the mpm points and grains are already registered (as points) in the
// uniform grids.
void ResamplingTools::determineDEMQuantitiesForHybridZones(
    RigidBody2DState &discrete_state, const int numGrainsBeforeSampling,
    const int numGrainsAfterSampling, const SimulationState &continuum_state,
    const int numPointsBeforeSampling, const CUniformGridH2D &ugrd_continuum,
    const CUniformGridH2D &ugrd_discrete) {
  const scalar sigma = m_dem_r_mean;
  const scalar max_dist = 6.0 * m_dem_r_mean;
  // const int numBodies = discrete_state.nbodies();
  // reassign v, rho
  for (unsigned bdy_idx = numGrainsBeforeSampling;
       bdy_idx < unsigned(numGrainsAfterSampling); ++bdy_idx) {
    const Vector2s x = discrete_state.q().segment<2>(3 * bdy_idx);
    const Vector2s min_c = x.array() - max_dist;
    const Vector2s max_c = x.array() + max_dist;

    Vector3s v_sum = Vector3s(0.0, 0.0, 0.0);
    int nNeighbor = 0;
    Vector3s q0_sum = Vector3s(0.0, 0.0, 0.0);
    scalar weight_sum = 0.0;

    Vector2u min_grid_idx, max_grid_idx;
    ugrd_discrete.getGridIDRange(min_c, max_c, min_grid_idx, max_grid_idx);
    for (unsigned b = min_grid_idx.y(); b <= max_grid_idx.y(); b++) {
      for (unsigned a = min_grid_idx.x(); a <= max_grid_idx.x(); a++) {
        Vector2u cell_idx = Vector2u{a, b};
        int nData;
        const int *IDs;
        ugrd_discrete.getIDs(cell_idx, &nData, &IDs);

        for (int k = 0; k < nData; k++) {
          if (IDs[k] >= numGrainsBeforeSampling)
            continue;

          Vector2s p = discrete_state.q().segment<2>(3 * IDs[k]);
          scalar dist = (p - x).norm();
          if (dist <= max_dist) {
            scalar weight = exp(-(dist * dist) / (2.0 * sigma * sigma));
            Vector3s v = discrete_state.v().segment<3>(3 * IDs[k]);
            Vector3s q0 = discrete_state.q0().segment<3>(3 * IDs[k]);
            v_sum += v * weight;
            q0_sum += q0 * weight;
            weight_sum += weight;
            nNeighbor++;
          }
        }
      }
    }

    ugrd_continuum.getGridIDRange(min_c, max_c, min_grid_idx, max_grid_idx);
    for (unsigned b = min_grid_idx.y(); b <= max_grid_idx.y(); b++) {
      for (unsigned a = min_grid_idx.x(); a <= max_grid_idx.x(); a++) {
        Vector2u cell_idx = Vector2u{a, b};
        int nData;
        const int *IDs;
        ugrd_continuum.getIDs(cell_idx, &nData, &IDs);

        for (int k = 0; k < nData; k++) {
          if (IDs[k] >= numPointsBeforeSampling)
            continue;

          Vector2s p = continuum_state.material_points.q.col(IDs[k]);
          scalar dist = (p - x).norm();
          if (dist <= max_dist) {
            scalar weight = exp(-(dist * dist) / (2.0 * sigma * sigma));
            Vector2s v = continuum_state.material_points.v.col(IDs[k]);
            Vector2s q0 = continuum_state.material_points.q0.col(IDs[k]);
            v_sum += Vector3s(v.x(), v.y(), 0.0) * weight;
            q0_sum += Vector3s(q0.x(), q0.y(), 0.0) * weight;
            weight_sum += weight;
          }
        }
      }
    }

    if (weight_sum == 0.0) {
      discrete_state.v().segment<3>(3 * bdy_idx) = Vector3s::Zero();
      std::cout << "WARNING: grain " << bdy_idx
                << " has no material point or grain in its neighbor. This "
                   "might be a bug. "
                << std::endl;
      // potentially wrong q0
      discrete_state.q0().segment<3>(3 * bdy_idx) = Vector3s(x.x(), x.y(), 0.0);
    } else {
      discrete_state.v().segment<3>(3 * bdy_idx) = v_sum / weight_sum;
      discrete_state.q0().segment<3>(3 * bdy_idx) = q0_sum / weight_sum;
    }
  }
}

void ResamplingTools::clearHomogenizedStress(
    const SimulationState &continuum_state) {
  m_homogenized_stress.clear();
  for (int i = 0; i < int(continuum_state.physics_grid.numGridPoints()); i++) {
    m_homogenized_stress.push_back(Matrix22sc::Zero());
  }
}

void ResamplingTools::homogenizeStress(const RigidBody2DState &discrete_state,
                                       const SimulationState &continuum_state,
                                       PenaltyImpactFrictionMap *pifmap) {
  // create a uniform grid, where each cell of the uniform grid is centered at
  // each mpm grid *corner*.
  scalar cell_width = continuum_state.physics_grid.cell_width;

  // std::cout << "Homogenizing stress..." << std::endl;
  // std::cout << "#grid points: " <<
  // continuum_state.physics_grid.numGridPoints() << std::endl;

  struct Collision {
    Vector2s p;
    Vector2s f;
    Vector2s r;
  };

  std::vector<Collision> collisions;

  CollisionCache &cache = pifmap->getCollisionCache();

  // Circle vs circle collisions
  for (int bdy_idx = 0; bdy_idx < int(discrete_state.numBodies()); bdy_idx++) {
    for (std::pair<int, CircleCircleCollision> &circle_vs_circle :
         cache.circle_circle_collisions.collisions()[bdy_idx]) {
      int bdy_idx1 = circle_vs_circle.first;

      Collision collision;
      collision.p = circle_vs_circle.second.contact_point;
      collision.r = discrete_state.q().segment<2>(bdy_idx * 3) -
                    discrete_state.q().segment<2>(bdy_idx1 * 3);
      collision.f = circle_vs_circle.second.normal_force +
                    circle_vs_circle.second.friction_force;

      collisions.push_back(collision);
    }
  }

  const scalar cell_area = cell_width * cell_width;

  const scalar hl = cell_width * 0.25;

  for (int i = 0; i < int(collisions.size()); i++) {
    Matrix22sc stressAtContactPoint =
        -0.5 * (collisions[i].f * collisions[i].r.transpose() +
                collisions[i].r * collisions[i].f.transpose());

    const std::pair<Array2u, Array2u> stencil{
        continuum_state.basis_functions->computeStencil(
            collisions[i].p, hl, continuum_state.physics_grid)};
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(continuum_state.physics_grid.nodeIndicesValid(x_idx, y_idx));
        const scalar weight{continuum_state.basis_functions->weight(
            collisions[i].p, hl, {x_idx, y_idx}, continuum_state.physics_grid)};
        assert(weight >= 0.0);

        const unsigned int grid_index{
            continuum_state.physics_grid.flatNodeIndex(x_idx, y_idx)};

        // WARNING: This will only work for linear basis functions. cell_area
        // should account for the stencil size.
        m_homogenized_stress[grid_index] +=
            weight * stressAtContactPoint / cell_area;
      }
    }
  }
}

void ResamplingTools::finalizeHomogenizedStress(
    const SimulationState &continuum_state, const scalar &factor) {
  for (int i = 0; i < int(continuum_state.physics_grid.numGridPoints()); i++) {
    m_homogenized_stress[i] *= factor;
  }
}

void ResamplingTools::clearHomogenizedVelocity(
    const SimulationState &continuum_state) {
  m_homogenized_velocity.clear();
  for (int i = 0; i < int(continuum_state.physics_grid.numGridPoints()); i++) {
    m_homogenized_velocity.push_back(Vector2s::Zero());
  }
}

void ResamplingTools::homogenizeVelocity(
    const RigidBody2DState &discrete_state,
    const SimulationState &continuum_state) {
  m_homogenized_momentum.clear();
  m_homogenized_mass.clear();

  for (int i = 0; i < int(continuum_state.physics_grid.numGridPoints()); i++) {
    m_homogenized_momentum.push_back(Vector2s::Zero());
    m_homogenized_mass.push_back(0.0);
  }

  scalar cell_width = continuum_state.physics_grid.cell_width;
  const scalar hl = cell_width * 0.25;

  for (int i = 0; i < int(discrete_state.numBodies()); i++) {
    if (discrete_state.fixed(i))
      continue;

    const std::unique_ptr<RigidBody2DGeometry> &current_geo{
        discrete_state.bodyGeometry(i)};
    if (current_geo->type() != RigidBody2DGeometryType::CIRCLE)
      continue;

    const scalar mass = discrete_state.m(i);
    const Vector2s velocity = discrete_state.v().segment<2>(3 * i);
    const Vector2s pos = discrete_state.q().segment<2>(3 * i);

    const std::pair<Array2u, Array2u> stencil{
        continuum_state.basis_functions->computeStencil(
            pos, hl, continuum_state.physics_grid)};
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(continuum_state.physics_grid.nodeIndicesValid(x_idx, y_idx));
        const scalar weight{continuum_state.basis_functions->weight(
            pos, hl, {x_idx, y_idx}, continuum_state.physics_grid)};
        assert(weight >= 0.0);

        const unsigned int grid_index{
            continuum_state.physics_grid.flatNodeIndex(x_idx, y_idx)};

        m_homogenized_momentum[grid_index] += weight * mass * velocity;
        m_homogenized_mass[grid_index] += weight * mass;
      }
    }
  }

  for (int i = 0; i < int(continuum_state.physics_grid.numGridPoints()); i++) {
    if (m_homogenized_mass[i] > 0.0) {
      m_homogenized_velocity[i] +=
          m_homogenized_momentum[i] / m_homogenized_mass[i];
    }
  }
}

void ResamplingTools::finalizeHomogenizedVelocity(
    const SimulationState &continuum_state, const scalar &factor) {
  for (int i = 0; i < int(continuum_state.physics_grid.numGridPoints()); i++) {
    m_homogenized_velocity[i] *= factor;
  }
}

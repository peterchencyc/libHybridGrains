#include "Coupling.h"
#include "mpmgrains2d/MaterialPoints.h"
#include "mpmgrains2d/PhysicsGrid.h"
#include "rigidbody2d/RigidBody2DSim.h"
#include <cmath>

#ifdef OPENMP_ENABLED
#include <algorithm>
#include <omp.h>
#endif

#ifndef NDEBUG
static Vector2s computeTotalGridMomentum(const Matrix2Xsc rasterized_momentum) {
  Vector2s total_mom = Vector2s::Zero();
  for (int node_idx = 0; node_idx < rasterized_momentum.cols(); node_idx++) {
    total_mom += rasterized_momentum.col(node_idx);
  }
  return total_mom;
}
#endif

void HybridCoupling2D::rasterizeGrainMass(
    SimulationState &continuum_sim, RigidBody2DSim &discrete_sim,
    const bool allowCouplingPartiallyFilledCell, const bool use_pre_position) {
#ifndef OPENMP_ENABLED
  const std::unique_ptr<BasisFunctions> &in_SF = continuum_sim.basis_functions;
  PhysicsGrid &grid = continuum_sim.physics_grid;
  const MaterialPoints &points = continuum_sim.material_points;
  RigidBody2DState &dstate = discrete_sim.state();

  const scalar default_mpm_hl = points.hl(0);

  const unsigned nnodes{grid.numGridPoints()};
  const unsigned nbodies{dstate.nbodies()};

  VectorXi dem_coupled{VectorXi::Zero(nbodies)};
  dem_mass_grid.setZero(nnodes);
  grid.hybridized.setZero(nnodes);

#ifndef NDEBUG
  scalar total_mass = 0;
#endif

  // For each discrete body
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    if (discrete_sim.isKinematicallyScripted(bdy_idx)) {
      continue;
    }

    const Vector2s q_bdy =
        use_pre_position ? Vector2s(dem_pos_before.segment<2>(3 * bdy_idx))
                         : Vector2s(dstate.q().segment<2>(3 * bdy_idx));
    const std::pair<Array2u, Array2u> stencil{
        in_SF->computeStencil(q_bdy, default_mpm_hl, grid)};

    bool all_masses_positive = true;

    // For each DEM body, check if it's surrounded by MPM nodes with positive
    // masses
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        if (!grid.nodeIndicesValid(x_idx, y_idx)) {
          continue;
        }

        const int p = grid.flatNodeIndex(x_idx, y_idx);
        if (grid.rasterized_mass(p) <= 1.0e-6) {
          all_masses_positive = false;
        }
      }
    }

    const bool coupled_grain =
        all_masses_positive || allowCouplingPartiallyFilledCell;

    if (!coupled_grain) {
      continue;
    }

    if (coupled_grain && !all_masses_positive) {
      std::cout << "Grain #" << bdy_idx
                << " will be coupled to a partially filled cell" << std::endl;
    }

    dem_coupled(bdy_idx) = 1;

    const scalar pnt_mass{dstate.getTotalMass(bdy_idx)};

#ifndef NDEBUG
    total_mass += pnt_mass;
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        if (!grid.nodeIndicesValid(x_idx, y_idx)) {
          continue;
        }
        const scalar w{
            in_SF->weight(q_bdy, default_mpm_hl, {x_idx, y_idx}, grid)};
#ifndef NDEBUG
        weight_sum += w;
#endif
        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < dem_mass_grid.size());
        dem_mass_grid(flat_node_idx) += w * pnt_mass;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

  packed_to_unpacked_dem = VectorXi::LinSpaced(nbodies, 0, nbodies - 1);
  packed_to_unpacked_dem.conservativeResize(
      std::stable_partition(
          packed_to_unpacked_dem.data(),
          packed_to_unpacked_dem.data() + packed_to_unpacked_dem.size(),
          [&dem_coupled](int i) { return dem_coupled(i) > 0; }) -
      packed_to_unpacked_dem.data());

  packed_to_unpacked_mpm = VectorXi::LinSpaced(nnodes, 0, nnodes - 1);
  packed_to_unpacked_mpm.conservativeResize(
      std::stable_partition(packed_to_unpacked_mpm.data(),
                            packed_to_unpacked_mpm.data() +
                                packed_to_unpacked_mpm.size(),
                            [this](int i) { return dem_mass_grid(i) > 0.0; }) -
      packed_to_unpacked_mpm.data());

#ifndef NDEBUG
  {
    // Postcondition: total grid mass is equal to total point mass
    scalar val1 = fabs(total_mass - dem_mass_grid.sum());
    ;
    if (val1 > 1.0e-6) {
      std::exit(EXIT_FAILURE);
    }
  }
#endif
#else
  const std::unique_ptr<BasisFunctions> &in_SF = continuum_sim.basis_functions;
  PhysicsGrid &grid = continuum_sim.physics_grid;
  const MaterialPoints &points = continuum_sim.material_points;
  RigidBody2DState &dstate = discrete_sim.state();

  const scalar default_mpm_hl = points.hl(0);

  const unsigned nnodes{grid.numGridPoints()};
  const unsigned nbodies{dstate.nbodies()};

  VectorXi dem_coupled{VectorXi::Zero(nbodies)};
  dem_mass_grid.setZero(nnodes);
  grid.hybridized.setZero(nnodes);

  // Allocate per-thread space to rasterize masses into
  static std::vector<VectorXs> per_thread_masses;
  if (per_thread_masses.empty()) {
    per_thread_masses.resize(omp_get_max_threads());
  }
  assert(int(per_thread_masses.size()) == omp_get_max_threads());

  // Zero out the storage
#pragma omp parallel for
  for (std::vector<VectorXs>::size_type idx = 0; idx < per_thread_masses.size();
       idx++) {
    per_thread_masses[idx].setZero(dem_mass_grid.size());
  }

#ifndef NDEBUG
  scalar total_mass = 0;
#endif

  // For each discrete body
#pragma omp parallel for
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; bdy_idx++) {
    const int tid{omp_get_thread_num()};
    assert(tid >= 0);
    assert(tid < int(per_thread_masses.size()));

    if (discrete_sim.isKinematicallyScripted(bdy_idx)) {
      continue;
    }

    const Vector2s q_bdy =
        use_pre_position ? Vector2s(dem_pos_before.segment<2>(3 * bdy_idx))
                         : Vector2s(dstate.q().segment<2>(3 * bdy_idx));
    const std::pair<Array2u, Array2u> stencil{
        in_SF->computeStencil(q_bdy, default_mpm_hl, grid)};

    bool all_masses_positive = true;

    // For each DEM body, check if it's surrounded by MPM nodes with positive
    // masses
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        if (!grid.nodeIndicesValid(x_idx, y_idx)) {
          continue;
        }

        const int p = grid.flatNodeIndex(x_idx, y_idx);
        if (grid.rasterized_mass(p) <= 1.0e-6) {
          all_masses_positive = false;
        }
      }
    }

    const bool coupled_grain =
        all_masses_positive || allowCouplingPartiallyFilledCell;

    if (!coupled_grain) {
      continue;
    }

    if (coupled_grain && !all_masses_positive) {
      std::cout << "Grain #" << bdy_idx
                << " will be coupled to a partially filled cell" << std::endl;
    }

    dem_coupled(bdy_idx) = 1;

    const scalar pnt_mass{dstate.getTotalMass(bdy_idx)};

#ifndef NDEBUG
#pragma omp critical
    total_mass += pnt_mass;
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        if (!grid.nodeIndicesValid(x_idx, y_idx)) {
          continue;
        }
        const scalar w{
            in_SF->weight(q_bdy, default_mpm_hl, {x_idx, y_idx}, grid)};
#ifndef NDEBUG
        weight_sum += w;
#endif
        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < dem_mass_grid.size());
        per_thread_masses[tid](flat_node_idx) += w * pnt_mass;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }
#pragma omp parallel for
  for (unsigned node_idx = 0; node_idx < dem_mass_grid.size(); node_idx++) {
    for (std::vector<VectorXs>::size_type idx = 0;
         idx < per_thread_masses.size(); idx++) {
      dem_mass_grid(node_idx) += per_thread_masses[idx](node_idx);
    }
  }

  packed_to_unpacked_dem = VectorXi::LinSpaced(nbodies, 0, nbodies - 1);
  packed_to_unpacked_dem.conservativeResize(
      std::stable_partition(
          packed_to_unpacked_dem.data(),
          packed_to_unpacked_dem.data() + packed_to_unpacked_dem.size(),
          [&dem_coupled](int i) { return dem_coupled(i) > 0; }) -
      packed_to_unpacked_dem.data());

  packed_to_unpacked_mpm = VectorXi::LinSpaced(nnodes, 0, nnodes - 1);
  packed_to_unpacked_mpm.conservativeResize(
      std::stable_partition(packed_to_unpacked_mpm.data(),
                            packed_to_unpacked_mpm.data() +
                                packed_to_unpacked_mpm.size(),
                            [this](int i) { return dem_mass_grid(i) > 0.0; }) -
      packed_to_unpacked_mpm.data());

#ifndef NDEBUG
  {
    // Postcondition: total grid mass is equal to total point mass
    scalar val1 = fabs(total_mass - dem_mass_grid.sum());
    ;
    if (val1 > 1.0e-6) {
      std::exit(EXIT_FAILURE);
    }
  }
#endif
#endif
}

void HybridCoupling2D::rasterizeGrainMom(
    const SimulationState &continuum_sim, const RigidBody2DSim &discrete_sim,
    const bool allowCouplingPartiallyFilledCell, const bool use_pre_position) {
#ifndef OPENMP_ENABLED
  const std::unique_ptr<BasisFunctions> &in_SF = continuum_sim.basis_functions;
  const PhysicsGrid &grid = continuum_sim.physics_grid;
  const MaterialPoints &points = continuum_sim.material_points;
  const RigidBody2DState &dstate = discrete_sim.state();

  const scalar default_mpm_hl = points.hl(0);

  const unsigned nnodes{grid.numGridPoints()};

  dem_momentum_grid.setZero(2, nnodes);

  const int num_packed_dem{static_cast<int>(packed_to_unpacked_dem.size())};

#ifndef NDEBUG
  Vector2s total_mom = Vector2s::Zero();
#endif

  // For each discrete body
  for (int i = 0; i < num_packed_dem; i++) {
    const int grain_flat_idx{packed_to_unpacked_dem(i)};

    if (discrete_sim.isKinematicallyScripted(grain_flat_idx)) {
      continue;
    }

    const Vector2s q_bdy =
        use_pre_position
            ? Vector2s(dem_pos_before.segment<2>(3 * grain_flat_idx))
            : Vector2s(dstate.q().segment<2>(3 * grain_flat_idx));
    const std::pair<Array2u, Array2u> stencil{
        in_SF->computeStencil(q_bdy, default_mpm_hl, grid)};

    const Vector2s pnt_mom{dstate.getTotalMass(grain_flat_idx) *
                           dstate.v().segment<2>(3 * grain_flat_idx)};

#ifndef NDEBUG
    total_mom += pnt_mom;
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        if (!grid.nodeIndicesValid(x_idx, y_idx)) {
          continue;
        }
        const scalar w{
            in_SF->weight(q_bdy, default_mpm_hl, {x_idx, y_idx}, grid)};
#ifndef NDEBUG
        weight_sum += w;
#endif
        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < dem_mass_grid.size());
        dem_momentum_grid.col(flat_node_idx) += w * pnt_mom;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

#ifndef NDEBUG
  {

    // Postcondition: grid momentum is equal to point momentum
    const Vector2s p0{total_mom};
    const Vector2s p1{computeTotalGridMomentum(dem_momentum_grid)};
    scalar val2 = (p0 - p1).lpNorm<Eigen::Infinity>();

    std::cout << "momentum difference: " << val2 << std::endl;
    if (val2 > 1.0e-6) {
      std::exit(EXIT_FAILURE);
    }
  }
#endif
#else
  const std::unique_ptr<BasisFunctions> &in_SF = continuum_sim.basis_functions;
  const PhysicsGrid &grid = continuum_sim.physics_grid;
  const MaterialPoints &points = continuum_sim.material_points;
  const RigidBody2DState &dstate = discrete_sim.state();

  const scalar default_mpm_hl = points.hl(0);

  const unsigned nnodes{grid.numGridPoints()};

  dem_momentum_grid.setZero(2, nnodes);

  // Allocate per-thread space to rasterize momentum into
  static std::vector<Matrix2Xsc> per_thread_momentum;
  if (per_thread_momentum.empty()) {
    per_thread_momentum.resize(omp_get_max_threads());
  }
  assert(int(per_thread_momentum.size()) == omp_get_max_threads());

  // Zero out the storage
#pragma omp parallel for
  for (std::vector<Matrix2Xsc>::size_type idx = 0;
       idx < per_thread_momentum.size(); idx++) {
    per_thread_momentum[idx].setZero(dem_momentum_grid.rows(),
                                     dem_momentum_grid.cols());
  }

  const int num_packed_dem{static_cast<int>(packed_to_unpacked_dem.size())};

#ifndef NDEBUG
  Vector2s total_mom = Vector2s::Zero();
#endif

// For each discrete body
#pragma omp parallel for
  for (int i = 0; i < num_packed_dem; i++) {
    const int tid{omp_get_thread_num()};
    assert(tid >= 0);
    assert(tid < int(per_thread_momentum.size()));

    const int grain_flat_idx{packed_to_unpacked_dem(i)};

    if (discrete_sim.isKinematicallyScripted(grain_flat_idx)) {
      continue;
    }

    const Vector2s q_bdy =
        use_pre_position
            ? Vector2s(dem_pos_before.segment<2>(3 * grain_flat_idx))
            : Vector2s(dstate.q().segment<2>(3 * grain_flat_idx));
    const std::pair<Array2u, Array2u> stencil{
        in_SF->computeStencil(q_bdy, default_mpm_hl, grid)};

    const Vector2s pnt_mom{dstate.getTotalMass(grain_flat_idx) *
                           dstate.v().segment<2>(3 * grain_flat_idx)};

#ifndef NDEBUG
    total_mom += pnt_mom;
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        if (!grid.nodeIndicesValid(x_idx, y_idx)) {
          continue;
        }
        const scalar w{
            in_SF->weight(q_bdy, default_mpm_hl, {x_idx, y_idx}, grid)};
#ifndef NDEBUG
        weight_sum += w;
#endif
        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < dem_mass_grid.size());
        per_thread_momentum[tid].col(flat_node_idx) += w * pnt_mom;
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }

#pragma omp parallel for
  for (unsigned node_idx = 0; node_idx < dem_mass_grid.size(); node_idx++) {
    for (std::vector<Matrix2Xsc>::size_type idx = 0;
         idx < per_thread_momentum.size(); idx++) {
      dem_momentum_grid.col(node_idx) += per_thread_momentum[idx].col(node_idx);
    }
  }
#ifndef NDEBUG
  {

    // Postcondition: grid momentum is equal to point momentum
    const Vector2s p0{total_mom};
    const Vector2s p1{computeTotalGridMomentum(dem_momentum_grid)};
    scalar val2 = (p0 - p1).lpNorm<Eigen::Infinity>();

    std::cout << "momentum difference: " << val2 << std::endl;
    if (val2 > 1.0e-6) {
      std::exit(EXIT_FAILURE);
    }
  }
#endif
#endif
}

void HybridCoupling2D::transferToGrain(
    const std::unique_ptr<BasisFunctions> &in_SF, const PhysicsGrid &grid,
    const MaterialPoints &points, RigidBody2DState &dstate,
    RigidBody2DSim &discrete_sim, const bool use_pre_position) {
  constexpr scalar alpha_phase_1 = 1.0; // 1.0: flip; 0.0: pic

  const int num_packed_dem{static_cast<int>(packed_to_unpacked_dem.size())};

  // Assign an h value to the discrete grains
  const scalar default_mpm_hl{points.hl(0)};
#pragma omp parallel for
  for (int i = 0; i < num_packed_dem; i++) {
    const int grain_flat_idx{packed_to_unpacked_dem(i)};

    if (discrete_sim.isKinematicallyScripted(grain_flat_idx)) {
      continue;
    }

    const Vector2s q_bdy =
        use_pre_position
            ? Vector2s(dem_pos_before.segment<2>(3 * grain_flat_idx))
            : Vector2s(dstate.q().segment<2>(3 * grain_flat_idx));

    const std::pair<Array2u, Array2u> stencil{
        in_SF->computeStencil(q_bdy, default_mpm_hl, grid)};

    Vector2s vpic{Vector2s::Zero()};
    Vector2s at{Vector2s::Zero()};

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        if (!grid.nodeIndicesValid(x_idx, y_idx)) {
          continue;
        }
        const scalar w{
            in_SF->weight(q_bdy, default_mpm_hl, {x_idx, y_idx}, grid)};
#ifndef NDEBUG
        weight_sum += w;
#endif
        const unsigned flat_node_idx{grid.flatNodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < dem_mass_grid.size());

        vpic += w * v_grid_after.col(flat_node_idx);
        at += w * at_grid_after.col(flat_node_idx);
      }
    }

    assert(fabs(weight_sum - 1.0) <= 1.0e-6);

    const Vector2s vflip{dstate.v().segment<2>(3 * grain_flat_idx) + at};

    dstate.v().segment<2>(3 * grain_flat_idx) =
        vpic * (1.0 - alpha_phase_1) + vflip * alpha_phase_1;
  }
}

void HybridCoupling2D::copyPrePositionsFull(RigidBody2DSim &discrete_sim) {
  dem_pos_before = discrete_sim.state().q();
}

void HybridCoupling2D::updateVelocitiesOnGrid(RigidBody2DSim &discrete_sim,
                                              SimulationState &continuum_sim) {
  const unsigned nnodes{continuum_sim.physics_grid.numGridPoints()};
  v_grid_after.setZero(2, nnodes);
  at_grid_after.setZero(2, nnodes);

  const int num_packed_mpm{static_cast<int>(packed_to_unpacked_mpm.size())};
#pragma omp parallel for
  for (int i = 0; i < num_packed_mpm; i++) {
    const int node_flat_idx{packed_to_unpacked_mpm(i)};

    const scalar m_mpm =
        continuum_sim.physics_grid.rasterized_mass(node_flat_idx);
    const scalar m_dem = dem_mass_grid(node_flat_idx);

    // regular
    {
      const scalar mass_sum = m_mpm + m_dem;
      const Vector2s lambda =
          (m_dem * continuum_sim.physics_grid.rasterized_momentum_new.col(
                       node_flat_idx) -
           m_mpm * dem_momentum_grid.col(node_flat_idx)) /
          mass_sum;

      continuum_sim.physics_grid.rasterized_momentum_new.col(node_flat_idx) -=
          lambda;
      const Vector2s dem_momentum_grid_after =
          dem_momentum_grid.col(node_flat_idx) + lambda;

      v_grid_after.col(node_flat_idx) = dem_momentum_grid_after / m_dem;
      at_grid_after.col(node_flat_idx) = lambda / (/*dt **/ m_dem);
    }

#ifndef NDEBUG
    const scalar val =
        (continuum_sim.physics_grid.rasterized_momentum_new.col(node_flat_idx) /
             m_mpm -
         v_grid_after.col(node_flat_idx))
            .lpNorm<Eigen::Infinity>();
    if (val > 1.0e-6) {
      std::cerr << "mpm & dem velocity difference after coupling: " << val
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
#endif
  }
}

void HybridCoupling2D::updatePositions_node_node(RigidBody2DSim &discrete_sim,
                                                 scalar dt) {
  const int num_packed_dem{static_cast<int>(packed_to_unpacked_dem.size())};
#pragma omp parallel for
  for (int i = 0; i < num_packed_dem; i++) {
    const int grain_flat_idx{packed_to_unpacked_dem(i)};
    discrete_sim.state().q().segment<2>(grain_flat_idx * 3) =
        dem_pos_before.segment<2>(grain_flat_idx * 3) +
        dt * discrete_sim.state().v().segment<2>(grain_flat_idx * 3);
  }
}

const VectorXs &HybridCoupling2D::getDEMGridMasses() const {
  return dem_mass_grid;
}

const VectorXi &HybridCoupling2D::getPackedToUnpackedMPM() const {
  return packed_to_unpacked_mpm;
}
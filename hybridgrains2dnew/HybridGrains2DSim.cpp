#include "HybridGrains2DSim.h"

#include <iomanip>
#include <iostream>

#include "mpmgrains2d/ExplicitIntegrator.h"

#include "rigidbody2d/PythonScripting.h"

#include "scisim/ConstrainedMaps/ImpactFrictionMap.h"
#include "scisim/HDF5File.h"
#include "scisim/Utilities.h"

#include "DiscreteIntegrator.h"
#include "Grid.h"

#include "scisim/Timer/TimeUtils.h"

HybridGrains2DSim::HybridGrains2DSim() : m_initialized(false) {}

HybridGrains2DSim::HybridGrains2DSim(
    const HybridIntegratorState &hybrid_integrator_state,
    const RigidBody2DState &discrete_state,
    const SimulationState &continuum_state)
    : m_avoid_a_void_enabled(hybrid_integrator_state.poormanSettings().enabled),
      m_avoid_a_void_step_counter(0),
      m_avoid_a_void_freq(
          hybrid_integrator_state.poormanSettings().avoidAVoidFreq),
      m_phi_threshold(hybrid_integrator_state.poormanSettings().phi_threshold),
      m_rzone_level_set(
          hybrid_integrator_state.poormanSettings().rzone_level_set),
      m_phi_window_size(
          hybrid_integrator_state.poormanSettings().phi_window_size),
      m_phi_samples_per_cell_side(
          hybrid_integrator_state.poormanSettings().phi_samples_per_cell_side),
      m_newly_inserted_continuum_stress_free(
          hybrid_integrator_state.poormanSettings()
              .newly_converted_continuum_stress_free),
      m_homogenize_stress(
          hybrid_integrator_state.poormanSettings().homogenize_stress),
      m_grid_smoothing_homogenized_stress(
          hybrid_integrator_state.poormanSettings()
              .grid_smoothing_homogenized_stress),
      m_homogenize_velocity(
          hybrid_integrator_state.poormanSettings().homogenize_velocity),
      m_grid_smoothing_homogenized_velocity(
          hybrid_integrator_state.poormanSettings()
              .grid_smoothing_homogenized_velocity),
      m_initialized(true), m_initial_resampling_performed(false),
      m_integrator_state(hybrid_integrator_state), m_discrete_sim(),
      m_continuum_state(continuum_state) {

  // clear up continuum data so avoid a void can be used to initialize continuum
  // data
  m_continuum_state.material_points.clear();

  m_discrete_sim.getState() = discrete_state;
  m_discrete_sim.clearConstraintCache();

  if (m_avoid_a_void_enabled) {
    m_ZoneTools.init(
        continuum_state,
        hybrid_integrator_state.poormanSettings().level_set_cell_width,
        hybrid_integrator_state.poormanSettings()
            .allow_direct_transitions_between_discrete_and_continuum);
    m_ResamplingTools.init(hybrid_integrator_state.poormanSettings().dem_r_mean,
                           hybrid_integrator_state.poormanSettings().dem_r_std);
  }

  m_continuum_state.physics_grid.clearRasterizedData();
  m_kinematically_locked_mpm_particles.clear();
  m_kinematically_locked_grains.clear();
}

void HybridGrains2DSim::initializeAvoidAVoid() {
  if (m_avoid_a_void_enabled) {
    m_ResamplingTools.clearHomogenizedStress(m_continuum_state);
    m_ResamplingTools.clearHomogenizedVelocity(m_continuum_state);

    updateHybridState(
        true,
        m_integrator_state.discreteIntegrator().impactFrictionMapPointer());
    m_avoid_a_void_step_counter = 0;
    updateHybridState(
        false,
        m_integrator_state.discreteIntegrator().impactFrictionMapPointer());
  }
}

HybridIntegratorState &HybridGrains2DSim::integratorState() {
  return m_integrator_state;
}

const HybridIntegratorState &HybridGrains2DSim::integratorState() const {
  return m_integrator_state;
}

RigidBody2DSim &HybridGrains2DSim::discreteSim() { return m_discrete_sim; }

const RigidBody2DSim &HybridGrains2DSim::discreteSim() const {
  return m_discrete_sim;
}

const RigidBody2DState &HybridGrains2DSim::discreteState() const {
  return m_discrete_sim.state();
}

const SimulationState &HybridGrains2DSim::continuumState() const {
  return m_continuum_state;
}

const ZoneTools &HybridGrains2DSim::zoneTools() const { return m_ZoneTools; }

bool HybridGrains2DSim::initialized() const { return m_initialized; }

void HybridGrains2DSim::
    stepSystem_constraintSetUp_prediction_correction_advection_node_node() {

  TimingTools timing_tools;

  if (!m_initial_resampling_performed) {
    updateHybridState(
        true,
        m_integrator_state.discreteIntegrator().impactFrictionMapPointer());
    updateHybridState(
        false,
        m_integrator_state.discreteIntegrator().impactFrictionMapPointer());
    m_initial_resampling_performed = true;
  }

  constexpr bool use_pre_position = true;
  // std::cout << "Saving pre-step discrete positions FULL...NN" << std::endl;
  m_hybrid_coupling.copyPrePositionsFull(m_discrete_sim);

  timing_tools.start();
  std::cout << "Stepping MPM Phase 0...NN" << std::endl;
  m_integrator_state.continuumIntegrator().step_zero_phase(m_continuum_state);
  timing_tools.stop("");

  timing_tools.start();
  std::cout << "rasterizing grain mass...NN" << std::endl;
  m_hybrid_coupling.rasterizeGrainMass(m_continuum_state, m_discrete_sim, false,
                                       use_pre_position);
  timing_tools.stop("   ");

  timing_tools.start();
  std::cout << "Stepping discrete integrator...NN" << std::endl;
  m_integrator_state.discreteIntegrator().step(m_discrete_sim);
  timing_tools.stop("");

  std::cout << "N: " << m_discrete_sim.state().nbodies() << std::endl;

  timing_tools.start();
  std::cout << "Stepping MPM Phase 1...NN" << std::endl;
  m_integrator_state.continuumIntegrator().step_first_phase(m_continuum_state);
  timing_tools.stop("");

  timing_tools.start();
  std::cout << "rasterizing grain momentum...NN" << std::endl;
  m_hybrid_coupling.rasterizeGrainMom(m_continuum_state, m_discrete_sim, false,
                                      use_pre_position);
  timing_tools.stop("   ");

  timing_tools.start();
  std::cout << "Updating velocities after coupling...NN" << std::endl;
  m_hybrid_coupling.updateVelocitiesOnGrid(m_discrete_sim, m_continuum_state);
  timing_tools.stop("   ");

  // !!!! TODO: This code is no longer valid. Just save out the weights now.
  // Save out the coupling weights for the continuum grid
  {
    /*
    const VectorXi& packed_to_unpacked_mpm =
    m_hybrid_coupling.getPackedToUnpackedMPM();
    m_continuum_state.physics_grid.homog_factor.setConstant(1.0);
    for( int idx = 0; idx < packed_to_unpacked_mpm.size(); idx++ )
    {
      m_continuum_state.physics_grid.homog_factor[packed_to_unpacked_mpm[idx]] =
    0.5;
    }
    */

    const VectorXi &packed_to_unpacked_mpm =
        m_hybrid_coupling.getPackedToUnpackedMPM();
    const int num_nodes = m_continuum_state.physics_grid.numGridPoints();
    for (int node_idx = 0; node_idx < num_nodes; node_idx++) {
      if (m_continuum_state.physics_grid.rasterized_mass(node_idx) > 0.0) {
        m_continuum_state.physics_grid.homog_factor[node_idx] = 1.0;
      } else {
        m_continuum_state.physics_grid.homog_factor[node_idx] = 0.0;
      }
    }

    for (int idx = 0; idx < packed_to_unpacked_mpm.size(); idx++) {
      m_continuum_state.physics_grid.homog_factor[packed_to_unpacked_mpm[idx]] =
          0.5;
    }
  }
  timing_tools.start();
  std::cout << "transfering back to grain...NN" << std::endl;
  m_hybrid_coupling.transferToGrain(
      m_continuum_state.basis_functions, m_continuum_state.physics_grid,
      m_continuum_state.material_points, m_discrete_sim.state(), m_discrete_sim,
      use_pre_position);
  timing_tools.stop("   ");

  timing_tools.start();
  std::cout << "Stepping continuum integrator (phase 2)...NN" << std::endl;
  m_integrator_state.continuumIntegrator().step_second_phase(m_continuum_state);
  timing_tools.stop("");

  timing_tools.start();
  std::cout << "Discrete position update (overwrite the positions of the "
               "coupled grains)..."
            << std::endl;
  m_hybrid_coupling.updatePositions_node_node(
      m_discrete_sim,
      scalar(m_integrator_state.discreteIntegrator().timestep()));
  timing_tools.stop("   ");

  // enrichment
  if (m_avoid_a_void_enabled > 0) {
    updateHybridState(
        false,
        m_integrator_state.discreteIntegrator().impactFrictionMapPointer());
  }
}

void HybridGrains2DSim::stepSystem() {
  assert(m_integrator_state.discreteIntegrator().computeTime() ==
         m_integrator_state.continuumIntegrator().computeTime());

  if ((m_integrator_state.overallTimestep() !=
       m_integrator_state.continuumIntegrator().timestep()) ||
      (m_integrator_state.overallTimestep() !=
       m_integrator_state.discreteIntegrator().timestep())) {
    std::cerr << "Currently only support [discrete dt == continuum dt] in 2d "
                 "hybrid sim."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Cache the start of step material point and discrete body locations
  m_mpm_q_start_of_step = m_continuum_state.material_points.q;
  m_grain_q_start = m_discrete_sim.getState().q();

  switch (m_integrator_state.integratorStyle()) {
  case HybridIntegratorState::IntegratorStyle::
      PREDICTION_CORRECTION_ADVECTION_NODE_NODE: {
    stepSystem_constraintSetUp_prediction_correction_advection_node_node();
    break;
  }
  default: {
    std::cerr << "unsupported integrator style" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }

  m_integrator_state.overallIteration()++;
}

#ifdef USE_HDF5
void HybridGrains2DSim::writeBinaryState(HDF5File &output_file) const {
#ifdef USE_HDF5
  // m_discrete_sim.writeBinaryState( "discrete", output_file );
  m_discrete_sim.writeBinaryState(output_file);
  m_continuum_state.writeBinaryState("continuum", output_file);
  if (m_hybrid_coupling.getDEMGridMasses().size() != 0) {
    assert(m_continuum_state.physics_grid.numGridPoints() ==
           m_hybrid_coupling.getDEMGridMasses().size());
    output_file.writeMatrix("discrete", "grid_node_masses",
                            m_hybrid_coupling.getDEMGridMasses());
  } else {
    VectorXs zeroMasses =
        VectorXs::Zero(m_continuum_state.physics_grid.numGridPoints());
    output_file.writeMatrix("discrete", "grid_node_masses", zeroMasses);
  }
#else
  std::cerr
      << "Error, HybridGrains2DSim::writeBinaryState requires HDF5 support."
      << std::endl;
  std::exit(EXIT_FAILURE);
#endif
}
#endif

void HybridGrains2DSim::serialize(std::ostream &output_stream) const {
  Utilities::serializeBuiltInType(m_initialized, output_stream);
  m_integrator_state.serialize(output_stream);
  m_discrete_sim.serialize(output_stream);
  m_continuum_state.serialize(output_stream);
}

void HybridGrains2DSim::deserialize(std::istream &input_stream) {
  m_initialized = Utilities::deserialize<bool>(input_stream);
  m_integrator_state.deserialize(input_stream);
  m_discrete_sim.deserialize(input_stream);
  m_continuum_state.deserialize(input_stream);
}

void HybridGrains2DSim::updateHybridStates_poorman(bool isInitializationStep,
                                                   ImpactFrictionMap *ifmap) {
  TimingTools timing_tools;

  if (m_avoid_a_void_step_counter == 0) {
    std::cout << "entering 2d poorman enrichment." << std::endl;

    std::cout << "updating packing fraction..." << std::endl;
    timing_tools.start();
    m_ZoneTools.updatePackingFraction(m_discrete_sim.state(), m_continuum_state,
                                      m_phi_window_size,
                                      m_phi_samples_per_cell_side);
    timing_tools.stop("");
    std::cout << "updating distance field..." << std::endl;
    timing_tools.start();
    m_ZoneTools.updateDistanceField(m_phi_threshold);
    timing_tools.stop("");
    std::cout << "computing reconciliation zonesupdating distance field..."
              << std::endl;
    timing_tools.start();
    m_ZoneTools.computeRZones(
        m_integrator_state.poormanSettings().rzone_level_set,
        m_integrator_state.poormanSettings().rzone_half_thickness);
    timing_tools.stop("");

    std::cout << "#Rzones: " << m_ZoneTools.getRZones().size() << std::endl;

    const int n_dem = m_discrete_sim.state().numBodies();
    const int n_mpm = m_continuum_state.material_points.npoints;

    int nGrainsBeforeSampling = 0, nGrainsAfterSampling = 0;
    std::cout << "(re)sampling DEM @ reconciliation zones..." << std::endl;
    timing_tools.start();
    m_ResamplingTools.avoidAVoidDEMparticlesInsideZonesWithTentativeQuantities(
        m_discrete_sim.state(), m_ZoneTools.getRZones(),
        m_integrator_state.poormanSettings().rho_dem, ifmap,
        nGrainsBeforeSampling, nGrainsAfterSampling);

    timing_tools.stop("");

    int nPointsBeforeSampling = 0, nPointsAfterSampling = 0;
    int _AfterSampling = 0, _BeforeSampling = 0;
    std::cout << "(re)sampling MPM @ reconciliation zones..." << std::endl;
    timing_tools.start();
    if (m_ZoneTools.allowDirectTransitionsBetweenDiscreteAndContinuum()) {
      m_ResamplingTools
          .avoidAVoidMPMparticlesInsideZonesWithTentativeQuantities(
              m_continuum_state, m_ZoneTools.getRZones(), nPointsBeforeSampling,
              _AfterSampling);
      m_ResamplingTools
          .avoidAVoidMPMparticlesInsideZonesWithTentativeQuantities(
              m_continuum_state,
              m_ZoneTools.getDirectTransitionZonesFromDiscreteToContinuum(),
              _BeforeSampling, nPointsAfterSampling);
    } else {
      m_ResamplingTools
          .avoidAVoidMPMparticlesInsideZonesWithTentativeQuantities(
              m_continuum_state, m_ZoneTools.getRZones(), nPointsBeforeSampling,
              nPointsAfterSampling);
    }
    timing_tools.stop("");

    // std::cout << "S1 mpm vel norm: " <<
    // m_continuum_state.material_points.v.norm() << ", #points: " <<
    // m_continuum_state.material_points.npoints << std::endl;

    int nPointsBeforeSampling2 = 0, nPointsAfterSampling2 = 0;
    std::cout << "(re)sampling MPM @ inner continuum zones..." << std::endl;
    timing_tools.start();
    std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>>
        inner_continuum_zones;
    std::vector<Vector2u> inner_continuum_zone_indices;
    if (m_ZoneTools.allowDirectTransitionsBetweenDiscreteAndContinuum()) {
      m_ZoneTools.identifyInnerContinuumRegionExcludingDirectTransitionZones(
          inner_continuum_zones, inner_continuum_zone_indices);
    } else {
      m_ZoneTools.identifyInnerContinuumRegion(inner_continuum_zones,
                                               inner_continuum_zone_indices);
    }
    if (inner_continuum_zones.size() == 0) {
      nPointsBeforeSampling2 = nPointsAfterSampling;
      nPointsAfterSampling2 = nPointsAfterSampling;
    } else {
      m_ResamplingTools
          .avoidAVoidMPMparticlesInsideZonesWithTentativeQuantities(
              m_continuum_state, inner_continuum_zones, nPointsBeforeSampling2,
              nPointsAfterSampling2);
      // std::cout << "S2 mpm vel norm: " <<
      // m_continuum_state.material_points.v.norm() << ", #points: " <<
      // m_continuum_state.material_points.npoints << std::endl;
    }
    timing_tools.stop("");

    CUniformGridH2D *ugrd_continuum =
        m_ZoneTools.setupUniformGridElementsAsPoints(m_continuum_state);
    CUniformGridH2D *ugrd_discrete =
        m_ZoneTools.setupUniformGridElementsAsPoints(m_discrete_sim.state());

    std::cout
        << "determining quantities for resampled MPM @ reconciliation zones..."
        << std::endl;
    timing_tools.start();
    const scalar h = m_continuum_state.physics_grid.cell_width;
    const scalar cell_volume = h * h;

    if (m_homogenize_stress) {
      std::cout << "finalizing homogenized stress..." << std::endl;
      m_ResamplingTools.homogenizeStress(
          m_discrete_sim.state(), m_continuum_state,
          (PenaltyImpactFrictionMap *)m_integrator_state.discreteIntegrator()
              .impactFrictionMapPointer());
      if (m_grid_smoothing_homogenized_stress) {
        m_ResamplingTools.finalizeHomogenizedStress(m_continuum_state,
                                                    1.0 / m_avoid_a_void_freq);
      }
    }

    if (m_homogenize_velocity) {
      std::cout << "finalizing homogenized velocity..." << std::endl;
      m_ResamplingTools.homogenizeVelocity(m_discrete_sim.state(),
                                           m_continuum_state);
      if (m_grid_smoothing_homogenized_velocity) {
        m_ResamplingTools.finalizeHomogenizedVelocity(
            m_continuum_state, 1.0 / m_avoid_a_void_freq);
      }
    }

    if (m_ZoneTools.allowDirectTransitionsBetweenDiscreteAndContinuum()) {
      m_ResamplingTools.determineMPMQuantitiesForHybridZones(
          m_continuum_state, m_ZoneTools.getRZoneIndices(),
          nPointsBeforeSampling, _AfterSampling, m_discrete_sim.state(),
          nGrainsBeforeSampling, *ugrd_continuum, *ugrd_discrete, cell_volume,
          m_newly_inserted_continuum_stress_free, m_homogenize_stress,
          (PenaltyImpactFrictionMap *)m_integrator_state.discreteIntegrator()
              .impactFrictionMapPointer(),
          m_homogenize_velocity, m_continuum_state.bulk_modulus,
          m_continuum_state.shear_modulus);
      // std::cout << "S3 mpm vel norm: " <<
      // m_continuum_state.material_points.v.norm() << ", #points: " <<
      // m_continuum_state.material_points.npoints << std::endl;
      m_ResamplingTools.determineMPMQuantitiesForHybridZones(
          m_continuum_state,
          m_ZoneTools.getDirectTransitionZonesFromDiscreteToContinuumIndices(),
          _BeforeSampling, nPointsAfterSampling, m_discrete_sim.state(),
          nGrainsBeforeSampling, *ugrd_continuum, *ugrd_discrete, cell_volume,
          m_newly_inserted_continuum_stress_free, m_homogenize_stress,
          (PenaltyImpactFrictionMap *)m_integrator_state.discreteIntegrator()
              .impactFrictionMapPointer(),
          m_homogenize_velocity, m_continuum_state.bulk_modulus,
          m_continuum_state.shear_modulus);
      // std::cout << "S4 mpm vel norm: " <<
      // m_continuum_state.material_points.v.norm() << ", #points: " <<
      // m_continuum_state.material_points.npoints << std::endl;
    } else {
      m_ResamplingTools.determineMPMQuantitiesForHybridZones(
          m_continuum_state, m_ZoneTools.getRZoneIndices(),
          nPointsBeforeSampling, nPointsAfterSampling, m_discrete_sim.state(),
          nGrainsBeforeSampling, *ugrd_continuum, *ugrd_discrete, cell_volume,
          m_newly_inserted_continuum_stress_free, m_homogenize_stress,
          (PenaltyImpactFrictionMap *)m_integrator_state.discreteIntegrator()
              .impactFrictionMapPointer(),
          m_homogenize_velocity, m_continuum_state.bulk_modulus,
          m_continuum_state.shear_modulus);
      // std::cout << "S5 mpm vel norm: " <<
      // m_continuum_state.material_points.v.norm() << ", #points: " <<
      // m_continuum_state.material_points.npoints << std::endl;
    }
    timing_tools.stop("");

    std::cout
        << "determining quantities for resampled DEM @ reconciliation zones..."
        << std::endl;
    timing_tools.start();
    m_ResamplingTools.determineDEMQuantitiesForHybridZones(
        m_discrete_sim.state(), nGrainsBeforeSampling, nGrainsAfterSampling,
        m_continuum_state, nPointsBeforeSampling, *ugrd_continuum,
        *ugrd_discrete);
    timing_tools.stop("");

    std::cout
        << "determining quantities for resampled MPM @ inner continuum zones..."
        << std::endl;
    timing_tools.start();
    m_ResamplingTools.determineMPMQuantitiesForInnerContinuumZones(
        m_continuum_state, nPointsBeforeSampling2, nPointsAfterSampling2,
        *ugrd_continuum, m_newly_inserted_continuum_stress_free);
    // std::cout << "S6 mpm vel norm: " <<
    // m_continuum_state.material_points.v.norm() << ", #points: " <<
    // m_continuum_state.material_points.npoints << std::endl;
    timing_tools.stop("");

    delete ugrd_continuum;
    delete ugrd_discrete;

    const int n_dem_post = m_discrete_sim.state().numBodies();
    const int n_mpm_post = m_continuum_state.material_points.npoints;

    std::cout << "deleting particles..." << std::endl;
    timing_tools.start();
    m_ResamplingTools.deleteParticlesUsingZoneIndicators(
        m_ZoneTools, m_discrete_sim.state(), m_continuum_state, ifmap);
    timing_tools.stop("");

    const int n_dem_post_del = m_discrete_sim.state().numBodies();
    const int n_mpm_post_del = m_continuum_state.material_points.npoints;

    std::cout << "avoid_a_void_dem: " << n_dem << "->" << n_dem_post << "->"
              << n_dem_post_del << std::endl;
    std::cout << "avoid_a_void_mpm: " << n_mpm << "->" << n_mpm_post << "->"
              << n_mpm_post_del << std::endl;

    std::function<scalar(const VectorXs &)> dem_weight_func = std::bind(
        &ZoneTools::demWeight, std::ref(m_ZoneTools), std::placeholders::_1);
    std::function<scalar(const VectorXs &)> mpm_weight_func = std::bind(
        &ZoneTools::mpmWeight, std::ref(m_ZoneTools), std::placeholders::_1);

    std::cout << "updating discrete mass weight..." << std::endl;
    std::cout << "discrete mass before update: "
              << m_discrete_sim.state().totalSimulatedMass() << std::endl;
    m_discrete_sim.state().updateMassWeights(dem_weight_func);
    m_discrete_sim.state().updateCurrentMassUsingMassWeights();
    std::cout << "discrete mass after update: "
              << m_discrete_sim.state().totalSimulatedMass() << std::endl;

    std::cout << "updating continuum mass weight..." << std::endl;
    std::cout << "continuum mass before update: "
              << m_continuum_state.material_points.totalMass() << std::endl;
    m_continuum_state.material_points.updateWeights(mpm_weight_func);
    m_continuum_state.material_points.updateMassAndVolumeUsingWeights();
    std::cout << "continuum mass after update: "
              << m_continuum_state.material_points.totalMass() << std::endl;

    if (m_homogenize_stress) {
      std::cout << "clear homogenized stress..." << std::endl;
      m_ResamplingTools.clearHomogenizedStress(m_continuum_state);
    }
    if (m_homogenize_velocity) {
      std::cout << "clear homogenized velocity..." << std::endl;
      m_ResamplingTools.clearHomogenizedVelocity(m_continuum_state);
    }
  } else if (m_avoid_a_void_step_counter >= m_avoid_a_void_freq) {
    m_avoid_a_void_step_counter = -1;
  }

  if (m_homogenize_stress && m_grid_smoothing_homogenized_stress) {
    std::cout << "accumulating homogenized stress..." << std::endl;
    m_ResamplingTools.homogenizeStress(
        m_discrete_sim.state(), m_continuum_state,
        (PenaltyImpactFrictionMap *)m_integrator_state.discreteIntegrator()
            .impactFrictionMapPointer());
  }
  if (m_homogenize_velocity && m_grid_smoothing_homogenized_velocity) {
    std::cout << "accumulating homogenized velocity..." << std::endl;
    m_ResamplingTools.homogenizeVelocity(m_discrete_sim.state(),
                                         m_continuum_state);
  }

  m_avoid_a_void_step_counter++;
}

void HybridGrains2DSim::updateHybridState(bool isInitializationStep,
                                          ImpactFrictionMap *ifmap) {
  if (m_avoid_a_void_enabled) {
    updateHybridStates_poorman(isInitializationStep, ifmap);
  }
}

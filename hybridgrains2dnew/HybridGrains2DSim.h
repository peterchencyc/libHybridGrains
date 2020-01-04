#ifndef HYBRID_GRAINS_2D_SIM_H
#define HYBRID_GRAINS_2D_SIM_H

#include "mpmgrains2d/SimulationState.h"

#include "rigidbody2d/RigidBody2DSim.h"

#include "Coupling.h"
#include "HybridIntegratorState.h"
#include "ResamplingTools.h"
#include "ZoneTools.h"

class ImpactFrictionMap;

class HybridGrains2DSim final {

public:
  HybridGrains2DSim();
  HybridGrains2DSim(const HybridIntegratorState &hybrid_integrator_state,
                    const RigidBody2DState &discrete_state,
                    const SimulationState &continuum_state);

  void initializeAvoidAVoid();

  HybridIntegratorState &integratorState();
  const HybridIntegratorState &integratorState() const;

  RigidBody2DSim &discreteSim();
  const RigidBody2DSim &discreteSim() const;
  const RigidBody2DState &discreteState() const;

  const SimulationState &continuumState() const;

  const ZoneTools &zoneTools() const;

  bool initialized() const;

  void stepSystem();

#ifdef USE_HDF5
  void writeBinaryState(HDF5File &output_file) const;
#endif

  void serialize(std::ostream &output_stream) const;
  void deserialize(std::istream &input_stream);

  // void executeTest();

private:
  bool m_avoid_a_void_enabled;
  int m_avoid_a_void_step_counter;
  int m_avoid_a_void_freq;
  double m_phi_threshold;
  double m_rzone_level_set;
  double m_phi_window_size;
  int m_phi_samples_per_cell_side;
  bool m_newly_inserted_continuum_stress_free;
  bool m_homogenize_stress;
  bool m_grid_smoothing_homogenized_stress;
  bool m_homogenize_velocity;
  bool m_grid_smoothing_homogenized_velocity;

  void updateHybridStates_poorman(bool isInitializationStep,
                                  ImpactFrictionMap *ifmap);
  void updateHybridState(bool isInitializationStep, ImpactFrictionMap *ifmap);

  void stepSystem_constraintSetUp_prediction_correction_advection_node_node();

  bool m_initialized;
  bool m_initial_resampling_performed;

  // Current state of the hybrid integration
  HybridIntegratorState m_integrator_state;

  // Discrete state
  RigidBody2DSim m_discrete_sim;

  // Continuum state
  SimulationState m_continuum_state;

  HybridCoupling2D m_hybrid_coupling;
  Eigen::VectorXd m_LMLambda;

  std::vector<unsigned> m_bodies_to_stabilize;

  Matrix2Xsc m_mpm_q_start_of_step;
  VectorXs m_grain_q_start;

  ZoneTools m_ZoneTools;
  ResamplingTools m_ResamplingTools;

  std::vector<int> m_kinematically_locked_mpm_particles;
  std::vector<int> m_kinematically_locked_grains;
};

#endif

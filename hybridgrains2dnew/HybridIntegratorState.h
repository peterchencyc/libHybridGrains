#ifndef HYBRID_INTEGRATOR_STATE_H
#define HYBRID_INTEGRATOR_STATE_H

#include "hybridgrains2dnew/DiscreteIntegrator.h"
#include "hybridgrains2dnew/MPMIntegrator.h"
#include "scisim/Math/Rational.h"

struct PoormanEnrichmentSettings final {
  PoormanEnrichmentSettings()
      : enabled(false), rzone_half_thickness(-1.0), dem_r_mean(-1.0),
        dem_r_std(-1.0), avoidAVoidFreq(-1), phi_threshold(-1.0),
        rzone_level_set(-1.0), phi_window_size(-1.0),
        level_set_cell_width(-1.0), phi_samples_per_cell_side(-1),
        target_dense_packing_fraction(-1.0), rho_dem(-1.0),
        allow_direct_transitions_between_discrete_and_continuum(false),
        newly_converted_continuum_stress_free(false), homogenize_stress(false),
        grid_smoothing_homogenized_stress(false), homogenize_velocity(false),
        grid_smoothing_homogenized_velocity(false) {}

  bool enabled;
  scalar rzone_half_thickness;
  scalar dem_r_mean;
  scalar dem_r_std;
  int avoidAVoidFreq;
  scalar phi_threshold;
  scalar rzone_level_set;
  scalar phi_window_size;
  scalar level_set_cell_width;
  int phi_samples_per_cell_side;
  scalar target_dense_packing_fraction;
  scalar rho_dem;
  bool allow_direct_transitions_between_discrete_and_continuum;
  bool newly_converted_continuum_stress_free;
  bool homogenize_stress;
  bool grid_smoothing_homogenized_stress;
  bool homogenize_velocity;
  bool grid_smoothing_homogenized_velocity;

  // bool LMCoupling;
  // bool LM_pos_correction;
};

class HybridIntegratorState final {

public:
  enum class IntegratorStyle {
    OLD_VERSION,
    OLD_VERSION_NODE_NODE,
    PREDICTION_CORRECTION_ADVECTION,
    PREDICTION_CORRECTION_ADVECTION_NODE_NODE,
    ITERATIVE
  };

  enum class HybridKinematicScript { NONE, SCRIPT_HYBRID_FRONTS };

  HybridIntegratorState();

  HybridIntegratorState(
      const Rational<std::intmax_t> &overall_dt,
      const Rational<std::intmax_t> &end_time,
      const DiscreteIntegrator &discrete_integrator,
      const MPMIntegrator &mpm_integrator,
      const IntegratorStyle integrator_style,
      const HybridKinematicScript kinematically_scripted_hybrid_fronts,
      PoormanEnrichmentSettings poorman_settings);

  scalar computeCurrentTime() const;

  void setEndTime(const Rational<std::intmax_t> &end_time);
  const Rational<std::intmax_t> &endTime() const;

  unsigned &overallIteration();

  const unsigned &overallIteration() const;

  const Rational<std::intmax_t> &overallTimestep() const;

  DiscreteIntegrator &discreteIntegrator();
  const DiscreteIntegrator &discreteIntegrator() const;

  MPMIntegrator &continuumIntegrator();

  void serialize(std::ostream &output_stream) const;

  void deserialize(std::istream &input_stream);

  IntegratorStyle integratorStyle() const;

  HybridKinematicScript kinematicallyScriptedHybridFronts() const;

  const PoormanEnrichmentSettings &poormanSettings() const;

private:
  unsigned m_overall_iteration;
  Rational<std::intmax_t> m_overall_dt;
  Rational<std::intmax_t> m_end_time;
  DiscreteIntegrator m_discrete_integrator;
  MPMIntegrator m_mpm_integrator;
  IntegratorStyle m_integrator_style;
  HybridKinematicScript m_kinematically_scripted_hybrid_fronts;
  PoormanEnrichmentSettings m_poorman_settings;
};

#endif

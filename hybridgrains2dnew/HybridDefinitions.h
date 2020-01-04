#ifndef HYBRID_DEFINITIONS_H
#define HYBRID_DEFINITIONS_H

#include "HybridIntegratorState.h"
#include "scisim/Math/MathDefines.h"
#include <bitset>
#include <iostream>

struct HybridIntegratorSettings final {
  std::string discrete_file_name;
  std::string continuum_file_name;
  Rational<std::intmax_t> overall_dt;
  std::string time_step_string;
  Rational<std::intmax_t> end_time;
  HybridIntegratorState::IntegratorStyle integrator_style;
  HybridIntegratorState::HybridKinematicScript
      kinematically_scripted_hybrid_fronts;

  // All bodies that initially intersect these rectangular regions are deleted
  std::vector<std::pair<Array2s, Array2s>> discrete_body_masks;

  PoormanEnrichmentSettings poorman_settings;
};

#endif

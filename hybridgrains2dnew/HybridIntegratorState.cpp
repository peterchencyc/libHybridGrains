#include "HybridIntegratorState.h"

HybridIntegratorState::HybridIntegratorState()
    : m_overall_iteration(0), m_overall_dt(0.0), m_end_time(0, 1),
      m_discrete_integrator(), m_mpm_integrator(),
      m_integrator_style(IntegratorStyle::PREDICTION_CORRECTION_ADVECTION),
      m_kinematically_scripted_hybrid_fronts(HybridKinematicScript::NONE),
      m_poorman_settings() {}

HybridIntegratorState::HybridIntegratorState(
    const Rational<std::intmax_t> &overall_dt,
    const Rational<std::intmax_t> &end_time,
    const DiscreteIntegrator &discrete_integrator,
    const MPMIntegrator &mpm_integrator, const IntegratorStyle integrator_style,
    const HybridKinematicScript kinematically_scripted_hybrid_fronts,
    PoormanEnrichmentSettings poorman_settings)
    : m_overall_iteration(0), m_overall_dt(overall_dt), m_end_time(end_time),
      m_discrete_integrator(discrete_integrator),
      m_mpm_integrator(mpm_integrator), m_integrator_style(integrator_style),
      m_kinematically_scripted_hybrid_fronts(
          kinematically_scripted_hybrid_fronts),
      m_poorman_settings(poorman_settings) {
  assert(m_overall_dt.positive());
  assert(m_end_time.positive());
  assert(m_discrete_integrator.computeTime() == m_mpm_integrator.computeTime());
  assert(m_discrete_integrator.computeTime() == 0.0);
}

scalar HybridIntegratorState::computeCurrentTime() const {
  return m_overall_iteration * scalar(m_overall_dt);
}

void HybridIntegratorState::setEndTime(
    const Rational<std::intmax_t> &end_time) {
  assert(end_time.positive());
  m_end_time = end_time;
}

const Rational<std::intmax_t> &HybridIntegratorState::endTime() const {
  return m_end_time;
}

unsigned &HybridIntegratorState::overallIteration() {
  return m_overall_iteration;
}

const unsigned &HybridIntegratorState::overallIteration() const {
  return m_overall_iteration;
}

const Rational<std::intmax_t> &HybridIntegratorState::overallTimestep() const {
  return m_overall_dt;
}

DiscreteIntegrator &HybridIntegratorState::discreteIntegrator() {
  return m_discrete_integrator;
}

const DiscreteIntegrator &HybridIntegratorState::discreteIntegrator() const {
  return m_discrete_integrator;
}

MPMIntegrator &HybridIntegratorState::continuumIntegrator() {
  return m_mpm_integrator;
}

void HybridIntegratorState::serialize(std::ostream &output_stream) const {
  std::cerr << "Update HybridIntegratorState::serialize" << std::endl;
  std::exit(EXIT_FAILURE);
  // Utilities::serializeBuiltInType( m_overall_iteration, output_stream );
  // RationalTools::serialize( m_overall_dt, output_stream );
  // RationalTools::serialize( m_end_time, output_stream );
  // m_discrete_integrator.serialize( output_stream );
  // m_mpm_integrator.serialize( output_stream );
}

void HybridIntegratorState::deserialize(std::istream &input_stream) {
  std::cerr << "Update HybridIntegratorState::deserialize" << std::endl;
  std::exit(EXIT_FAILURE);
  // m_overall_iteration = Utilities::deserialize<unsigned>( input_stream );
  // RationalTools::deserialize( m_overall_dt, input_stream );
  // RationalTools::deserialize( m_end_time, input_stream );
  // m_discrete_integrator.deserialize(input_stream);
  // m_mpm_integrator.deserialize(input_stream);
}

HybridIntegratorState::IntegratorStyle
HybridIntegratorState::integratorStyle() const {
  return m_integrator_style;
}

HybridIntegratorState::HybridKinematicScript
HybridIntegratorState::kinematicallyScriptedHybridFronts() const {
  return m_kinematically_scripted_hybrid_fronts;
}

const PoormanEnrichmentSettings &
HybridIntegratorState::poormanSettings() const {
  return m_poorman_settings;
}

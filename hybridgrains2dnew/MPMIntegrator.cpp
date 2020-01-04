#include "MPMIntegrator.h"

#include "mpmgrains2d/ExplicitIntegrator.h"
#include "mpmgrains2d/SimulationState.h"

//#include <iostream>

MPMIntegrator::MPMIntegrator() : m_iteration(0), m_dt(0, 1) {}

MPMIntegrator::MPMIntegrator(const Rational<std::intmax_t> &dt)
    : m_iteration(0), m_dt(dt) {
  assert(m_dt.positive());
}

const Rational<std::intmax_t> &MPMIntegrator::timestep() const { return m_dt; }

unsigned MPMIntegrator::iteration() const { return m_iteration; }

void MPMIntegrator::step(SimulationState &state) {
  ExplicitIntegrator::flow(scalar(m_dt), state);

  m_iteration++;
}

void MPMIntegrator::step_zero_phase(SimulationState &state) {
  ExplicitIntegrator::flow_zero_phase(state);
}

void MPMIntegrator::step_first_phase(SimulationState &state) {
  ExplicitIntegrator::flow_first_phase(scalar(m_dt), state);
}

void MPMIntegrator::step_second_phase(SimulationState &state) {
  ExplicitIntegrator::flow_second_phase(scalar(m_dt), state);

  m_iteration++;
}

scalar MPMIntegrator::computeTime() const {
  return scalar(std::intmax_t(m_iteration) * m_dt);
}

void MPMIntegrator::serialize(std::ostream &output_stream) const {
  Utilities::serializeBuiltInType(m_iteration, output_stream);
  RationalTools::serialize(m_dt, output_stream);
}

void MPMIntegrator::deserialize(std::istream &input_stream) {
  m_iteration = Utilities::deserialize<unsigned>(input_stream);
  RationalTools::deserialize(m_dt, input_stream);
}

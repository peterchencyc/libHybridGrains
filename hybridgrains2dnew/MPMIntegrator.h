#ifndef MPM_INTEGRATOR_H
#define MPM_INTEGRATOR_H

//#include <memory>

#include "mpmgrains2d/SimulationState.h"
#include "scisim/Math/Rational.h"

struct SimulationState;

class MPMIntegrator final {

public:
  MPMIntegrator();
  explicit MPMIntegrator(const Rational<std::intmax_t> &dt);

  const Rational<std::intmax_t> &timestep() const;

  unsigned iteration() const;

  void step(SimulationState &state);
  void step_zero_phase(SimulationState &state);
  void step_first_phase(SimulationState &state);
  void step_second_phase(SimulationState &state);

  scalar computeTime() const;

  void serialize(std::ostream &output_stream) const;
  void deserialize(std::istream &input_stream);

private:
  unsigned m_iteration;
  Rational<std::intmax_t> m_dt;
};

#endif

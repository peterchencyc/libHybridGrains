#ifndef EXPLICIT_INTEGRATOR_2D_H
#define EXPLICIT_INTEGRATOR_2D_H

#include "scisim/Math/MathDefines.h"

struct SimulationState;

namespace ExplicitIntegrator {
void flow(const scalar &dt, SimulationState &state);
void flow_zero_phase(SimulationState &state);
void flow_first_phase(const scalar &dt, SimulationState &state);
void flow_second_phase(const scalar &dt, SimulationState &state);
} // namespace ExplicitIntegrator

namespace MPMGrains2DSim {
void flow(const scalar &dt, SimulationState &state);
void flow_zero_phase(SimulationState &state);
void flow_first_phase(const scalar &dt, SimulationState &state);
void flow_second_phase(const scalar &dt, SimulationState &state);
} // namespace MPMGrains2DSim

#endif

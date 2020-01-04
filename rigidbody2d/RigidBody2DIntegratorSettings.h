#ifndef RIGID_BODY_2D_INTEGRATOR_SETTINGS_H
#define RIGID_BODY_2D_INTEGRATOR_SETTINGS_H

#include "scisim/ConstrainedMaps/FrictionSolver.h"
#include "scisim/ConstrainedMaps/ImpactFrictionMap.h"
#include "scisim/ConstrainedMaps/ImpactMaps/ImpactOperator.h"
#include "scisim/UnconstrainedMaps/UnconstrainedMap.h"

struct RigidBody2DIntegratorSettings final {
  Rational<std::intmax_t> dt;
  scalar end_time;

  std::unique_ptr<UnconstrainedMap> unconstrained_map;
  std::unique_ptr<ImpactOperator> impact_operator;
  std::unique_ptr<FrictionSolver> friction_solver;
  std::unique_ptr<ImpactFrictionMap> if_map;
  std::unique_ptr<ImpactMap> impact_map;

  scalar spatial_grid_scale;

  bool reduce_bandwidth;

  // TODO: Move these out of here
  scalar CoR;
  scalar mu;
};

#endif

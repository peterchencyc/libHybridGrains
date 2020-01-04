#ifndef DISCRETE_INTEGRATOR_H
#define DISCRETE_INTEGRATOR_H

#include <memory>

#include "scisim/ConstrainedMaps/FrictionSolver.h"
#include "scisim/ConstrainedMaps/ImpactFrictionMap.h"
#include "scisim/ConstrainedMaps/ImpactMaps/ImpactMap.h"
#include "scisim/ConstrainedMaps/ImpactMaps/ImpactOperator.h"
#include "scisim/Math/MathDefines.h"
#include "scisim/Math/Rational.h"
#include "scisim/UnconstrainedMaps/UnconstrainedMap.h"

#include "rigidbody2d/PythonScripting.h"

class RigidBody2DSim;
class PythonScripting;

class DiscreteIntegrator final {

public:
  DiscreteIntegrator();
  DiscreteIntegrator(
      const Rational<std::intmax_t> &dt,
      const std::unique_ptr<UnconstrainedMap> &unconstrained_map,
      const std::unique_ptr<ImpactOperator> &impact_operator,
      const std::unique_ptr<FrictionSolver> &friction_solver,
      const std::unique_ptr<ImpactFrictionMap> &impact_friction_map,
      const std::unique_ptr<ImpactMap> &impact_map, const scalar &cor,
      const scalar &mu, const bool &reduce_bandwidth);

  DiscreteIntegrator(const DiscreteIntegrator &other);
  DiscreteIntegrator(DiscreteIntegrator &&) = default;

  DiscreteIntegrator &operator=(DiscreteIntegrator other);
  DiscreteIntegrator &operator=(DiscreteIntegrator &&) = default;

  ~DiscreteIntegrator() = default;

  void setPythonCallback(const std::string &path,
                         const std::string &module_name);
  void pythonStartOfSim(RigidBody2DSim &sim);
  void pythonEndOfSim(RigidBody2DSim &sim);

  const Rational<std::intmax_t> &timestep() const;

  void step(RigidBody2DSim &sim);

  scalar computeTime() const;

  std::unique_ptr<ImpactFrictionMap> &impactFrictionMap();
  ImpactFrictionMap *impactFrictionMapPointer();

  bool frictionIsEnabled() const;

  void serialize(std::ostream &output_stream) const;
  void deserialize(std::istream &input_stream);

private:
  unsigned m_iteration;
  Rational<std::intmax_t> m_dt;
  std::unique_ptr<UnconstrainedMap> m_unconstrained_map;
  std::unique_ptr<ImpactOperator> m_impact_operator;
  std::unique_ptr<FrictionSolver> m_friction_solver;
  std::unique_ptr<ImpactFrictionMap> m_impact_friction_map;
  std::unique_ptr<ImpactMap> m_impact_map;
  scalar m_cor;
  scalar m_mu;
  bool m_reduce_bandwidth;
  PythonScripting m_python_scripting;
};

#endif

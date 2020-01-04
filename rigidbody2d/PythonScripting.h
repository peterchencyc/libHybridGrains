// PythonScripting.h
//
// Breannan Smith
// Last updated: 09/10/2015

#ifndef PYTHON_SCRIPTING_H
#define PYTHON_SCRIPTING_H

#ifdef USE_PYTHON
#include "scisim/PythonObject.h"
#include <Python.h>
#endif

#include "scisim/ScriptingCallback.h"

class RigidBody2DState;

class PythonScripting final : public ScriptingCallback {

public:
  PythonScripting();
  PythonScripting(const std::string &path, const std::string &module_name);
  explicit PythonScripting(std::istream &input_stream);

  virtual ~PythonScripting() override = default;

  friend void swap(PythonScripting &first, PythonScripting &second);

#ifdef USE_PYTHON
  static void initializeCallbacks();
#endif

  void setState(RigidBody2DState &state);
  // void setInitialIterate( unsigned& initial_iterate );
  void forgetState();

  void serialize(std::ostream &output_stream) const;

  virtual std::string name() const override;
  virtual std::string path() const override;

private:
  void intializePythonCallbacks();

  virtual void restitutionCoefficient(
      const std::vector<std::unique_ptr<Constraint>> &active_set,
      VectorXs &cor) override;

  virtual void frictionCoefficient(
      const std::vector<std::unique_ptr<Constraint>> &active_set,
      VectorXs &mu) override;

  virtual void startOfSim() override;

  virtual void endOfSim() override;

  virtual void startOfStep(const unsigned next_iteration,
                           const Rational<std::intmax_t> &dt) override;

  virtual void endOfStep(const unsigned next_iteration,
                         const Rational<std::intmax_t> &dt) override;

  std::string m_path;
  std::string m_module_name;
#ifdef USE_PYTHON
  PythonObject m_loaded_module;
  PythonObject m_loaded_start_of_sim_callback;
  PythonObject m_loaded_end_of_sim_callback;
  PythonObject m_loaded_start_of_step_callback;
  PythonObject m_loaded_end_of_step_callback;
  PythonObject m_loaded_friction_coefficient_callback;
  PythonObject m_loaded_restitution_coefficient_callback;
#endif
};

#endif

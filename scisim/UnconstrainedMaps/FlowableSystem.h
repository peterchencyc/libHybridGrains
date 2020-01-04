#ifndef FLOWABLE_SYSTEM_H
#define FLOWABLE_SYSTEM_H

#include "scisim/Math/MathDefines.h"

class FlowableSystem {

public:
  // TODO: Change to unsigned
  // Returns the number of degrees of freedom in this system
  virtual int nqdofs() const = 0;
  virtual int nvdofs() const = 0;
  virtual unsigned numVelDoFsPerBody() const = 0;
  virtual unsigned ambientSpaceDimensions() const = 0;
  unsigned numBodies() const;

  // True if the ith object in the simulation is 'kinematically scripted'
  // TODO: Change parameter to unsigned
  virtual bool isKinematicallyScripted(const int i) const = 0;

  // Given positions and velocities, computes the force acting on the DoFs,
  // overwriting the contents of F
  virtual void computeForce(const VectorXs &q, const VectorXs &v,
                            const scalar &t, VectorXs &F) = 0;

  // Sets forces on fixed degrees of freedom to 0
  virtual void zeroOutForcesOnFixedBodies(VectorXs &F) const;

  // Updates the configuration at q0 to q1 using v0 and dt with a linear update.
  virtual void linearInertialConfigurationUpdate(const VectorXs &q0,
                                                 const VectorXs &v0,
                                                 const scalar &dt,
                                                 VectorXs &q1) const = 0;

  // Returns the mass matrix as a sparse matrix
  virtual const SparseMatrixsc &M() const = 0;
  // Returns the inverse mass matrix as a sparse matrix
  virtual const SparseMatrixsc &Minv() const = 0;

  // Returns the reference configuration mass matrix as a sparse matrix
  virtual const SparseMatrixsc &M0() const = 0;
  // Returns the reference configuration inverse mass matrix as a sparse matrix
  virtual const SparseMatrixsc &Minv0() const = 0;

  // For the given velocity and the system's current configuration and mass,
  // computes the momentum
  virtual void computeMomentum(const VectorXs &v, VectorXs &p) const = 0;
  // For the given velocity and the system's current configuration and mass,
  // computes the angular momentum
  virtual void computeAngularMomentum(const VectorXs &v, VectorXs &L) const = 0;

  virtual const std::vector<bool> &fixed() const;

protected:
  FlowableSystem() = default;
  FlowableSystem(const FlowableSystem &) = default;
  FlowableSystem &operator=(const FlowableSystem &other) = default;
  FlowableSystem(FlowableSystem &&) noexcept = default;
  FlowableSystem &operator=(FlowableSystem &&other) noexcept = default;
  virtual ~FlowableSystem() = 0;
};

#endif

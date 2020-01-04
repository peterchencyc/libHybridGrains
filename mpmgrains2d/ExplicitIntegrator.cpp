#include "ExplicitIntegrator.h"

#include "SimulationState.h"

void ExplicitIntegrator::flow_zero_phase(SimulationState &state) {
  // Wipe the grid velocity, deformation gradient, velocity, and mass
  state.physics_grid.clearRasterizedData();

  // Rasterize the point masses to the grid
  state.physics_grid.rasterizePointMasses(state.material_points,
                                          state.basis_functions);

  // Rasterize the point momenta to the grid
  state.physics_grid.rasterizePointMomentum(state.material_points,
                                            state.basis_functions);
}

void ExplicitIntegrator::flow_first_phase(const scalar &dt,
                                          SimulationState &state) {
  // Compute the stress at material points
  state.material_points.computeHyperelasticCauchyStress(state.shear_modulus,
                                                        state.bulk_modulus);

  // Compute the internal elastic forces on the grid
  state.physics_grid.computeForces(state.material_points, state.basis_functions,
                                   state.near_earth_gravity);

  // Update the grid momentum
  state.physics_grid.updateMomentum(dt);

  // Grid based collision response (Stomakhin "A material point method for snow
  // simulation" Section 8)
  state.physics_grid.resolvePlaneCollisions(state.static_planes);
}

void ExplicitIntegrator::flow_second_phase(const scalar &dt,
                                           SimulationState &state) {
  // this step is intentionally put in the second phase rather than the first
  // phase (for Lagrange Multiplier Coupling)
  state.physics_grid.lumpedMassVelocityAndAccelerationUpdate(dt);

  state.material_points.computeVelGrad(
      state.basis_functions, state.physics_grid, state.physics_grid.velocity);
  state.material_points.elasticPrediction(dt, state.material_points.vel_grad);
  state.material_points.plasticCorrection(
      state.bulk_modulus, state.shear_modulus, state.Drucker_Prager_alpha);

  // Transfer the update from the grid to the point velocities
  state.material_points.updateVelocities(state.basis_functions,
                                         state.physics_grid, state.alpha, dt);

  // Update point positions
  state.material_points.updatePositions(state.basis_functions,
                                        state.physics_grid, dt);

  // Point based collision response (Stomakhin "A material point method for snow
  // simulation" Section 8)
  state.material_points.resolvePlaneCollisions(state.static_planes);

  assert(state.isConsistent());
}

void ExplicitIntegrator::flow(const scalar &dt, SimulationState &state) {
  flow_zero_phase(state);
  flow_first_phase(dt, state);
  flow_second_phase(dt, state);
}

void MPMGrains2DSim::flow_zero_phase(SimulationState &state) {
  ExplicitIntegrator::flow_zero_phase(state);
}

void MPMGrains2DSim::flow_first_phase(const scalar &dt,
                                      SimulationState &state) {
  ExplicitIntegrator::flow_first_phase(dt, state);
}

void MPMGrains2DSim::flow_second_phase(const scalar &dt,
                                       SimulationState &state) {
  ExplicitIntegrator::flow_second_phase(dt, state);
}

void MPMGrains2DSim::flow(const scalar &dt, SimulationState &state) {
  ExplicitIntegrator::flow(dt, state);
}

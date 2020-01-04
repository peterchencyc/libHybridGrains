#ifndef COUPLING_2D_H
#define COUPLING_2D_H

#include <map>
#include <memory>
#include <vector>

#include "DiscreteIntegrator.h"
#include "MPMIntegrator.h"
#include "ZoneTools.h"

#include "scisim/Math/MathDefines.h"
#include "scisim/Math/Rational.h"
#include "scisim/Timer/TimeUtils.h"

class HybridCoupling2D final {

public:
  //  grain-in-cell and particle-in-cell: rasterize DEM grain info onto MPM grid
  void rasterizeGrainMass(SimulationState &continuum_sim,
                          RigidBody2DSim &discrete_sim,
                          const bool allowCouplingPartiallyFilledCell,
                          const bool use_pre_position);
  void rasterizeGrainMom(const SimulationState &continuum_sim,
                         const RigidBody2DSim &discrete_sim,
                         const bool allowCouplingPartiallyFilledCell,
                         const bool use_pre_position);
  void transferToGrain(const std::unique_ptr<BasisFunctions> &in_SF,
                       const PhysicsGrid &grid, const MaterialPoints &points,
                       RigidBody2DState &dstate, RigidBody2DSim &discrete_sim,
                       const bool use_pre_position);

  // Resizes the coupling system. The size changes according to the number of
  // elements (which would change once we have resampling), and the size of the
  // packed tables.
  void resizeSystem(RigidBody2DSim &discrete_sim,
                    SimulationState &continuum_sim, VectorXs &lambda);

  // Copies the pre positions (excluding the rotation related quantity R) of the
  // grains.
  void copyPrePositionsFull(RigidBody2DSim &discrete_sim);

  void updateVelocitiesOnGrid(RigidBody2DSim &discrete_sim,
                              SimulationState &continuum_sim);

  // Perform a full update for the grain positions.
  void updatePositions_node_node(RigidBody2DSim &discrete_sim, scalar dt);

  const VectorXs &getDEMGridMasses() const;

  const VectorXi &getPackedToUnpackedMPM() const;

private:
  VectorXi packed_to_unpacked_dem;
  VectorXi unpacked_to_packed_dem;
  VectorXi packed_to_unpacked_mpm;
  VectorXi unpacked_to_packed_mpm;

  TimingTools m_timing_tools;

  // Data rasterized from DEM grains
  VectorXs dem_mass_grid;
  Matrix2Xsc dem_momentum_grid;

  Matrix2Xsc v_grid_after;
  Matrix2Xsc at_grid_after;

  // for future PREDICTION_CORRECTION_ADVECTION_NODE_NODE
  VectorXs dem_pos_before;
};

#endif

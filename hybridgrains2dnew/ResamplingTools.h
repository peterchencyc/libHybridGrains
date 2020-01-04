#ifndef RESAMPLING_TOOLS_H
#define RESAMPLING_TOOLS_H

#include "scisim/Math/MathDefines.h"

#include "AvoidAvoid.h"
#include "Grid.h"
#include "ZoneTools.h"
#include "mpmgrains2d/BasisFunctions.h"
#include "mpmgrains2d/SimulationState.h"
#include "rigidbody2d/PenaltyImpactFrictionMap.h"
#include "rigidbody2d/RigidBody2DSim.h"
#include <random>

class ImpactFrictionMap;

class ResamplingTools final {

public:
  void init(const scalar &dem_r_mean, const scalar &dem_r_std);

  void deleteParticlesUsingZoneIndicators(ZoneTools &zone_tools,
                                          RigidBody2DState &discrete_state,
                                          SimulationState &continuum_state,
                                          ImpactFrictionMap *ifmap);

  void avoidAVoidMPMparticlesInsideZonesWithTentativeQuantities(
      SimulationState &continuum_state,
      const std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
      int &numPointsBeforeSampling, int &numPointsAfterSampling);

  void avoidAVoidDEMparticlesInsideZonesWithTentativeQuantities(
      RigidBody2DState &discrete_state,
      const std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
      const scalar &in_rho_dem, ImpactFrictionMap *ifmap,
      int &numGrainsBeforeSampling, int &numGrainsAfterSampling);

  void clearHomogenizedStress(const SimulationState &continuum_state);
  void homogenizeStress(const RigidBody2DState &discrete_state,
                        const SimulationState &continuum_state,
                        PenaltyImpactFrictionMap *pifmap);
  void finalizeHomogenizedStress(const SimulationState &continuum_state,
                                 const scalar &factor);

  void clearHomogenizedVelocity(const SimulationState &continuum_state);
  void homogenizeVelocity(const RigidBody2DState &discrete_state,
                          const SimulationState &continuum_state);
  void finalizeHomogenizedVelocity(const SimulationState &continuum_state,
                                   const scalar &factor);

  void determineMPMQuantitiesForHybridZones(
      SimulationState &continuum_state,
      const std::vector<Vector2u> &zone_indices,
      const int numPointsBeforeSampling, const int numPointsAfterSampling,
      const RigidBody2DState &discrete_state, const int numGrainsBeforeSampling,
      const CUniformGridH2D &ugrd_continuum,
      const CUniformGridH2D &ugrd_discrete, const scalar &cell_volume,
      bool newly_inserted_continuum_stress_free, bool homogenize_stress,
      PenaltyImpactFrictionMap *pifmap, bool homogenize_velocity,
      const scalar kappa, const scalar mu);

  void determineDEMQuantitiesForHybridZones(
      RigidBody2DState &discrete_state, const int numGrainsBeforeSampling,
      const int numGrainsAfterSampling, const SimulationState &continuum_state,
      const int numPointsBeforeSampling, const CUniformGridH2D &ugrd_continuum,
      const CUniformGridH2D &ugrd_discrete);

  void determineMPMQuantitiesForInnerContinuumZones(
      SimulationState &continuum_state, const int numPointsBeforeSampling,
      const int numPointsAfterSampling, const CUniformGridH2D &ugrd,
      bool newly_inserted_continuum_stress_free);

private:
  void deleteMPMparticlesUsingZoneIndicators(
      ZoneTools &zone_tools, SimulationState &continuum_state,
      std::vector<int> &mpm_old_to_new_id_map);
  void deleteDEMparticlesUsingZoneIndicators(
      ZoneTools &zone_tools, RigidBody2DState &discrete_state,
      std::vector<int> &dem_old_to_new_id_map, ImpactFrictionMap *ifmap);
  void reassignMassAndVolumesForHybridAndNewlyGeneratedContinuumZone(
      SimulationState &continuum_state,
      const std::vector<Vector2u> &zone_indices,
      const CUniformGridH2D &ugrd_continuum, const scalar cell_volume);

  void reassignJVelAndBebarForHybridMaterialPoint(
      SimulationState &continuum_state, const int pnt_idx, const scalar sigma,
      const scalar max_dist, const int numPointsBeforeSampling,
      const int numPointsAfterSampling, const RigidBody2DState &discrete_state,
      const int numGrainsBeforeSampling, const CUniformGridH2D &ugrd_continuum,
      const CUniformGridH2D &ugrd_discrete,
      bool newly_inserted_continuum_stress_free, bool homogenize_stress,
      PenaltyImpactFrictionMap *pifmap, bool homogenize_velocity,
      const scalar kappa, const scalar mu);

  void computeFieldsForInterpolation(SimulationState &continuum_state,
                                     const int numPointsBeforeSampling);

  void reassignJVelAndBebarForInternalMaterialPoint(
      SimulationState &continuum_state, const int pnt_idx,
      bool newly_inserted_continuum_stress_free);

  double constantRadiusSampler(const double &r);
  double randomRadiusSampler();

  CAvoidAVoid m_avoid_a_void;

  std::default_random_engine m_generator;
  std::normal_distribution<double> m_dist;
  scalar m_dem_r_mean;
  scalar m_dem_r_std;

  std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>
      m_homogenized_stress;
  std::vector<Vector2s> m_homogenized_momentum;
  std::vector<scalar> m_homogenized_mass;
  std::vector<Vector2s> m_homogenized_velocity;

  Grid *m_J_field;
  Grid *m_Vel_field;
  Grid *m_Bebar_field;
};

#endif

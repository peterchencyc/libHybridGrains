#ifndef ZONE_TOOLS_H
#define ZONE_TOOLS_H

#include <Eigen/StdVector>

#include "scisim/Math/MathDefines.h"

#include "mpmgrains2d/BasisFunctions.h"
#include "mpmgrains2d/SimulationState.h"
#include "rigidbody2d/RigidBody2DSim.h"

#include "CircleRectangleIntersection.h"
#include "distgrid.h"
#include "levelsetutil.h"
#include "uniformgrid.h"

enum ZoneType { ZT_CONTINUUM, ZT_DISCRETE, ZT_HYBRID };

class ZoneTools {
public:
  void init(const SimulationState &continuum_state,
            const scalar level_set_cell_width,
            bool allow_direct_transitions_between_discrete_and_continuum);

  scalar demWeight(const VectorXs &q);
  scalar mpmWeight(const VectorXs &q);

  // void updatePackingFraction( RigidBody2DState& discrete_state );
  void updatePackingFraction(RigidBody2DState &discrete_state,
                             SimulationState &continuum_state,
                             const scalar phi_window_size,
                             const int phi_samples_per_cell_side);
  void updateDistanceField(scalar phi_threshold);
  void computeRZones(scalar rzone_level_set, scalar rzone_half_thickness);

  bool isParticleInDiscreteZone_MPMIndicator(const Vector2s &x) const;
  bool isParticleInContinuumZone_MPMIndicator(const Vector2s &x) const;
  bool isParticleInHybridZone_MPMIndicator(const Vector2s &x) const;

  std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &getRZones();
  std::vector<Vector2u> &getRZoneIndices();

  std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &
  getDirectTransitionZonesFromDiscreteToContinuum();
  std::vector<Vector2u> &
  getDirectTransitionZonesFromDiscreteToContinuumIndices();

  void determineMPMParticlesToDelete(SimulationState &continuum_state,
                                     std::vector<unsigned> &id_list_to_del);
  void determineDEMGrainsToDelete(RigidBody2DState &discrete_state,
                                  std::vector<unsigned> &id_list_to_del);

  void identifyInnerContinuumRegion(
      std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
      std::vector<Vector2u> &zone_indices);
  void identifyInnerContinuumRegionExcludingDirectTransitionZones(
      std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
      std::vector<Vector2u> &zone_indices);

  CUniformGridH2D *
  setupUniformGridElementsAsPoints(RigidBody2DState &discrete_state);
  CUniformGridH2D *
  setupUniformGridElementsAsPoints(SimulationState &continuum_state);

  bool allowDirectTransitionsBetweenDiscreteAndContinuum() const;

private:
  scalar computePackingFractionForBox(const Vector2s &box_min,
                                      const Vector2s &box_max, Matrix2Xsc &pos,
                                      VectorXs &r, CUniformGridH2D &ugrd,
                                      VectorXi &mail_box, const int test_id);
  void estimateLevelSetLevelZoneIndicators(scalar rzone_level_set,
                                           scalar rzone_half_thickness);
  void sampleMPMLevelZoneIndicators();

  DistGrid m_dist_grid;
  Eigen::VectorXd m_packing_fractions;
  CrossingPoints m_boundary_samples;

  Eigen::VectorXi m_zone_indicators_level_set;
  Eigen::VectorXi m_zone_indicators_mpm;
  Eigen::VectorXi m_prev_zone_indicators_mpm;

  std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> m_RZones;
  std::vector<Vector2u> m_RZoneIndices;
  std::vector<scalar> m_RZonePackingFractions;

  std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>>
      m_DirectTransitionZonesFromDiscreteToContinuum;
  std::vector<Vector2u> m_DirectTransitionZonesFromDiscreteToContinuumIndices;

  std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>>
      m_RZones_level_set_res;

  Array2u m_mpm_cell_count;
  Vector2s m_mpm_grid_min;
  scalar m_mpm_cell_width;

  bool m_allow_direct_transitions_between_discrete_and_continuum;
};

#endif

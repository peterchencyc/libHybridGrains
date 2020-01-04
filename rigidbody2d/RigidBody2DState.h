#ifndef RIGID_BODY_2D_STATE
#define RIGID_BODY_2D_STATE

#include "PlanarPortal.h"
#include "RigidBody2DForce.h"
#include "RigidBody2DGeometry.h"
#include "RigidBody2DStaticDrum.h"
#include "RigidBody2DStaticPlane.h"
#include "scisim/Math/MathDefines.h"

class ImpactFrictionMap;

class RigidBody2DState final {

public:
  RigidBody2DState();
  RigidBody2DState(
      const VectorXs &q, const VectorXs &v, const VectorXs &m,
      const std::vector<bool> &fixed, const VectorXu &geometry_indices,
      const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
      const std::vector<std::unique_ptr<RigidBody2DForce>> &forces,
      const std::vector<RigidBody2DStaticPlane> &planes,
      const std::vector<PlanarPortal> &planar_portals,
      const std::vector<RigidBody2DStaticDrum> &drums);

  RigidBody2DState(const RigidBody2DState &rhs);
  RigidBody2DState(RigidBody2DState &&) = default;

  RigidBody2DState &operator=(const RigidBody2DState &rhs);
  RigidBody2DState &operator=(RigidBody2DState &&) = default;

  unsigned nbodies() const;
  unsigned numBodies() const;
  unsigned numSimulatedBodies() const;
  unsigned numGeometryInstances() const;

  scalar computeTotalSimulatedMass() const;
  Vector2s computeTotalSimulatedMomentum() const;

  VectorXs &q();
  const VectorXs &q() const;

  VectorXs &q0();
  const VectorXs &q0() const;

  VectorXs &v();
  const VectorXs &v() const;

  const SparseMatrixsc &M() const;

  const SparseMatrixsc &Minv() const;

  const scalar &getTotalMass(const unsigned bdy_idx) const;

  const scalar &m(const unsigned bdy_idx) const;
  const scalar &I(const unsigned bdy_idx) const;

  bool fixed(const int idx) const;
  bool alwaysFixed(const int idx) const;

  void fixBody(const unsigned idx);
  void unfixBody(const unsigned idx);

  // Adds a new body at the end of the state vector
  void addBody(const Vector2s &x, const scalar &theta, const Vector2s &v,
               const scalar &omega, const scalar &rho, const unsigned geo_idx,
               const bool fixed, ImpactFrictionMap *ifmap);

  void addBody(const Vector2s &x, const Vector2s &x0, const scalar &theta,
               const Vector2s &v, const scalar &omega, const scalar &rho,
               const unsigned geo_idx, const bool fixed,
               ImpactFrictionMap *ifmap);
  void
  addBodies(const std::vector<Vector2s> &x, const std::vector<Vector2s> &x0,
            const std::vector<scalar> &theta, const std::vector<Vector2s> &v,
            const std::vector<scalar> &omega, const std::vector<scalar> &rho,
            const std::vector<unsigned> &geo_idx,
            const std::vector<bool> &fixed, ImpactFrictionMap *ifmap);

  // Removes bodies that intersect the given boxes
  void removeBodiesIntersectingBoxes(
      const std::vector<std::pair<Array2s, Array2s>> &boxes,
      ImpactFrictionMap *ifmap);

  // Removes bodies at the given indices from the simulation
  void removeBodies(const Eigen::Ref<const VectorXu> &indices,
                    ImpactFrictionMap *ifmap);
  void removeBodies(const std::vector<unsigned> &indices,
                    ImpactFrictionMap *ifmap);
  void removeBodies(const std::vector<unsigned> &indices,
                    std::vector<int> &old_to_new_idx_ref_map,
                    ImpactFrictionMap *ifmap);

  // Adds a new circle geometry instance to the back of the geometry vector, and
  // return the geometry index
  int addCircleGeometry(const scalar &r);

  // Removes geometry instances from the simulation
  void removeGeometry(const Eigen::Ref<const VectorXu> &indices);
  void removeGeometry(const std::vector<unsigned> &indices);

  // TODO: Remove these
  std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry();
  const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry() const;
  VectorXu &geometryIndices();
  const VectorXu &geometryIndices() const;

  const std::unique_ptr<RigidBody2DGeometry> &
  bodyGeometry(const unsigned bdy_idx) const;

  const std::vector<std::unique_ptr<RigidBody2DForce>> &forces() const;

  std::vector<RigidBody2DStaticPlane> &planes();
  const std::vector<RigidBody2DStaticPlane> &planes() const;
  void deleteStaticPlane(const unsigned plane_index);

  std::vector<PlanarPortal> &planarPortals();
  const std::vector<PlanarPortal> &planarPortals() const;

  const std::vector<RigidBody2DStaticDrum> &drums() const;
  std::vector<RigidBody2DStaticDrum> &drums();

  // Computes a bounding box around the system
  Array4s computeBoundingBox() const;

  void getAllCircleBodies(Matrix2Xsc &pos, VectorXs &radii,
                          VectorXu &bdy_idx) const;
  void getAllNonFixedCircleBodies(Matrix2Xsc &pos, VectorXs &radii,
                                  VectorXu &bdy_idx) const;

  // Serialization and deserialization
  void serialize(std::ostream &output_stream) const;
  void deserialize(std::istream &input_stream);

  const scalar &hybridFactor(const unsigned idx) const;
  scalar &hybridFactor(const unsigned idx);

  const VectorXs &hybridFactors() const;

  scalar totalMass() const;
  scalar totalSimulatedMass() const;

  const scalar &distanceField(const unsigned idx) const;
  scalar &distanceField(const unsigned idx);

  void modifyMass(const unsigned pnt_idx, const scalar &new_mass);

  void setCurrentMassAsInitialMass();
  void updateCurrentMassUsingMassWeights();
  void updateMassWeights(std::function<scalar(const VectorXs &)> &weight_func);

  const scalar &massWeight(const unsigned idx) const;
  scalar &massWeight(const unsigned idx);
  const VectorXs &massWeights() const;

  scalar computeMaximumSpeed() const;

  int getUniqueBodyIndex(const int idx) const;
  std::vector<int> &uniqueBodyIDs();

#ifndef NDEBUG
  void checkStateConsistency();
#endif

private:
  // Format: x0, y0, theta0, x1, y1, theta1, ...
  VectorXs m_q;
  VectorXs m_q0;
  // Format: vx0, vy0, omega0, vx1, vy1, omega1, ...
  VectorXs m_v;
  SparseMatrixsc m_M;
  SparseMatrixsc m_Minv;
  std::vector<bool> m_fixed;
  std::vector<bool> m_always_fixed;
  VectorXu m_geometry_indices;
  std::vector<std::unique_ptr<RigidBody2DGeometry>> m_geometry;
  std::vector<std::unique_ptr<RigidBody2DForce>> m_forces;
  std::vector<RigidBody2DStaticPlane> m_planes;
  std::vector<PlanarPortal> m_planar_portals;
  std::vector<RigidBody2DStaticDrum> m_drums;
  // Unique indices for bodies that are maintained across additions or deletions
  int m_next_unique_id;
  std::vector<int> m_uique_ids;

  // For hybridization
  VectorXs m_homog_factor;
  VectorXs m_distance_field;
  SparseMatrixsc m_M0;
  VectorXs m_mass_weights;
};

#endif

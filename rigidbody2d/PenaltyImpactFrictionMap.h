#ifndef PENALTY_IMPACT_FRICTION_MAP_H
#define PENALTY_IMPACT_FRICTION_MAP_H

#include "scisim/ConstrainedMaps/ImpactFrictionMap.h"
#include "scisim/Utilities.h"

// #ifdef USE_HDF5
// class HDF5File;
// #endif

struct CircleHalfPlaneCollision final {
  CircleHalfPlaneCollision(const scalar &pen_depth_in, const Vector2s &n_in,
                           const Vector2s &delta_s_in, const Vector2s &vrel_in,
                           const Vector2s &r_in, const scalar &weight_in)
      : pen_depth(pen_depth_in), n(n_in), delta_s(delta_s_in), vrel(vrel_in),
        r(r_in), weight(weight_in) {
    assert(pen_depth <= 0.0);
    assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
    assert(std::fabs(n.dot(delta_s)) <= 1.0e-6);
    assert(weight >= 0.0);
    assert(weight <= 1.0);
  }

  scalar pen_depth;
  Vector2s n;
  Vector2s delta_s;
  Vector2s vrel;
  Vector2s r;
  scalar weight;

  Vector2s normal_force;
  Vector2s friction_force;
};

struct CircleDrumCollision final {
  CircleDrumCollision(const scalar &pen_depth_in, const Vector2s &n_in,
                      const Vector2s &delta_s_in, const Vector2s &vrel_in,
                      const Vector2s &r_in, const scalar &weight_in)
      : pen_depth(pen_depth_in), n(n_in), delta_s(delta_s_in), vrel(vrel_in),
        r(r_in), weight(weight_in) {
    assert(pen_depth <= 0.0);
    assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
    assert(std::fabs(n.dot(delta_s)) <= 1.0e-6);
    assert(weight >= 0.0);
    assert(weight <= 1.0);
  }

  scalar pen_depth;
  Vector2s n;
  Vector2s delta_s;
  Vector2s vrel;
  Vector2s r;
  scalar weight;

  Vector2s normal_force;
  Vector2s friction_force;
};

struct CircleCircleCollision final {
  CircleCircleCollision(const scalar &pen_depth_in, const Vector2s &n_in,
                        const Vector2s &delta_s_in, const Vector2s &vrel_in,
                        const Vector2s &r0_in, const Vector2s &r1_in,
                        const scalar &weight_in,
                        const Vector2s &contact_point_in)
      : pen_depth(pen_depth_in), n(n_in), delta_s(delta_s_in), vrel(vrel_in),
        r0(r0_in), r1(r1_in), weight(weight_in),
        contact_point(contact_point_in) {
    assert(pen_depth <= 0.0);
    assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
    assert(std::fabs(n.dot(delta_s)) <= 1.0e-6);
    assert(weight >= 0.0);
    assert(weight <= 1.0);
  }

  scalar pen_depth;
  Vector2s n;
  Vector2s delta_s;
  Vector2s vrel;
  Vector2s r0;
  Vector2s r1;
  scalar weight;

  Vector2s contact_point;
  Vector2s normal_force;
  Vector2s friction_force;
};

struct KinematicCircleCircleCollision final {
  KinematicCircleCircleCollision(const scalar &pen_depth_in,
                                 const Vector2s &n_in,
                                 const Vector2s &delta_s_in,
                                 const Vector2s &vrel_in, const Vector2s &r0_in,
                                 const scalar &weight_in)
      : pen_depth(pen_depth_in), n(n_in), delta_s(delta_s_in), vrel(vrel_in),
        r0(r0_in), weight(weight_in) {
    assert(pen_depth <= 0.0);
    assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
    assert(std::fabs(n.dot(delta_s)) <= 1.0e-6);
    assert(weight >= 0.0);
    assert(weight <= 1.0);
  }

  scalar pen_depth;
  Vector2s n;
  Vector2s delta_s;
  Vector2s vrel;
  Vector2s r0;
  scalar weight;

  Vector2s normal_force;
  Vector2s friction_force;
};

struct KinematicCircleBoxCollision final {
  KinematicCircleBoxCollision(const scalar &pen_depth_in, const Vector2s &n_in,
                              const Vector2s &delta_s_in,
                              const Vector2s &vrel_in, const Vector2s &r0_in,
                              const scalar &weight_in)
      : pen_depth(pen_depth_in), n(n_in), delta_s(delta_s_in), vrel(vrel_in),
        r0(r0_in), weight(weight_in) {
    assert(pen_depth <= 0.0);
    assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
    assert(std::fabs(n.dot(delta_s)) <= 1.0e-6);
    assert(weight >= 0.0);
    assert(weight <= 1.0);
  }

  scalar pen_depth;
  Vector2s n;
  Vector2s delta_s;
  Vector2s vrel;
  Vector2s r0;
  scalar weight;

  Vector2s normal_force;
  Vector2s friction_force;
};

struct TeleportedCircleCircleCollision final {
  TeleportedCircleCircleCollision(const scalar &pen_depth_in,
                                  const Vector2s &n_in,
                                  const Vector2s &delta_s_in,
                                  const Vector2s &vrel_in,
                                  const Vector2s &r0_in, const Vector2s &r1_in,
                                  const scalar &weight_in)
      : pen_depth(pen_depth_in), n(n_in), delta_s(delta_s_in), vrel(vrel_in),
        r0(r0_in), r1(r1_in), weight(weight_in) {
    assert(pen_depth <= 0.0);
    assert(std::fabs(n.norm() - 1.0) <= 1.0e-6);
    assert(std::fabs(n.dot(delta_s)) <= 1.0e-6);
    assert(weight >= 0.0);
    assert(weight <= 1.0);
  }

  scalar pen_depth;
  Vector2s n;
  Vector2s delta_s;
  Vector2s vrel;
  Vector2s r0;
  Vector2s r1;
  scalar weight;

  Vector2s normal_force;
  Vector2s friction_force;
};

template <typename T> class CollisionTypeCache final {

public:
  CollisionTypeCache() : CollisionTypeCache(0) {}

  CollisionTypeCache(const int num_bodies)
      : m_num_bodies(num_bodies), m_num_collisions(0), m_entries(num_bodies) {}

  void clear() {
    for (int idx = 0; idx < m_num_bodies; idx++) {
      m_entries[idx].clear();
    }
    m_num_collisions = 0;
  }

  bool empty() const { return m_num_collisions == 0; }

  int numBodies() const { return m_num_bodies; }

  int numCollisions() const { return m_num_collisions; }

  T *find(const int idx0, const int idx1) {
    if (empty()) {
      return nullptr;
    }
    assert(idx0 >= 0);
    assert(idx0 < m_num_bodies);
    for (std::pair<int, T> &entry : m_entries[idx0]) {
      if (entry.first == idx1) {
        return &entry.second;
      }
    }
    return nullptr;
  }

  void insert(const int idx0, const int idx1, const T &entry) {
    assert(find(idx0, idx1) == nullptr);
    m_entries[idx0].emplace_back(idx1, entry);
    m_num_collisions++;
  }

  std::vector<std::vector<std::pair<int, T>>> &collisions() {
    return m_entries;
  }

  std::vector<std::pair<int, T>> &collisions(const int bdy_idx) {
    assert(bdy_idx >= 0);
    assert(bdy_idx < int(m_entries.size()));
    return m_entries[bdy_idx];
  }

  void resize(const int nbodies) {
    assert(nbodies >= 0);
    m_entries.resize(nbodies);
    m_num_bodies = nbodies;
  }

  void serialize(std::ostream &output_stream) const {
    Utilities::serializeBuiltInType(m_num_bodies, output_stream);
    Utilities::serializeBuiltInType(m_num_collisions, output_stream);
    assert(int(m_entries.size()) == m_num_bodies);
    for (const std::vector<std::pair<int, T>> &entries : m_entries) {
      Utilities::serializeBuiltInType(int(entries.size()), output_stream);
      for (const std::pair<int, T> &entry : entries) {
        Utilities::serializeBuiltInType(entry.first, output_stream);
        entry.second.serialize(output_stream);
      }
    }
  }

  void deserialize(std::istream &input_stream) {
    m_num_bodies = Utilities::deserialize<int>(input_stream);
    m_num_collisions = Utilities::deserialize<int>(input_stream);
    m_entries.resize(m_num_bodies);
    for (int bdy_idx = 0; bdy_idx < m_num_bodies; bdy_idx++) {
      const int num_entries = Utilities::deserialize<int>(input_stream);
      m_entries[bdy_idx].reserve(num_entries);
      for (int entry_idx = 0; entry_idx < num_entries; entry_idx++) {
        int idx1 = Utilities::deserialize<int>(input_stream);
        T val(input_stream);
        m_entries[bdy_idx].emplace_back(idx1, val);
      }
    }
  }

  void relabelEntries(const std::vector<int> &relabels) {
    m_num_collisions = 0;
    for (std::vector<std::pair<int, T>> &row : m_entries) {
      for (std::pair<int, T> &entry : row) {
        assert(entry.first >= 0);
        assert(entry.first < int(relabels.size()));
        assert(relabels[entry.first] <= entry.first);
        entry.first = relabels[entry.first];
      }
      row.erase(std::remove_if(
                    row.begin(), row.end(),
                    [](const std::pair<int, T> &v) { return v.first == -1; }),
                row.end());
      m_num_collisions += row.size();
    }
  }

private:
  int m_num_bodies;
  int m_num_collisions;
  std::vector<std::vector<std::pair<int, T>>> m_entries;
};

struct CollisionCache final {
  // First index: circle
  // Second index: half plane
  CollisionTypeCache<CircleHalfPlaneCollision> circle_half_plane_collisions;

  // First index: circle
  // Second index: drum
  CollisionTypeCache<CircleDrumCollision> circle_drum_collisions;

  // First index: circle 0
  // Second index: circle 1
  CollisionTypeCache<CircleCircleCollision> circle_circle_collisions;

  // First index: circle 0
  // Second index: circle 1
  CollisionTypeCache<KinematicCircleCircleCollision>
      kinematic_circle_circle_collisions;

  // First index: body 0
  // Second index: body1
  CollisionTypeCache<KinematicCircleBoxCollision>
      kinematic_box_circle_collisions;

  // First index: body 0
  // Second index: body1
  CollisionTypeCache<TeleportedCircleCircleCollision>
      teleported_circle_circle_collisions;

  void clear() {
    circle_half_plane_collisions.clear();
    circle_drum_collisions.clear();
    circle_circle_collisions.clear();
    kinematic_circle_circle_collisions.clear();
    kinematic_box_circle_collisions.clear();
    teleported_circle_circle_collisions.clear();
  }

  void resizeNumBodies(const int nbodies) {
    circle_half_plane_collisions.resize(nbodies);
    circle_drum_collisions.resize(nbodies);
    circle_circle_collisions.resize(nbodies);
    kinematic_circle_circle_collisions.resize(nbodies);
    kinematic_box_circle_collisions.resize(nbodies);
    teleported_circle_circle_collisions.resize(nbodies);
  }

  bool empty() const {
    return circle_half_plane_collisions.empty() &&
           circle_drum_collisions.empty() && circle_circle_collisions.empty() &&
           kinematic_circle_circle_collisions.empty() &&
           kinematic_box_circle_collisions.empty() &&
           teleported_circle_circle_collisions.empty();
  }

  int numBodies() const {
    assert(circle_half_plane_collisions.numBodies() ==
           circle_drum_collisions.numBodies());
    assert(circle_half_plane_collisions.numBodies() ==
           circle_circle_collisions.numBodies());
    assert(circle_half_plane_collisions.numBodies() ==
           kinematic_circle_circle_collisions.numBodies());
    assert(circle_half_plane_collisions.numBodies() ==
           kinematic_box_circle_collisions.numBodies());
    assert(circle_half_plane_collisions.numBodies() ==
           teleported_circle_circle_collisions.numBodies());
    return circle_half_plane_collisions.numBodies();
  }

  int numTotalCollisions() const {
    return circle_half_plane_collisions.numCollisions() +
           circle_drum_collisions.numCollisions() +
           circle_circle_collisions.numCollisions() +
           kinematic_circle_circle_collisions.numCollisions() +
           kinematic_box_circle_collisions.numCollisions() +
           teleported_circle_circle_collisions.numCollisions();
  }

  void relabelEntries(const std::vector<int> &relabels) {
    circle_circle_collisions.relabelEntries(relabels);
    kinematic_circle_circle_collisions.relabelEntries(relabels);
    kinematic_box_circle_collisions.relabelEntries(relabels);
    teleported_circle_circle_collisions.relabelEntries(relabels);
  }

  // TODO: Make an 'is consistent' method here
};

class PenaltyImpactFrictionMap final : public ImpactFrictionMap {

public:
  PenaltyImpactFrictionMap(const scalar &kn, const scalar &kt,
                           const scalar &gamman, const scalar &gammat,
                           const scalar &mu);

  // explicit PenaltyImpactFrictionMap( std::istream& input_stream );

  PenaltyImpactFrictionMap(const PenaltyImpactFrictionMap &) = delete;
  PenaltyImpactFrictionMap(PenaltyImpactFrictionMap &&) = delete;
  PenaltyImpactFrictionMap &
  operator=(const PenaltyImpactFrictionMap &) = delete;
  PenaltyImpactFrictionMap &operator=(PenaltyImpactFrictionMap &&) = delete;

  virtual ~PenaltyImpactFrictionMap() override = default;

  // NB: This integrator ignores most of these parameters
  virtual void flow(ScriptingCallback &call_back, FlowableSystem &fsys,
                    ConstrainedSystem &csys, UnconstrainedMap &umap,
                    FrictionSolver &friction_solver, const unsigned iteration,
                    const scalar &dt, const scalar &CoR_default,
                    const scalar &mu_default, const bool reduce_bandwidth,
                    const VectorXs &q0, const VectorXs &v0, VectorXs &q1,
                    VectorXs &v1) override;

  virtual void flowWithWeights(ScriptingCallback &call_back,
                               FlowableSystem &fsys, ConstrainedSystem &csys,
                               UnconstrainedMap &umap,
                               FrictionSolver &friction_solver,
                               const unsigned iteration, const scalar &dt,
                               const scalar &CoR, const scalar &mu_default,
                               const bool reduce_bandwidth, const VectorXs &q0,
                               const VectorXs &v0, const VectorXs &w,
                               VectorXs &q1, VectorXs &v1) override;

  virtual void resetCachedData() override;

  virtual void serialize(std::ostream &output_stream) const override;

  virtual std::string name() const override;

  virtual void exportForcesNextStep(HDF5File &output_file) override;
  virtual void saveImpulsesNextStep(CollisionImpulses &impulses) override;

  virtual std::unique_ptr<ImpactFrictionMap> clone() const override;

  virtual bool stabilize() const override;

  virtual void enlargeCache(const int new_size) override;
  virtual void deleteCacheEntries(const VectorXs &q) override;

  CollisionCache &getCollisionCache();

private:
  void updateCollisionCache(const bool reduce_bandwidth, const scalar &dt,
                            ConstrainedSystem &csys, const VectorXs &q,
                            const VectorXs &v);
  void updateCollisionCache(const bool reduce_bandwidth, const scalar &dt,
                            ConstrainedSystem &csys, const VectorXs &q,
                            const VectorXs &v, const VectorXs &w);
  void updateAnchorsAndAccumulateCollisionForces(const VectorXs &v,
                                                 VectorXs &F);
  // void updateAnchorsAndAccumulateCollisionForces( std::function< scalar(const
  // VectorXs&) >& weight_func, const VectorXs& v, VectorXs& F );

  scalar m_kn;
  scalar m_kt;
  scalar m_gamman;
  scalar m_gammat;
  scalar m_mu;

  CollisionCache m_collision_cache;
  CollisionCache m_old_cache;
};

#endif

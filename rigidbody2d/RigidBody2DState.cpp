#include "RigidBody2DState.h"

#include "scisim/ConstrainedMaps/ImpactFrictionMap.h"
#include "scisim/Math/MathUtilities.h"
#include "scisim/StringUtilities.h"
#include "scisim/Utilities.h"

#include "AnnulusGeometry.h"
#include "BoxGeometry.h"
#include "CircleBoxTools.h"
#include "CircleGeometry.h"
#include "NearEarthGravityForce.h"

#include <iostream>

static SparseMatrixsc generateM(const VectorXs &m) {
  SparseMatrixsc M{SparseMatrixsc::Index(m.size()),
                   SparseMatrixsc::Index(m.size())};
  M.reserve(SparseMatrixsc::Index(m.size()));
  for (int col = 0; col < m.size(); ++col) {
    M.startVec(col);
    const int row{col};
    M.insertBack(row, col) = m(col);
  }
  M.finalize();
  M.makeCompressed();
  return M;
}

static SparseMatrixsc generateMinv(const VectorXs &m) {
  SparseMatrixsc Minv{SparseMatrixsc::Index(m.size()),
                      SparseMatrixsc::Index(m.size())};
  Minv.reserve(SparseMatrixsc::Index(m.size()));
  for (int col = 0; col < m.size(); ++col) {
    Minv.startVec(col);
    const int row{col};
    Minv.insertBack(row, col) = 1.0 / m(col);
  }
  Minv.finalize();
  Minv.makeCompressed();
  return Minv;
}

#ifndef NDEBUG
void RigidBody2DState::checkStateConsistency() {
  assert(m_q.size() % 3 == 0);
  assert(m_q0.size() == m_q.size());
  assert(m_q.size() == m_v.size());

  assert(m_M.rows() == m_q.size());
  assert(m_M.cols() == m_q.size());
  assert(
      (Eigen::Map<const ArrayXs>{m_M.valuePtr(), m_M.nonZeros()} > 0.0).all());
  for (unsigned bdy_idx = 0; bdy_idx < m_q.size() / 3; ++bdy_idx) {
    assert(m_M.valuePtr()[3 * bdy_idx] == m_M.valuePtr()[3 * bdy_idx + 1]);
  }

  assert(m_Minv.rows() == m_q.size());
  assert(m_Minv.cols() == m_q.size());
  assert((Eigen::Map<const ArrayXs>{m_Minv.valuePtr(), m_Minv.nonZeros()} > 0.0)
             .all());
  for (unsigned bdy_idx = 0; bdy_idx < m_q.size() / 3; ++bdy_idx) {
    assert(m_Minv.valuePtr()[3 * bdy_idx] ==
           m_Minv.valuePtr()[3 * bdy_idx + 1]);
  }

  // Check that the product of M and M^{-1} is the identity
  {
    const SparseMatrixsc prod{m_M * m_Minv};
    assert(prod.nonZeros() == m_M.rows());
    assert(((Eigen::Map<const ArrayXs>{prod.valuePtr(), prod.nonZeros()} - 1.0)
                .abs() <= 1.0e-6)
               .all());
  }

  assert(static_cast<int>(m_fixed.size()) == m_q.size() / 3);
  assert(static_cast<int>(m_always_fixed.size()) == m_q.size() / 3);

  assert(static_cast<int>(m_geometry_indices.size()) == m_q.size() / 3);
  assert((m_geometry_indices.array() < unsigned(m_geometry.size())).all());

  assert(m_homog_factor.size() == m_q.size() / 3);
  assert(m_distance_field.size() == m_q.size() / 3);

  assert(m_M0.rows() == m_M.rows());
  assert(m_M0.cols() == m_M.cols());

  assert(m_mass_weights.size() == m_q.size() / 3);

  assert(m_uique_ids.size() == numBodies());
}
#endif

RigidBody2DState::RigidBody2DState()
    : m_q(), m_q0(), m_v(), m_M(), m_Minv(), m_fixed(), m_always_fixed(),
      m_geometry_indices(), m_geometry(), m_forces(), m_planes(),
      m_planar_portals(), m_drums(), m_next_unique_id(0), m_uique_ids(),
      m_homog_factor(), m_distance_field(), m_M0(), m_mass_weights() {}

RigidBody2DState::RigidBody2DState(
    const VectorXs &q, const VectorXs &v, const VectorXs &m,
    const std::vector<bool> &fixed, const VectorXu &geometry_indices,
    const std::vector<std::unique_ptr<RigidBody2DGeometry>> &geometry,
    const std::vector<std::unique_ptr<RigidBody2DForce>> &forces,
    const std::vector<RigidBody2DStaticPlane> &planes,
    const std::vector<PlanarPortal> &planar_portals,
    const std::vector<RigidBody2DStaticDrum> &drums)
    : m_q(q), m_q0(q), m_v(v), m_M(generateM(m)), m_Minv(generateMinv(m)),
      m_fixed(fixed), m_always_fixed(fixed),
      m_geometry_indices(geometry_indices),
      m_geometry(Utilities::cloneVector(geometry)),
      m_forces(Utilities::cloneVector(forces)), m_planes(planes),
      m_planar_portals(planar_portals), m_drums(drums), m_next_unique_id(),
      m_uique_ids(), m_homog_factor(VectorXs::Ones(m_q.size() / 3)),
      m_distance_field(VectorXs::Zero(m_q.size() / 3)), m_M0(generateM(m)),
      m_mass_weights(VectorXs::Zero(m_q.size() / 3)) {

  m_next_unique_id = 0;
  m_uique_ids.resize(numBodies());
  for (int bidx = 0; bidx < int(m_uique_ids.size()); bidx++) {
    m_uique_ids[bidx] = m_next_unique_id;
    m_next_unique_id++;
  }

#ifndef NDEBUG
  checkStateConsistency();
#endif

  //  std::cout << "DEM masses: " << std::endl << m_M << std::endl;
}

RigidBody2DState::RigidBody2DState(const RigidBody2DState &rhs)
    : m_q(rhs.m_q), m_q0(rhs.m_q0), m_v(rhs.m_v), m_M(rhs.m_M),
      m_Minv(rhs.m_Minv), m_fixed(rhs.m_fixed),
      m_always_fixed(rhs.m_always_fixed),
      m_geometry_indices(rhs.m_geometry_indices),
      m_geometry(Utilities::cloneVector(rhs.m_geometry)),
      m_forces(Utilities::cloneVector(rhs.m_forces)), m_planes(rhs.m_planes),
      m_planar_portals(rhs.m_planar_portals), m_drums(rhs.m_drums),
      m_next_unique_id(rhs.m_next_unique_id), m_uique_ids(rhs.m_uique_ids),
      m_homog_factor(rhs.m_homog_factor),
      m_distance_field(rhs.m_distance_field), m_M0(rhs.m_M0),
      m_mass_weights(rhs.m_mass_weights) {
#ifndef NDEBUG
  checkStateConsistency();
#endif
}

RigidBody2DState &RigidBody2DState::operator=(const RigidBody2DState &rhs) {
  RigidBody2DState copy{rhs};
  using std::swap;
  swap(*this, copy);
  return *this;
}

unsigned RigidBody2DState::nbodies() const { return unsigned(m_q.size()) / 3; }

unsigned RigidBody2DState::numBodies() const {
  return unsigned(m_q.size()) / 3;
}

unsigned RigidBody2DState::numSimulatedBodies() const {
  return static_cast<unsigned>(
      std::count(std::begin(m_fixed), std::end(m_fixed), false));
}

unsigned RigidBody2DState::numGeometryInstances() const {
  return unsigned(m_geometry.size());
}

scalar RigidBody2DState::computeTotalSimulatedMass() const {
  scalar M{0.0};
  const unsigned nbodies{numBodies()};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    if (!m_fixed[bdy_idx]) {
      assert(m_M.valuePtr()[3 * bdy_idx] == m_M.valuePtr()[3 * bdy_idx + 1]);
      M += m_M.valuePtr()[3 * bdy_idx];
    }
  }
  return M;
}

Vector2s RigidBody2DState::computeTotalSimulatedMomentum() const {
  Vector2s p{Vector2s::Zero()};
  const unsigned nbodies{numBodies()};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    if (!m_fixed[bdy_idx]) {
      assert(m_M.valuePtr()[3 * bdy_idx] == m_M.valuePtr()[3 * bdy_idx + 1]);
      p += m_M.valuePtr()[3 * bdy_idx] * m_v.segment<2>(3 * bdy_idx);
    }
  }
  return p;
}

VectorXs &RigidBody2DState::q() { return m_q; }

const VectorXs &RigidBody2DState::q() const { return m_q; }

VectorXs &RigidBody2DState::q0() { return m_q0; }

const VectorXs &RigidBody2DState::q0() const { return m_q0; }

VectorXs &RigidBody2DState::v() { return m_v; }

const VectorXs &RigidBody2DState::v() const { return m_v; }

const SparseMatrixsc &RigidBody2DState::M() const { return m_M; }

const SparseMatrixsc &RigidBody2DState::Minv() const { return m_Minv; }

const scalar &RigidBody2DState::getTotalMass(const unsigned bdy_idx) const {
  assert(bdy_idx < m_q.size() / 3);
  assert(m_M.valuePtr()[3 * bdy_idx] == m_M.valuePtr()[3 * bdy_idx + 1]);
  assert(m_M.valuePtr()[3 * bdy_idx] > 0.0);
  return m_M.valuePtr()[3 * bdy_idx];
}

const scalar &RigidBody2DState::m(const unsigned bdy_idx) const {
  assert(bdy_idx < m_q.size() / 3);
  assert(m_M.valuePtr()[3 * bdy_idx] == m_M.valuePtr()[3 * bdy_idx + 1]);
  assert(m_M.valuePtr()[3 * bdy_idx] > 0.0);
  return m_M.valuePtr()[3 * bdy_idx];
}

const scalar &RigidBody2DState::I(const unsigned bdy_idx) const {
  assert(bdy_idx < m_q.size() / 3);
  assert(m_M.valuePtr()[3 * bdy_idx + 2] > 0.0);
  return m_M.valuePtr()[3 * bdy_idx + 2];
}

bool RigidBody2DState::fixed(const int idx) const {
  assert(idx >= 0);
  assert(idx < static_cast<int>(m_fixed.size()));
  return m_fixed[idx];
}

void RigidBody2DState::addBody(const Vector2s &x, const scalar &theta,
                               const Vector2s &v, const scalar &omega,
                               const scalar &rho, const unsigned geo_idx,
                               const bool fixed, ImpactFrictionMap *ifmap) {
  assert(rho > 0.0);
  assert(geo_idx < m_geometry.size());

  assert(m_q.size() % 3 == 0);
  const unsigned original_num_bodies{unsigned(m_q.size() / 3)};
  const unsigned new_num_bodies{original_num_bodies + 1};

  // Format: x0, y0, theta0, x1, y1, theta1, ...
  m_q.conservativeResize(3 * new_num_bodies);
  m_q.segment<3>(3 * original_num_bodies) << x.x(), x.y(), theta;

  m_q0.conservativeResize(3 * new_num_bodies);
  m_q0.segment<3>(3 * original_num_bodies) << x.x(), x.y(), theta;

  // Format: vx0, vy0, omega0, vx1, vy1, omega1, ...
  m_v.conservativeResize(3 * new_num_bodies);
  m_v.segment<3>(3 * original_num_bodies) << v.x(), v.y(), omega;

  // Update the geometry references
  m_geometry_indices.conservativeResize(new_num_bodies);
  m_geometry_indices(original_num_bodies) = geo_idx;

  // Update fixed body tags
  m_fixed.push_back(fixed);
  m_always_fixed.push_back(false);

  m_homog_factor.conservativeResize(new_num_bodies);
  m_homog_factor(original_num_bodies) = 1.0;

  m_distance_field.conservativeResize(new_num_bodies);
  m_distance_field(original_num_bodies) = 0.0;

  m_mass_weights.conservativeResize(new_num_bodies);
  m_mass_weights(original_num_bodies) = 1.0;

  scalar m;
  scalar I;
  m_geometry[geo_idx]->computeMassAndInertia(rho, m, I);

  // Update the mass matrix
  {
    SparseMatrixsc M{SparseMatrixsc::Index(3 * new_num_bodies),
                     SparseMatrixsc::Index(3 * new_num_bodies)};
    M.reserve(3 * new_num_bodies);
    // Copy the old masses
    for (unsigned col = 0; col < 3 * original_num_bodies; ++col) {
      M.startVec(col);
      const unsigned row = col;
      M.insertBack(row, col) = m_M.valuePtr()[col];
    }
    // Insert the new masses
    M.startVec(3 * original_num_bodies);
    M.insertBack(3 * original_num_bodies, 3 * original_num_bodies) = m;
    M.startVec(3 * original_num_bodies + 1);
    M.insertBack(3 * original_num_bodies + 1, 3 * original_num_bodies + 1) = m;
    M.startVec(3 * original_num_bodies + 2);
    M.insertBack(3 * original_num_bodies + 2, 3 * original_num_bodies + 2) = I;
    M.finalize();
    m_M.swap(M);
    m_M0 = m_M;
  }
  // Update the inverse mass matrix
  {
    SparseMatrixsc Minv{SparseMatrixsc::Index(3 * new_num_bodies),
                        SparseMatrixsc::Index(3 * new_num_bodies)};
    Minv.reserve(3 * new_num_bodies);
    // Copy the old inverse masses
    for (unsigned col = 0; col < 3 * original_num_bodies; ++col) {
      Minv.startVec(col);
      const unsigned row = col;
      Minv.insertBack(row, col) = m_Minv.valuePtr()[col];
    }
    // Insert the new inverse masses
    Minv.startVec(3 * original_num_bodies);
    Minv.insertBack(3 * original_num_bodies, 3 * original_num_bodies) = 1.0 / m;
    Minv.startVec(3 * original_num_bodies + 1);
    Minv.insertBack(3 * original_num_bodies + 1, 3 * original_num_bodies + 1) =
        1.0 / m;
    Minv.startVec(3 * original_num_bodies + 2);
    Minv.insertBack(3 * original_num_bodies + 2, 3 * original_num_bodies + 2) =
        1.0 / I;
    Minv.finalize();
    m_Minv.swap(Minv);
  }

  m_uique_ids.emplace_back(m_next_unique_id);
  m_next_unique_id++;

#ifndef NDEBUG
  checkStateConsistency();
#endif

  if (ifmap != nullptr) {
    ifmap->enlargeCache(nbodies());
  }
}

void RigidBody2DState::addBody(const Vector2s &x, const Vector2s &x0,
                               const scalar &theta, const Vector2s &v,
                               const scalar &omega, const scalar &rho,
                               const unsigned geo_idx, const bool fixed,
                               ImpactFrictionMap *ifmap) {
  assert(rho > 0.0);
  assert(geo_idx < m_geometry.size());

  assert(m_q.size() % 3 == 0);
  const unsigned original_num_bodies{unsigned(m_q.size() / 3)};
  const unsigned new_num_bodies{original_num_bodies + 1};

  // Format: x0, y0, theta0, x1, y1, theta1, ...
  m_q.conservativeResize(3 * new_num_bodies);
  m_q.segment<3>(3 * original_num_bodies) << x.x(), x.y(), theta;

  m_q0.conservativeResize(3 * new_num_bodies);
  m_q0.segment<3>(3 * original_num_bodies) << x0.x(), x0.y(), theta;

  // Format: vx0, vy0, omega0, vx1, vy1, omega1, ...
  m_v.conservativeResize(3 * new_num_bodies);
  m_v.segment<3>(3 * original_num_bodies) << v.x(), v.y(), omega;

  // Update the geometry references
  m_geometry_indices.conservativeResize(new_num_bodies);
  m_geometry_indices(original_num_bodies) = geo_idx;

  // Update fixed body tags
  m_fixed.push_back(fixed);
  m_always_fixed.push_back(false);

  m_homog_factor.conservativeResize(new_num_bodies);
  m_homog_factor(original_num_bodies) = 1.0;

  m_distance_field.conservativeResize(new_num_bodies);
  m_distance_field(original_num_bodies) = 0.0;

  m_mass_weights.conservativeResize(new_num_bodies);
  m_mass_weights(original_num_bodies) = 1.0;

  scalar m;
  scalar I;
  m_geometry[geo_idx]->computeMassAndInertia(rho, m, I);

  // Update the mass matrix
  {
    SparseMatrixsc M{SparseMatrixsc::Index(3 * new_num_bodies),
                     SparseMatrixsc::Index(3 * new_num_bodies)};
    M.reserve(3 * new_num_bodies);
    // Copy the old masses
    for (unsigned col = 0; col < 3 * original_num_bodies; ++col) {
      M.startVec(col);
      const unsigned row = col;
      M.insertBack(row, col) = m_M.valuePtr()[col];
    }
    // Insert the new masses
    M.startVec(3 * original_num_bodies);
    M.insertBack(3 * original_num_bodies, 3 * original_num_bodies) = m;
    M.startVec(3 * original_num_bodies + 1);
    M.insertBack(3 * original_num_bodies + 1, 3 * original_num_bodies + 1) = m;
    M.startVec(3 * original_num_bodies + 2);
    M.insertBack(3 * original_num_bodies + 2, 3 * original_num_bodies + 2) = I;
    M.finalize();
    m_M.swap(M);
    m_M0 = m_M;
  }
  // Update the inverse mass matrix
  {
    SparseMatrixsc Minv{SparseMatrixsc::Index(3 * new_num_bodies),
                        SparseMatrixsc::Index(3 * new_num_bodies)};
    Minv.reserve(3 * new_num_bodies);
    // Copy the old inverse masses
    for (unsigned col = 0; col < 3 * original_num_bodies; ++col) {
      Minv.startVec(col);
      const unsigned row = col;
      Minv.insertBack(row, col) = m_Minv.valuePtr()[col];
    }
    // Insert the new inverse masses
    Minv.startVec(3 * original_num_bodies);
    Minv.insertBack(3 * original_num_bodies, 3 * original_num_bodies) = 1.0 / m;
    Minv.startVec(3 * original_num_bodies + 1);
    Minv.insertBack(3 * original_num_bodies + 1, 3 * original_num_bodies + 1) =
        1.0 / m;
    Minv.startVec(3 * original_num_bodies + 2);
    Minv.insertBack(3 * original_num_bodies + 2, 3 * original_num_bodies + 2) =
        1.0 / I;
    Minv.finalize();
    m_Minv.swap(Minv);
  }

  m_uique_ids.emplace_back(m_next_unique_id);
  m_next_unique_id++;

#ifndef NDEBUG
  checkStateConsistency();
#endif

  if (ifmap != nullptr) {
    ifmap->enlargeCache(nbodies());
  }
}

void RigidBody2DState::addBodies(
    const std::vector<Vector2s> &x, const std::vector<Vector2s> &x0,
    const std::vector<scalar> &theta, const std::vector<Vector2s> &v,
    const std::vector<scalar> &omega, const std::vector<scalar> &rho,
    const std::vector<unsigned> &geo_idx, const std::vector<bool> &fixed,
    ImpactFrictionMap *ifmap) {
  // assert( rho > 0.0 );
  // assert( geo_idx < m_geometry.size() );

  assert(m_q.size() % 3 == 0);
  const unsigned original_num_bodies{unsigned(m_q.size() / 3)};
  const unsigned new_num_bodies{original_num_bodies + unsigned(x.size())};

  // Format: x0, y0, theta0, x1, y1, theta1, ...
  m_q.conservativeResize(3 * new_num_bodies);
  m_q0.conservativeResize(3 * new_num_bodies);

  // Format: vx0, vy0, omega0, vx1, vy1, omega1, ...
  m_v.conservativeResize(3 * new_num_bodies);

  // Update the geometry references
  m_geometry_indices.conservativeResize(new_num_bodies);

  m_homog_factor.conservativeResize(new_num_bodies);

  m_distance_field.conservativeResize(new_num_bodies);

  m_mass_weights.conservativeResize(new_num_bodies);

  std::cout << "RigidBody2DState::addBodies(): " << std::endl;
  std::cout << "  original_num_bodies: " << original_num_bodies << std::endl;
  std::cout << "  new_num_bodies: " << new_num_bodies << std::endl;
  std::cout << "  m_fixed_size: " << m_fixed.size() << std::endl;

  for (size_t pnt_idx = 0; pnt_idx < x.size(); pnt_idx++) {
    m_q.segment<3>(3 * (original_num_bodies + pnt_idx)) << x[pnt_idx].x(),
        x[pnt_idx].y(), theta[pnt_idx];
    m_q0.segment<3>(3 * (original_num_bodies + pnt_idx)) << x0[pnt_idx].x(),
        x0[pnt_idx].y(), theta[pnt_idx];
    m_v.segment<3>(3 * (original_num_bodies + pnt_idx)) << v[pnt_idx].x(),
        v[pnt_idx].y(), omega[pnt_idx];
    m_geometry_indices((original_num_bodies + pnt_idx)) = geo_idx[pnt_idx];

    // Update fixed body tags
    m_fixed.push_back(fixed[pnt_idx]);
    m_always_fixed.push_back(false);

    m_homog_factor(original_num_bodies + pnt_idx) = 1.0;

    m_distance_field(original_num_bodies + pnt_idx) = 0.0;

    m_mass_weights(original_num_bodies + pnt_idx) = 1.0;
  }

  // Update the initial mass matrix
  SparseMatrixsc M0{SparseMatrixsc::Index(3 * new_num_bodies),
                    SparseMatrixsc::Index(3 * new_num_bodies)};
  M0.reserve(3 * new_num_bodies);
  {
    // Copy the old masses
    for (unsigned col = 0; col < 3 * original_num_bodies; ++col) {
      M0.startVec(col);
      const unsigned row = col;
      M0.insertBack(row, col) = m_M0.valuePtr()[col];
    }
  }

  for (size_t pnt_idx = 0; pnt_idx < x.size(); pnt_idx++) {
    scalar m;
    scalar I;
    m_geometry[geo_idx[pnt_idx]]->computeMassAndInertia(rho[pnt_idx], m, I);

    unsigned idx = original_num_bodies + pnt_idx;

    // Insert the new masses
    M0.startVec(3 * idx);
    M0.insertBack(3 * idx, 3 * idx) = m;
    M0.startVec(3 * idx + 1);
    M0.insertBack(3 * idx + 1, 3 * idx + 1) = m;
    M0.startVec(3 * idx + 2);
    M0.insertBack(3 * idx + 2, 3 * idx + 2) = I;
  }

  M0.finalize();
  m_M0.swap(M0);

  // Update the mass matrix
  SparseMatrixsc M{SparseMatrixsc::Index(3 * new_num_bodies),
                   SparseMatrixsc::Index(3 * new_num_bodies)};
  M.reserve(3 * new_num_bodies);
  {
    // Copy the old masses
    for (unsigned col = 0; col < 3 * original_num_bodies; ++col) {
      M.startVec(col);
      const unsigned row = col;
      M.insertBack(row, col) = m_M.valuePtr()[col];
    }
  }

  for (size_t pnt_idx = 0; pnt_idx < x.size(); pnt_idx++) {
    scalar m;
    scalar I;
    m_geometry[geo_idx[pnt_idx]]->computeMassAndInertia(rho[pnt_idx], m, I);

    unsigned idx = original_num_bodies + pnt_idx;

    // Insert the new masses
    M.startVec(3 * idx);
    M.insertBack(3 * idx, 3 * idx) = m;
    M.startVec(3 * idx + 1);
    M.insertBack(3 * idx + 1, 3 * idx + 1) = m;
    M.startVec(3 * idx + 2);
    M.insertBack(3 * idx + 2, 3 * idx + 2) = I;
  }

  M.finalize();
  m_M.swap(M);

  // Update the inverse mass matrix
  SparseMatrixsc Minv{SparseMatrixsc::Index(3 * new_num_bodies),
                      SparseMatrixsc::Index(3 * new_num_bodies)};
  Minv.reserve(3 * new_num_bodies);
  {
    // Copy the old inverse masses
    for (unsigned col = 0; col < 3 * original_num_bodies; ++col) {
      Minv.startVec(col);
      const unsigned row = col;
      Minv.insertBack(row, col) = m_Minv.valuePtr()[col];
    }
  }

  for (size_t pnt_idx = 0; pnt_idx < x.size(); pnt_idx++) {
    scalar m;
    scalar I;
    m_geometry[geo_idx[pnt_idx]]->computeMassAndInertia(rho[pnt_idx], m, I);

    unsigned idx = original_num_bodies + pnt_idx;

    // Insert the new inverse masses
    Minv.startVec(3 * idx);
    Minv.insertBack(3 * idx, 3 * idx) = 1.0 / m;
    Minv.startVec(3 * idx + 1);
    Minv.insertBack(3 * idx + 1, 3 * idx + 1) = 1.0 / m;
    Minv.startVec(3 * idx + 2);
    Minv.insertBack(3 * idx + 2, 3 * idx + 2) = 1.0 / I;
  }

  Minv.finalize();
  m_Minv.swap(Minv);

  for (int pnt_idx = 0; pnt_idx < int(x.size()); pnt_idx++) {
    m_uique_ids.emplace_back(m_next_unique_id);
    m_next_unique_id++;
  }

#ifndef NDEBUG
  checkStateConsistency();
#endif

  if (ifmap != nullptr) {
    ifmap->enlargeCache(nbodies());
  }
}

int RigidBody2DState::addCircleGeometry(const scalar &r) {
  int geo_idx = int(m_geometry.size());
  m_geometry.emplace_back(new CircleGeometry{r});
  return geo_idx;
}

void RigidBody2DState::removeBodiesIntersectingBoxes(
    const std::vector<std::pair<Array2s, Array2s>> &boxes,
    ImpactFrictionMap *ifmap) {
  if (boxes.size() > 1) {
    // Need a set to filter out case where body intersects multiple input boxes
    std::cerr << "Error, multiple boxes in removeBodiesIntersectingBoxes not "
                 "supported, yet."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::vector<unsigned> bodies_to_delete;
  const unsigned num_bodies = nbodies();

  for (const std::pair<Array2s, Array2s> &box : boxes) {
    const Vector2s x1 = 0.5 * (box.first + box.second);
    const scalar theta1 = 0.0;
    const Vector2s r1 = 0.5 * (box.second - box.first);

    for (unsigned idx = 0; idx < num_bodies; idx++) {
      const std::unique_ptr<RigidBody2DGeometry> &geo = bodyGeometry(idx);
      if (geo->type() != RigidBody2DGeometryType::CIRCLE) {
        continue;
      }
      const CircleGeometry &circle = *static_cast<CircleGeometry *>(geo.get());
      const scalar &r0 = circle.r();
      const Vector2s x0 = m_q.segment<2>(3 * idx);
      if (CircleBoxTools::isActive(x0, r0, x1, theta1, r1)) {
        bodies_to_delete.emplace_back(idx);
      }
    }
  }

  removeBodies(bodies_to_delete, ifmap);
}

void RigidBody2DState::removeBodies(const Eigen::Ref<const VectorXu> &indices,
                                    ImpactFrictionMap *ifmap) {
  if (indices.size() == 0) {
    return;
  }

  // Use q to mark dofs for deletion
  for (unsigned delete_idx = 0; delete_idx < indices.size(); ++delete_idx) {
    m_q(3 * indices[delete_idx]) = SCALAR_NAN;
  }

  if (ifmap != nullptr) {
    ifmap->deleteCacheEntries(m_q);
  }

  const unsigned nbodies_initial{static_cast<unsigned>(m_q.size()) / 3};

  Eigen::Map<VectorXs> M0_flat{m_M0.valuePtr(), m_M0.nonZeros()};
  Eigen::Map<VectorXs> M_flat{m_M.valuePtr(), m_M.nonZeros()};
  Eigen::Map<VectorXs> Minv_flat{m_Minv.valuePtr(), m_Minv.nonZeros()};

  unsigned copy_to = 0;
  unsigned copy_from = 0;
  for (; copy_from < nbodies_initial; copy_from++, copy_to++) {
    while (copy_from < nbodies_initial && std::isnan(m_q(3 * copy_from))) {
      copy_from++;
    }
    if (copy_from == nbodies_initial) {
      break;
    }
    m_q.segment<3>(3 * copy_to) = m_q.segment<3>(3 * copy_from);
    m_q0.segment<3>(3 * copy_to) = m_q0.segment<3>(3 * copy_from);
    m_v.segment<3>(3 * copy_to) = m_v.segment<3>(3 * copy_from);
    M0_flat.segment<3>(3 * copy_to) = M0_flat.segment<3>(3 * copy_from);
    M_flat.segment<3>(3 * copy_to) = M_flat.segment<3>(3 * copy_from);
    Minv_flat.segment<3>(3 * copy_to) = Minv_flat.segment<3>(3 * copy_from);
    m_fixed[copy_to] = m_fixed[copy_from];
    m_always_fixed[copy_to] = m_always_fixed[copy_from];
    m_geometry_indices(copy_to) = m_geometry_indices(copy_from);
    m_uique_ids[copy_to] = m_uique_ids[copy_from];
    m_homog_factor(copy_to) = m_homog_factor(copy_from);
    m_distance_field(copy_to) = m_distance_field(copy_from);
    m_mass_weights(copy_to) = m_mass_weights(copy_from);
  }

  const unsigned new_num_dofs{3 * copy_to};

  m_q.conservativeResize(new_num_dofs);
  m_q0.conservativeResize(new_num_dofs);
  m_v.conservativeResize(new_num_dofs);
  // m_M.conservativeResize( new_num_dofs, new_num_dofs );
  // m_M.makeCompressed();
  // m_Minv.conservativeResize( new_num_dofs, new_num_dofs );
  // m_Minv.makeCompressed();
  m_fixed.resize(copy_to);
  m_fixed.shrink_to_fit();
  m_always_fixed.resize(copy_to);
  m_always_fixed.shrink_to_fit();
  m_geometry_indices.conservativeResize(copy_to);
  m_uique_ids.resize(copy_to);
  m_homog_factor.conservativeResize(copy_to);
  m_distance_field.conservativeResize(copy_to);
  m_mass_weights.conservativeResize(copy_to);

  // Note: Conservative resize on sparse matrix seems to cause issues...
  // Update the initial mass matrix
  {
    SparseMatrixsc M0{SparseMatrixsc::Index(new_num_dofs),
                      SparseMatrixsc::Index(new_num_dofs)};
    M0.reserve(new_num_dofs);
    // Copy the updated masses over
    for (unsigned col = 0; col < new_num_dofs; ++col) {
      M0.startVec(col);
      M0.insertBack(col, col) = M0_flat(col);
    }
    M0.finalize();
    M0.makeCompressed();
    m_M0.swap(M0);
  }

  // Update the mass matrix
  {
    SparseMatrixsc M{SparseMatrixsc::Index(new_num_dofs),
                     SparseMatrixsc::Index(new_num_dofs)};
    M.reserve(new_num_dofs);
    // Copy the updated masses over
    for (unsigned col = 0; col < new_num_dofs; ++col) {
      M.startVec(col);
      M.insertBack(col, col) = M_flat(col);
    }
    M.finalize();
    M.makeCompressed();
    m_M.swap(M);
  }

  // Update the inverse mass matrix
  {
    SparseMatrixsc Minv{SparseMatrixsc::Index(new_num_dofs),
                        SparseMatrixsc::Index(new_num_dofs)};
    Minv.reserve(new_num_dofs);
    // Copy the old inverse masses
    for (unsigned col = 0; col < new_num_dofs; ++col) {
      Minv.startVec(col);
      Minv.insertBack(col, col) = Minv_flat(col);
    }
    Minv.finalize();
    Minv.makeCompressed();
    m_Minv.swap(Minv);
  }

#ifndef NDEBUG
  checkStateConsistency();
#endif
}

void RigidBody2DState::removeBodies(const std::vector<unsigned> &indices,
                                    ImpactFrictionMap *ifmap) {
  const Eigen::Map<const VectorXu> index_map{indices.data(),
                                             long(indices.size())};
  removeBodies(index_map, ifmap);
}

void RigidBody2DState::removeBodies(const std::vector<unsigned> &indices,
                                    std::vector<int> &old_to_new_idx_ref_map,
                                    ImpactFrictionMap *ifmap) {
  const unsigned nbodies_initial{static_cast<unsigned>(m_q.size()) / 3};
  const Eigen::Map<const VectorXu> index_map{indices.data(),
                                             long(indices.size())};
  removeBodies(index_map, ifmap);

  old_to_new_idx_ref_map.clear();
  unsigned pt = 0;
  unsigned dest = 0;
  for (unsigned idx = 0; idx < nbodies_initial; idx++) {
    bool erase = false;
    if (pt < indices.size()) {
      if (indices[pt] == idx) {
        erase = true;
        pt++;
      }
    }

    if (erase)
      old_to_new_idx_ref_map.push_back(-1);
    else {
      old_to_new_idx_ref_map.push_back(dest);
      dest++;
    }
  }
}

void RigidBody2DState::removeGeometry(
    const Eigen::Ref<const VectorXu> &indices) {
  if (indices.size() == 0) {
    return;
  }

  // Mark geometry slated for deletion
  for (unsigned delete_idx = 0; delete_idx < indices.size(); ++delete_idx) {
    m_geometry[indices[delete_idx]] = nullptr;
  }

  const unsigned ngeo_initial{static_cast<unsigned>(m_geometry.size())};

  // Map from old geo indices to new geo indices. count, ngeo
  VectorXu index_map(ngeo_initial);

  unsigned copy_to = 0;
  unsigned copy_from = 0;
  for (; copy_from < ngeo_initial; copy_from++, copy_to++) {
    while (copy_from < ngeo_initial && m_geometry[copy_from] == nullptr) {
#ifndef NDEBUG
      index_map(copy_from) = std::numeric_limits<unsigned>::max();
#endif
      copy_from++;
    }
    if (copy_from == ngeo_initial) {
      break;
    }
    m_geometry[copy_to] = std::move(m_geometry[copy_from]);
    index_map[copy_from] = copy_to;
  }

  m_geometry.resize(copy_to);
  m_geometry.shrink_to_fit();

  for (unsigned bdy_idx = 0; bdy_idx < m_geometry_indices.size(); ++bdy_idx) {
    m_geometry_indices(bdy_idx) = index_map(m_geometry_indices(bdy_idx));
  }

#ifndef NDEBUG
  checkStateConsistency();
#endif
}

void RigidBody2DState::removeGeometry(const std::vector<unsigned> &indices) {
  const Eigen::Map<const VectorXu> index_map{indices.data(),
                                             long(indices.size())};
  removeGeometry(index_map);
}

bool RigidBody2DState::alwaysFixed(const int idx) const {
  assert(idx >= 0);
  assert(idx < static_cast<int>(m_always_fixed.size()));
  return m_always_fixed[idx];
}

void RigidBody2DState::fixBody(const unsigned idx) {
  assert(idx < m_fixed.size());
  m_fixed[idx] = true;
}

void RigidBody2DState::unfixBody(const unsigned idx) {
  assert(idx < m_fixed.size());
  m_fixed[idx] = false;
}

std::vector<std::unique_ptr<RigidBody2DGeometry>> &
RigidBody2DState::geometry() {
  return m_geometry;
}

const std::vector<std::unique_ptr<RigidBody2DGeometry>> &
RigidBody2DState::geometry() const {
  return m_geometry;
}

VectorXu &RigidBody2DState::geometryIndices() { return m_geometry_indices; }

const VectorXu &RigidBody2DState::geometryIndices() const {
  return m_geometry_indices;
}

const std::unique_ptr<RigidBody2DGeometry> &
RigidBody2DState::bodyGeometry(const unsigned bdy_idx) const {
  assert(bdy_idx < m_geometry_indices.size());
  assert(m_geometry_indices(bdy_idx) < m_geometry.size());
  return m_geometry[m_geometry_indices(bdy_idx)];
}

const std::vector<std::unique_ptr<RigidBody2DForce>> &
RigidBody2DState::forces() const {
  return m_forces;
}

std::vector<RigidBody2DStaticPlane> &RigidBody2DState::planes() {
  return m_planes;
}

const std::vector<RigidBody2DStaticPlane> &RigidBody2DState::planes() const {
  return m_planes;
}

void RigidBody2DState::deleteStaticPlane(const unsigned plane_index) {
  assert(plane_index < m_planes.size());
  m_planes.erase(m_planes.begin() + plane_index);
}

std::vector<PlanarPortal> &RigidBody2DState::planarPortals() {
  return m_planar_portals;
}

const std::vector<PlanarPortal> &RigidBody2DState::planarPortals() const {
  return m_planar_portals;
}

const std::vector<RigidBody2DStaticDrum> &RigidBody2DState::drums() const {
  return m_drums;
}

std::vector<RigidBody2DStaticDrum> &RigidBody2DState::drums() {
  return m_drums;
}

Array4s RigidBody2DState::computeBoundingBox() const {
  const unsigned nbodies{static_cast<unsigned>(m_q.size() / 3)};
  assert(nbodies > 0);

  // For each body
  Array4s bounds;
  bounds.segment<2>(0).setConstant(SCALAR_INFINITY);
  bounds.segment<2>(2).setConstant(-SCALAR_INFINITY);
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    const Vector2s x{m_q.segment<2>(3 * bdy_idx)};
    const scalar theta{m_q(3 * bdy_idx + 2)};
    Array2s body_min;
    Array2s body_max;
    m_geometry[m_geometry_indices(bdy_idx)]->computeAABB(x, theta, body_min,
                                                         body_max);
    bounds.segment<2>(0) = bounds.segment<2>(0).min(body_min);
    bounds.segment<2>(2) = bounds.segment<2>(2).max(body_max);
  }
  assert((bounds.segment<2>(0) < bounds.segment<2>(2)).all());

  return bounds;
}

void RigidBody2DState::getAllCircleBodies(Matrix2Xsc &pos, VectorXs &radii,
                                          VectorXu &bdy_idx) const {
  const int nBodies = numBodies();
  pos.resize(2, nBodies);
  radii.resize(nBodies);
  bdy_idx.resize(nBodies);

  int np = 0;

  for (int idx = 0; idx < nBodies; idx++) {
    const std::unique_ptr<RigidBody2DGeometry> &current_geo{bodyGeometry(idx)};
    if (current_geo->type() == RigidBody2DGeometryType::CIRCLE) {
      const CircleGeometry *circle_geo =
          dynamic_cast<const CircleGeometry *>(current_geo.get());

      pos.col(np) = q().segment<2>(3 * idx);
      radii(np) = circle_geo->r();
      bdy_idx(np) = idx;
      np++;
    }
  }

  pos.conservativeResize(2, np);
  radii.conservativeResize(np);
  bdy_idx.conservativeResize(np);
}

void RigidBody2DState::getAllNonFixedCircleBodies(Matrix2Xsc &pos,
                                                  VectorXs &radii,
                                                  VectorXu &bdy_idx) const {
  const int nBodies = numBodies();
  pos.resize(2, nBodies);
  radii.resize(nBodies);
  bdy_idx.resize(nBodies);

  int np = 0;

  for (int idx = 0; idx < nBodies; idx++) {
    if (m_fixed[idx])
      continue;

    const std::unique_ptr<RigidBody2DGeometry> &current_geo{bodyGeometry(idx)};
    if (current_geo->type() == RigidBody2DGeometryType::CIRCLE) {
      const CircleGeometry *circle_geo =
          dynamic_cast<const CircleGeometry *>(current_geo.get());

      pos.col(np) = q().segment<2>(3 * idx);
      radii(np) = circle_geo->r();
      bdy_idx(np) = idx;
      np++;
    }
  }

  pos.conservativeResize(2, np);
  radii.conservativeResize(np);
  bdy_idx.conservativeResize(np);
}

void RigidBody2DState::serialize(std::ostream &output_stream) const {
  assert(output_stream.good());
  MathUtilities::serialize(m_q, output_stream);
  MathUtilities::serialize(m_q0, output_stream);
  MathUtilities::serialize(m_v, output_stream);
  MathUtilities::serialize(m_M, output_stream);
  MathUtilities::serialize(m_Minv, output_stream);
  Utilities::serializeVectorBuiltInType(m_fixed, output_stream);
  Utilities::serializeVectorBuiltInType(m_always_fixed, output_stream);
  MathUtilities::serialize(m_geometry_indices, output_stream);
  Utilities::serializeVectorCustomTypePointers(m_geometry, output_stream);
  Utilities::serializeVectorCustomTypePointers(m_forces, output_stream);
  Utilities::serializeVectorCustomType(m_planes, output_stream);
  Utilities::serializeVectorCustomType(m_planar_portals, output_stream);
  MathUtilities::serialize(m_homog_factor, output_stream);
  MathUtilities::serialize(m_distance_field, output_stream);
  MathUtilities::serialize(m_M0, output_stream);
  MathUtilities::serialize(m_mass_weights, output_stream);
}

static void
deserializeGeo(std::istream &input_stream,
               std::vector<std::unique_ptr<RigidBody2DGeometry>> &geo) {
  geo.clear();
  const std::vector<std::unique_ptr<RigidBody2DGeometry>>::size_type ngeo{
      Utilities::deserialize<
          std::vector<std::unique_ptr<RigidBody2DGeometry>>::size_type>(
          input_stream)};
  geo.resize(ngeo);
  for (std::vector<std::unique_ptr<RigidBody2DGeometry>>::size_type geo_idx = 0;
       geo_idx < ngeo; ++geo_idx) {
    // Read in the geometry type
    const RigidBody2DGeometryType geo_type{
        Utilities::deserialize<RigidBody2DGeometryType>(input_stream)};
    switch (geo_type) {
    case RigidBody2DGeometryType::CIRCLE: {
      geo[geo_idx].reset(new CircleGeometry{input_stream});
      break;
    }
    case RigidBody2DGeometryType::BOX: {
      geo[geo_idx].reset(new BoxGeometry{input_stream});
      break;
    }
    case RigidBody2DGeometryType::ANNULUS: {
      geo[geo_idx].reset(new AnnulusGeometry{input_stream});
      break;
    }
    }
  }
}

static void
deserializeForces(std::istream &input_stream,
                  std::vector<std::unique_ptr<RigidBody2DForce>> &forces) {
  forces.clear();
  const std::vector<std::unique_ptr<RigidBody2DForce>>::size_type nforces{
      Utilities::deserialize<
          std::vector<std::unique_ptr<RigidBody2DForce>>::size_type>(
          input_stream)};
  forces.resize(nforces);
  for (std::vector<std::unique_ptr<RigidBody2DForce>>::size_type force_idx = 0;
       force_idx < nforces; ++force_idx) {
    // Read in the force name
    const std::string force_name{
        StringUtilities::deserializeString(input_stream)};
    if ("near_earth_gravity" == force_name) {
      forces[force_idx].reset(new NearEarthGravityForce{input_stream});
    }
    // TODO: Check for hybrid drag force here
    else {
      std::cerr << "Invalid 2D force type " << force_name
                << " encountered in deserializeForces, this is a bug, exiting."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
}

void RigidBody2DState::deserialize(std::istream &input_stream) {
  assert(input_stream.good());

  m_q = MathUtilities::deserialize<VectorXs>(input_stream);
  m_q0 = MathUtilities::deserialize<VectorXs>(input_stream);
  m_v = MathUtilities::deserialize<VectorXs>(input_stream);
  MathUtilities::deserialize(m_M, input_stream);
  MathUtilities::deserialize(m_Minv, input_stream);
  Utilities::deserializeVectorBuiltInType(m_fixed, input_stream);
  Utilities::deserializeVectorBuiltInType(m_always_fixed, input_stream);
  m_geometry_indices = MathUtilities::deserialize<VectorXu>(input_stream);
  deserializeGeo(input_stream, m_geometry);
  deserializeForces(input_stream, m_forces);
  Utilities::deserializeVectorCustomType(m_planes, input_stream);
  Utilities::deserializeVectorCustomType(m_planar_portals, input_stream);
  m_homog_factor = MathUtilities::deserialize<VectorXs>(input_stream);
  m_distance_field = MathUtilities::deserialize<VectorXs>(input_stream);
  MathUtilities::deserialize(m_M0, input_stream);
  m_mass_weights = MathUtilities::deserialize<VectorXs>(input_stream);
#ifndef NDEBUG
  checkStateConsistency();
#endif
}

const scalar &RigidBody2DState::hybridFactor(const unsigned idx) const {
  assert(idx < m_homog_factor.size());
  return m_homog_factor(idx);
}

scalar &RigidBody2DState::hybridFactor(const unsigned idx) {
  assert(idx < m_homog_factor.size());
  return m_homog_factor(idx);
}

const VectorXs &RigidBody2DState::hybridFactors() const {
  return m_homog_factor;
}

const scalar &RigidBody2DState::massWeight(const unsigned idx) const {
  assert(idx < m_mass_weights.size());
  return m_mass_weights(idx);
}

scalar &RigidBody2DState::massWeight(const unsigned idx) {
  assert(idx < m_mass_weights.size());
  return m_mass_weights(idx);
}

const VectorXs &RigidBody2DState::massWeights() const { return m_mass_weights; }

scalar RigidBody2DState::totalMass() const {
  assert(m_q.size() % 3 == 0);
  assert(m_M.nonZeros() == m_q.size());

  const Eigen::Map<const VectorXs> M_flat{m_M.valuePtr(), m_M.nonZeros()};
  const unsigned nbodies{unsigned(m_q.size()) / 3};
  scalar total{0.0};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    assert(M_flat(3 * bdy_idx) == M_flat(3 * bdy_idx + 1));
    total += M_flat(3 * bdy_idx);
  }
  return total;
}

scalar RigidBody2DState::totalSimulatedMass() const {
  assert(m_q.size() % 3 == 0);
  assert(m_M.nonZeros() == m_q.size());

  const Eigen::Map<const VectorXs> M_flat{m_M.valuePtr(), m_M.nonZeros()};
  const unsigned nbodies{unsigned(m_q.size()) / 3};
  scalar total{0.0};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    assert(M_flat(3 * bdy_idx) == M_flat(3 * bdy_idx + 1));
    if (!fixed(bdy_idx)) {
      total += M_flat(3 * bdy_idx);
    }
  }
  return total;
}

const scalar &RigidBody2DState::distanceField(const unsigned idx) const {
  return m_distance_field(idx);
}

scalar &RigidBody2DState::distanceField(const unsigned idx) {
  return m_distance_field(idx);
}

void RigidBody2DState::modifyMass(const unsigned pnt_idx,
                                  const scalar &new_mass) {
  Eigen::Map<VectorXs> M_flat{m_M.valuePtr(), m_M.nonZeros()};
  const scalar percentage_change{new_mass / M_flat(3 * pnt_idx)};
  M_flat(3 * pnt_idx + 0) = new_mass;
  M_flat(3 * pnt_idx + 1) = new_mass;
  M_flat(3 * pnt_idx + 2) *= percentage_change;
  Eigen::Map<VectorXs> Minv_flat{m_Minv.valuePtr(), m_Minv.nonZeros()};
  Minv_flat(3 * pnt_idx + 0) = 1.0 / M_flat(3 * pnt_idx + 0);
  Minv_flat(3 * pnt_idx + 1) = 1.0 / M_flat(3 * pnt_idx + 1);
  Minv_flat(3 * pnt_idx + 2) = 1.0 / M_flat(3 * pnt_idx + 2);
}

void RigidBody2DState::setCurrentMassAsInitialMass() { m_M0 = m_M; }

void RigidBody2DState::updateCurrentMassUsingMassWeights() {
  const Eigen::Map<VectorXs> M0_flat{m_M0.valuePtr(), m_M0.nonZeros()};
  Eigen::Map<VectorXs> M_flat{m_M.valuePtr(), m_M.nonZeros()};
  Eigen::Map<VectorXs> Minv_flat{m_Minv.valuePtr(), m_Minv.nonZeros()};

  const unsigned nbodies{unsigned(m_q.size()) / 3};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    M_flat(3 * bdy_idx + 0) =
        M0_flat(3 * bdy_idx + 0) * m_mass_weights(bdy_idx);
    M_flat(3 * bdy_idx + 1) =
        M0_flat(3 * bdy_idx + 1) * m_mass_weights(bdy_idx);
    M_flat(3 * bdy_idx + 2) =
        M0_flat(3 * bdy_idx + 2) * m_mass_weights(bdy_idx);

    Minv_flat(3 * bdy_idx + 0) = 1.0 / M_flat(3 * bdy_idx + 0);
    Minv_flat(3 * bdy_idx + 1) = 1.0 / M_flat(3 * bdy_idx + 1);
    Minv_flat(3 * bdy_idx + 2) = 1.0 / M_flat(3 * bdy_idx + 2);
  }
}

void RigidBody2DState::updateMassWeights(
    std::function<scalar(const VectorXs &)> &weight_func) {
  const unsigned nbodies{unsigned(m_q.size()) / 3};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    const VectorXs q{m_q.segment<2>(3 * bdy_idx)};
    const scalar w = weight_func(q);
    m_mass_weights(bdy_idx) = w;
  }
}

scalar RigidBody2DState::computeMaximumSpeed() const {
  const unsigned nbodies{static_cast<unsigned>(m_q.size() / 3)};
  assert(nbodies > 0);
  scalar max_speed = 0.0;
  for (unsigned bdy_idx = 0; bdy_idx < nbodies; ++bdy_idx) {
    const Vector3s v{m_v.segment<3>(3 * bdy_idx)};
    max_speed = std::max<scalar>(max_speed, v.norm());
  }

  return max_speed;
}

int RigidBody2DState::getUniqueBodyIndex(const int idx) const {
  assert(idx >= 0);
  assert(idx < int(m_uique_ids.size()));
  return m_uique_ids[idx];
}

std::vector<int> &RigidBody2DState::uniqueBodyIDs() { return m_uique_ids; }

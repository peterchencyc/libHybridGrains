#include "ZoneTools.h"
#include "uniformgrid.h"
#include <iomanip>
#include <iostream>

void ZoneTools::init(
    const SimulationState &continuum_state, const scalar level_set_cell_width,
    bool allow_direct_transitions_between_discrete_and_continuum) {
  scalar cell_width_ratio =
      continuum_state.physics_grid.cell_width / level_set_cell_width;
  m_allow_direct_transitions_between_discrete_and_continuum =
      allow_direct_transitions_between_discrete_and_continuum;

  if (m_allow_direct_transitions_between_discrete_and_continuum)
    std::cout << "ALLOWING DIRECT TRANSITIONS BETWEEN DISCRETE AND CONTINUUM"
              << std::endl;
  else
    std::cout << "No Direct Transitions between discrete and continuum"
              << std::endl;

  if (cell_width_ratio < 1.0) {
    std::cerr << "ERROR: The level set cell width is larger than the mpm cell "
                 "width. This case is intentionally not supported. Exiting..."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (fabs(cell_width_ratio - floor(cell_width_ratio))) {
    std::cout << "WARNING: It seems the mpm cell width is not an integer "
                 "multiple of the level set cell width. Please consider to "
                 "choose a level set cell width such that the mpm cell width "
                 "is an integer multiple of it."
              << std::endl;
  }

  const Vector2u cell_count =
      (continuum_state.physics_grid.cell_count.cast<scalar>() *
       continuum_state.physics_grid.cell_width / level_set_cell_width)
          .unaryExpr([](const scalar &s) { return std::ceil(s); })
          .cast<unsigned>()
          .array();

  m_dist_grid.resize(continuum_state.physics_grid.min, cell_count,
                     level_set_cell_width);

  m_packing_fractions.setZero(m_dist_grid.num_cells());

  m_zone_indicators_level_set.resize(m_dist_grid.num_cells());
  m_zone_indicators_mpm.resize(continuum_state.physics_grid.numGridCells());
  m_prev_zone_indicators_mpm.resize(
      continuum_state.physics_grid.numGridCells());

  m_mpm_cell_count = continuum_state.physics_grid.cell_count;
  m_mpm_grid_min = continuum_state.physics_grid.min;
  m_mpm_cell_width = continuum_state.physics_grid.cell_width;

  for (int i = 0; i < int(continuum_state.physics_grid.numGridCells()); i++) {
    m_zone_indicators_mpm(i) = ZT_DISCRETE;
    m_prev_zone_indicators_mpm(i) = ZT_DISCRETE;
  }

  for (int i = 0; i < int(m_dist_grid.num_cells()); i++) {
    m_zone_indicators_level_set(i) = ZT_DISCRETE;
  }
}

bool ZoneTools::allowDirectTransitionsBetweenDiscreteAndContinuum() const {
  return m_allow_direct_transitions_between_discrete_and_continuum;
}

bool ZoneTools::isParticleInDiscreteZone_MPMIndicator(const Vector2s &x) const {
  Vector2i mpm_idx =
      ((x - m_mpm_grid_min) / m_mpm_cell_width)
          .unaryExpr([](const scalar &s) { return std::floor(s); })
          .cast<int>();

  if ((mpm_idx.array() < m_mpm_cell_count.cast<int>().array()).all() &&
      (mpm_idx.array() >= Vector2i::Zero().array()).all()) {
    const int flat_idx = mpm_idx.y() * m_mpm_cell_count.x() + mpm_idx.x();
    if (m_zone_indicators_mpm(flat_idx) == ZT_DISCRETE) {
      return true;
    }
  }

  return false;
}

bool ZoneTools::isParticleInContinuumZone_MPMIndicator(
    const Vector2s &x) const {
  Vector2i mpm_idx =
      ((x - m_mpm_grid_min) / m_mpm_cell_width)
          .unaryExpr([](const scalar &s) { return std::floor(s); })
          .cast<int>();

  if ((mpm_idx.array() < m_mpm_cell_count.cast<int>().array()).all() &&
      (mpm_idx.array() >= Vector2i::Zero().array()).all()) {
    const int flat_idx = mpm_idx.y() * m_mpm_cell_count.x() + mpm_idx.x();
    if (m_zone_indicators_mpm(flat_idx) == ZT_CONTINUUM) {
      return true;
    }
  }

  return false;
}

bool ZoneTools::isParticleInHybridZone_MPMIndicator(const Vector2s &x) const {
  Vector2i mpm_idx =
      ((x - m_mpm_grid_min) / m_mpm_cell_width)
          .unaryExpr([](const scalar &s) { return std::floor(s); })
          .cast<int>();

  if ((mpm_idx.array() < m_mpm_cell_count.cast<int>().array()).all() &&
      (mpm_idx.array() >= Vector2i::Zero().array()).all()) {
    const int flat_idx = mpm_idx.y() * m_mpm_cell_count.x() + mpm_idx.x();
    if (m_zone_indicators_mpm(flat_idx) == ZT_HYBRID) {
      return true;
    }
  }

  return false;
}

scalar ZoneTools::demWeight(const VectorXs &q) {
  if (isParticleInHybridZone_MPMIndicator(q.segment<2>(0))) {
    // std::cout << "H";
    return 0.5;
  } else if (isParticleInDiscreteZone_MPMIndicator(q.segment<2>(0))) {
    // std::cout << "D1";
    return 1.0;
  } else {
    // std::cout << "M0";
    return 0.0;
  }
}

scalar ZoneTools::mpmWeight(const VectorXs &q) {
  if (isParticleInHybridZone_MPMIndicator(q.segment<2>(0))) {
    // std::cout << "H";
    return 0.5;
  } else if (isParticleInContinuumZone_MPMIndicator(q.segment<2>(0))) {
    // std::cout << "M1";
    return 1.0;
  } else {
    // std::cout << "D0";
    return 0.0;
  }
}

scalar ZoneTools::computePackingFractionForBox(
    const Vector2s &box_min, const Vector2s &box_max, Matrix2Xsc &pos,
    VectorXs &r, CUniformGridH2D &ugrd, VectorXi &mail_box, const int test_id) {
  double intersection_volume_tot = 0.0;
  const Vector2s delta = box_max - box_min;
  double region_volume = delta.x() * delta.y();

  Vector2u min_grid_idx, max_grid_idx;
  ugrd.getGridIDRange(box_min, box_max, min_grid_idx, max_grid_idx);

  for (unsigned j = min_grid_idx.y(); j <= max_grid_idx.y(); j++) {
    for (unsigned i = min_grid_idx.x(); i <= max_grid_idx.x(); i++) {
      Vector2u cell_idx = Vector2u{i, j};
      int nData;
      const int *IDs;
      ugrd.getIDs(cell_idx, &nData, &IDs);

      for (int q = 0; q < nData; q++) {
        if (mail_box[IDs[q]] == test_id)
          continue;
        mail_box[IDs[q]] = test_id;

        Vector2s p = pos.col(IDs[q]);

        const double intersection_volume = circleRectangleIntersectionArea(
            box_min.x(), box_max.x(), box_min.y(), box_max.y(), p.x(), p.y(),
            r(IDs[q]));

        intersection_volume_tot += intersection_volume;
      }
    }
  }

  return intersection_volume_tot / region_volume;
}

void ZoneTools::updatePackingFraction(RigidBody2DState &discrete_state,
                                      SimulationState &continuum_state,
                                      const scalar phi_window_size,
                                      const int phi_samples_per_cell_side) {
  Matrix2Xsc pos;
  VectorXs r;
  VectorXu bdy_idx;
  discrete_state.getAllNonFixedCircleBodies(pos, r, bdy_idx);

  m_packing_fractions.setZero();

  VectorXi mailbox;
  mailbox.setZero(bdy_idx.size());

  const Vector2s region_min = m_dist_grid.min_coords;
  const Vector2u region_res = m_dist_grid.N;
  const scalar cell_width = m_dist_grid.h;

  CUniformGridH2D grd(region_min, region_res, cell_width);

  for (int i = 0; i < r.size(); i++) {

    Vector2s min_c = pos.col(i).array() - r(i);
    Vector2s max_c = pos.col(i).array() + r(i);

    grd.registerData(min_c, max_c, i);
  }

  int test_id = 1;

  for (int j = 0; j < int(region_res.y()); j++) {
    for (int i = 0; i < int(region_res.x()); i++) {
      scalar packing_fraction = 0.0;
      for (int v = 0; v < phi_samples_per_cell_side; v++) {
        for (int u = 0; u < phi_samples_per_cell_side; u++) {
          Vector2s _p{
              region_min.x() +
                  (i + (u + 0.5) / phi_samples_per_cell_side) * cell_width,
              region_min.y() +
                  (j + (v + 0.5) / phi_samples_per_cell_side) * cell_width};

          Vector2s min_c = _p.array() - phi_window_size * 0.5;
          Vector2s max_c = _p.array() + phi_window_size * 0.5;

          packing_fraction += computePackingFractionForBox(
              min_c, max_c, pos, r, grd, mailbox, test_id++);
        }
      }

      const int flat_idx = j * region_res.x() + i;

      m_packing_fractions(flat_idx) =
          packing_fraction /
          (phi_samples_per_cell_side * phi_samples_per_cell_side);
    }
  }

  // rasterize material points
  for (int p = 0; p < int(continuum_state.material_points.npoints); p++) {
    const Vector2s x = continuum_state.material_points.q.col(p) - region_min;
    const Vector2s x_min = x.array() - phi_window_size * 0.5;
    const Vector2s x_max = x.array() + phi_window_size * 0.5;

    const Vector2i gid_min =
        (x_min / cell_width)
            .unaryExpr([](const scalar &s) { return std::floor(s); })
            .cast<int>()
            .array();
    const Vector2i gid_max =
        (x_max / cell_width)
            .unaryExpr([](const scalar &s) { return std::floor(s); })
            .cast<int>()
            .array();

    for (int j = gid_min.y(); j <= gid_max.y(); j++) {
      for (int i = gid_min.x(); i <= gid_max.x(); i++) {
        const Vector2i gid = Vector2i{i, j};
        if ((gid.array() > 0).all() &&
            (gid.array() < region_res.cast<int>().array()).all()) {
          const int flat_idx = gid.y() * region_res.x() + gid.x();
          m_packing_fractions(flat_idx) = 1.0;
        }
      }
    }
  }
}

static void addSamplesAlongEdge(std::vector<Eigen::Vector2d> &samples,
                                const Eigen::Vector2d &p0,
                                const Eigen::Vector2d &p1,
                                const int n_samples_per_edge) {
  for (int i = 0; i < n_samples_per_edge; i++) {
    const double s = (i + 0.5) / n_samples_per_edge;
    const Eigen::Vector2d p = (1.0 - s) * p0 + s * p1;
    samples.push_back(p);
  }
}

void ZoneTools::updateDistanceField(scalar phi_threshold) {
  // Step1: generate boundary sampling points

  const int n_samples_per_edge = 4;
  std::vector<Eigen::Vector2d> boundary_samples;

  // int a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0, a6 = 0;

  for (int j = 0; j <= int(m_dist_grid.N.y()); j++) {
    for (int i = 0; i < int(m_dist_grid.N.x()); i++) {
      // Left side of the current cell
      const double cl = m_dist_grid.min_coords.x() + i * m_dist_grid.h;
      // Right side of the current cell
      const double cr = m_dist_grid.min_coords.x() + (i + 1) * m_dist_grid.h;
      // Top of the current cell
      const double ct = m_dist_grid.min_coords.y() + (j + 1) * m_dist_grid.h;

      const Eigen::Vector2d p_top_left(cl, ct);
      const Eigen::Vector2d p_top_right(cr, ct);

      const int cell_flat_idx = j * m_dist_grid.N.x() + i;
      const int cell_flat_idx_top = (j - 1) * m_dist_grid.N.x() + i;

      if (j == 0) {
        assert(cell_flat_idx < m_packing_fractions.size());
        if (m_packing_fractions(cell_flat_idx) > phi_threshold) {
          addSamplesAlongEdge(boundary_samples, p_top_left, p_top_right,
                              n_samples_per_edge);
          // a1++;
        }
      } else if (j == int(m_dist_grid.N.y())) {
        assert(cell_flat_idx_top < m_packing_fractions.size());
        if (m_packing_fractions(cell_flat_idx_top) > phi_threshold) {
          addSamplesAlongEdge(boundary_samples, p_top_left, p_top_right,
                              n_samples_per_edge);
          // a2++;
        }
      } else {
        if (cell_flat_idx_top >= m_packing_fractions.size()) {
          std::cerr << "Bad flat_idx_top!" << std::endl;
          std::cerr << cell_flat_idx_top << " of " << m_packing_fractions.size()
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (cell_flat_idx >= m_packing_fractions.size()) {
          std::cerr << "Bad flat_idx!" << std::endl;
          std::cerr << cell_flat_idx << " of " << m_packing_fractions.size()
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }

        assert(cell_flat_idx_top < m_packing_fractions.size());
        assert(cell_flat_idx < m_packing_fractions.size());
        if (((m_packing_fractions(cell_flat_idx_top) > phi_threshold) &&
             (m_packing_fractions(cell_flat_idx) <= phi_threshold)) ||
            ((m_packing_fractions(cell_flat_idx_top) <= phi_threshold) &&
             (m_packing_fractions(cell_flat_idx) > phi_threshold))) {
          addSamplesAlongEdge(boundary_samples, p_top_left, p_top_right,
                              n_samples_per_edge);
          // a3++;
        }
      }
    }
  }

  for (int j = 0; j < int(m_dist_grid.N.y()); j++) {
    for (int i = 0; i <= int(m_dist_grid.N.x()); i++) {
      Eigen::Vector2d p_bottom_left, p_top_left;
      const double cl = m_dist_grid.min_coords.x() + i * m_dist_grid.h;
      const double cb = m_dist_grid.min_coords.y() + j * m_dist_grid.h;
      const double ct = m_dist_grid.min_coords.y() + (j + 1) * m_dist_grid.h;

      p_bottom_left.x() = p_top_left.x() = cl;
      p_bottom_left.y() = cb;
      p_top_left.y() = ct;

      const int cell_flat_idx = j * (m_dist_grid.N.x()) + i;
      const int cell_flat_idx_left = j * (m_dist_grid.N.x()) + (i - 1);

      if (i == 0) {
        if (m_packing_fractions(cell_flat_idx) > phi_threshold) {
          addSamplesAlongEdge(boundary_samples, p_bottom_left, p_top_left,
                              n_samples_per_edge);
          // a4++;
        }
      } else if (i == int(m_dist_grid.N.x())) {
        if (m_packing_fractions(cell_flat_idx_left) > phi_threshold) {
          addSamplesAlongEdge(boundary_samples, p_bottom_left, p_top_left,
                              n_samples_per_edge);
          // a5++;
        }
      } else {
        if (((m_packing_fractions(cell_flat_idx_left) > phi_threshold) &&
             (m_packing_fractions(cell_flat_idx) <= phi_threshold)) ||
            ((m_packing_fractions(cell_flat_idx_left) <= phi_threshold) &&
             (m_packing_fractions(cell_flat_idx) > phi_threshold))) {
          addSamplesAlongEdge(boundary_samples, p_bottom_left, p_top_left,
                              n_samples_per_edge);
          // a6++;
        }
      }
    }
  }

  // std::cout << "a1 = " << a1 << ", a2 = " << a2 << ", a3 = " << a3
  //  << ", a4 = " << a4 << ", a5 = " << a5 << ", a6 = " << a6 << std::endl;

  m_boundary_samples.nPoints = unsigned(boundary_samples.size());
  m_boundary_samples.nAlloc = unsigned(boundary_samples.size());
  m_boundary_samples.x.setZero(2, boundary_samples.size());

  for (std::vector<Eigen::Vector2d>::size_type i = 0;
       i < boundary_samples.size(); i++) {
    m_boundary_samples.x.col(i) = boundary_samples[i];
  }

  // generate the distance grid:

  redistanceGrid(m_dist_grid, m_boundary_samples, m_dist_grid.h);

  // maybe also assign the signs:
  //   positive if (packing_fraction > T_phi), negative otherwise

  Eigen::VectorXi sign_flag;
  sign_flag.setZero(m_dist_grid.num_points());

  for (int j = 0; j < int(m_dist_grid.N.y()); j++) {
    for (int i = 0; i < int(m_dist_grid.N.x()); i++) {
      const int flat_idx_bl = j * (m_dist_grid.N.x() + 1) + i;
      const int flat_idx_br = j * (m_dist_grid.N.x() + 1) + i + 1;
      const int flat_idx_tl = (j + 1) * (m_dist_grid.N.x() + 1) + i;
      const int flat_idx_tr = (j + 1) * (m_dist_grid.N.x() + 1) + i + 1;

      const int cell_flat_idx_phi = j * (m_dist_grid.N.x()) + i;
      if (m_packing_fractions(cell_flat_idx_phi) > phi_threshold) {
        // positive sign
        m_dist_grid.minDist(flat_idx_bl) =
            fabs(m_dist_grid.minDist(flat_idx_bl));
        m_dist_grid.minDist(flat_idx_br) =
            fabs(m_dist_grid.minDist(flat_idx_br));
        m_dist_grid.minDist(flat_idx_tl) =
            fabs(m_dist_grid.minDist(flat_idx_tl));
        m_dist_grid.minDist(flat_idx_tr) =
            fabs(m_dist_grid.minDist(flat_idx_tr));

        sign_flag(flat_idx_bl) = 1;
        sign_flag(flat_idx_br) = 1;
        sign_flag(flat_idx_tl) = 1;
        sign_flag(flat_idx_tr) = 1;
      }
    }
  }

  for (int j = 0; j <= int(m_dist_grid.N.y()); j++) {
    for (int i = 0; i <= int(m_dist_grid.N.x()); i++) {
      const int flat_idx = j * (m_dist_grid.N.x() + 1) + i;
      if (sign_flag(flat_idx) == 0) {
        m_dist_grid.minDist(flat_idx) = -fabs(m_dist_grid.minDist(flat_idx));
      }
    }
  }
}

void ZoneTools::estimateLevelSetLevelZoneIndicators(
    scalar rzone_level_set, scalar rzone_half_thickness) {
  m_RZones_level_set_res.clear();

  const double interface_level_set = rzone_level_set;
  const double half_thickness = rzone_half_thickness;

  for (int j = 0; j < int(m_dist_grid.N.y()); j++) {
    for (int i = 0; i < int(m_dist_grid.N.x()); i++) {
      const int flat_idx_fl = j * (m_dist_grid.N.x() + 1) + i;
      const int flat_idx_fr = j * (m_dist_grid.N.x() + 1) + i + 1;
      const int flat_idx_bl = (j + 1) * (m_dist_grid.N.x() + 1) + i;
      const int flat_idx_br = (j + 1) * (m_dist_grid.N.x() + 1) + i + 1;

      const int cell_flat_idx_zone = j * (m_dist_grid.N.x()) + i;

      const double min_dist =
          std::min<double>(std::min<double>(m_dist_grid.minDist(flat_idx_fl),
                                            m_dist_grid.minDist(flat_idx_fr)),
                           std::min<double>(m_dist_grid.minDist(flat_idx_bl),
                                            m_dist_grid.minDist(flat_idx_br)));

      const double max_dist =
          std::max<double>(std::max<double>(m_dist_grid.minDist(flat_idx_fl),
                                            m_dist_grid.minDist(flat_idx_fr)),
                           std::max<double>(m_dist_grid.minDist(flat_idx_bl),
                                            m_dist_grid.minDist(flat_idx_br)));

      if (max_dist <= interface_level_set - half_thickness) {
        m_zone_indicators_level_set(cell_flat_idx_zone) = ZT_DISCRETE;
      } else if (min_dist >= interface_level_set + half_thickness) {
        m_zone_indicators_level_set(cell_flat_idx_zone) = ZT_CONTINUUM;
      } else {
        m_zone_indicators_level_set(cell_flat_idx_zone) = ZT_HYBRID;

        const Vector2s delta = (m_dist_grid.max_coords - m_dist_grid.min_coords)
                                   .cwiseQuotient(m_dist_grid.N.cast<scalar>());
        Vector2s cell_min = m_dist_grid.min_coords +
                            delta.cwiseProduct(Vector2s{scalar(i), scalar(j)});
        Vector2s cell_max = m_dist_grid.min_coords +
                            delta.cwiseProduct(Vector2s{i + 1.0, j + 1.0});

        Vector4s zone;
        zone << cell_min.x(), cell_min.y(), cell_max.x(), cell_max.y();
        m_RZones_level_set_res.push_back(zone);
      }
    }
  }
}

static void checkIndicatorConsistency(const VectorXi &indicator,
                                      Array2u cell_count) {
  int discont = 0;

  for (int j = 0; j < int(cell_count.y() - 1); j++) {
    for (int i = 0; i < int(cell_count.x() - 1); i++) {
      int type_c, type_x, type_y;
      const int c_flat_idx = j * cell_count.x() + i;
      const int x_flat_idx = j * cell_count.x() + i + 1;
      const int y_flat_idx = (j + 1) * cell_count.x() + i;

      type_c = indicator(c_flat_idx);
      type_x = indicator(x_flat_idx);
      type_y = indicator(y_flat_idx);

      if (type_c == ZT_CONTINUUM) {
        if (type_x == ZT_DISCRETE)
          discont++;
        if (type_y == ZT_DISCRETE)
          discont++;
      } else if (type_c == ZT_DISCRETE) {
        if (type_x == ZT_CONTINUUM)
          discont++;
        if (type_y == ZT_CONTINUUM)
          discont++;
      } else if (type_c == ZT_HYBRID) {
      }
    }
  }

  std::cout << "checkIndicatorConsistency(): " << discont << std::endl;
}

void ZoneTools::sampleMPMLevelZoneIndicators() {
  for (int i = 0; i < int(m_mpm_cell_count.x() * m_mpm_cell_count.y()); i++) {
    m_prev_zone_indicators_mpm(i) = m_zone_indicators_mpm(i);
  }

  m_RZones.clear();
  m_RZoneIndices.clear();
  m_DirectTransitionZonesFromDiscreteToContinuum.clear();
  m_DirectTransitionZonesFromDiscreteToContinuumIndices.clear();

  int nHybrid_zones_mpm = 0;
  int nContinuum_zones_mpm = 0;
  int nDiscrete_zones_mpm = 0;

  for (int j = 0; j < int(m_mpm_cell_count.y()); j++) {
    for (int i = 0; i < int(m_mpm_cell_count.x()); i++) {
      Vector2s c = Vector2s{scalar(i), scalar(j)};

      Vector2s cell_min = m_mpm_grid_min + c * m_mpm_cell_width;
      Vector2s cell_max =
          m_mpm_grid_min + Vector2s{c.array() + 1.0} * m_mpm_cell_width;

      Vector2i level_set_gid_min =
          ((cell_min - m_dist_grid.min_coords) / m_dist_grid.h)
              .unaryExpr([](const scalar &s) { return std::floor(s); })
              .cast<int>()
              .array();
      Vector2i level_set_gid_max =
          ((cell_max - m_dist_grid.min_coords) / m_dist_grid.h)
              .unaryExpr([](const scalar &s) { return std::floor(s); })
              .cast<int>()
              .array();

      level_set_gid_min =
          Vector2i{std::min<int>(m_dist_grid.N.x() - 1,
                                 std::max<int>(0, level_set_gid_min.x())),
                   std::min<int>(m_dist_grid.N.y() - 1,
                                 std::max<int>(0, level_set_gid_min.y()))};

      level_set_gid_max =
          Vector2i{std::min<int>(m_dist_grid.N.x() - 1,
                                 std::max<int>(0, level_set_gid_max.x())),
                   std::min<int>(m_dist_grid.N.y() - 1,
                                 std::max<int>(0, level_set_gid_max.y()))};

      int nHybrid = 0;
      int nContinuum = 0;
      int nDiscrete = 0;

      for (int v = level_set_gid_min.y(); v <= level_set_gid_max.y(); v++) {
        for (int u = level_set_gid_min.x(); u <= level_set_gid_max.x(); u++) {
          const int level_set_cell_flat_idx = v * m_dist_grid.N.x() + u;
          if (m_zone_indicators_level_set(level_set_cell_flat_idx) ==
              ZT_CONTINUUM)
            nContinuum++;
          else if (m_zone_indicators_level_set(level_set_cell_flat_idx) ==
                   ZT_DISCRETE)
            nDiscrete++;
          else if (m_zone_indicators_level_set(level_set_cell_flat_idx) ==
                   ZT_HYBRID)
            nHybrid++;
        }
      }

      const int mpm_cell_flat_idx = j * m_mpm_cell_count.x() + i;

      if (nHybrid > 0) {
        m_zone_indicators_mpm(mpm_cell_flat_idx) = ZT_HYBRID;
        Vector4s zone;
        zone << cell_min.x(), cell_min.y(), cell_max.x(), cell_max.y();
        m_RZones.push_back(zone);
        m_RZoneIndices.push_back(Vector2u{unsigned(i), unsigned(j)});
        nHybrid_zones_mpm++;
      } else if ((nDiscrete > 0) && (nContinuum > 0)) {
        std::cout
            << mpm_cell_flat_idx
            << ": nHybrid = 0, but nDiscrete > 0 and nContinuum > 0: nHybrid="
            << nHybrid << ", nDiscrete=" << nDiscrete
            << ", nContinuum=" << nContinuum << std::endl;
        std::cout << "level_set_gid_min = " << level_set_gid_min.x() << ", "
                  << level_set_gid_min.y() << std::endl;
        std::cout << "level_set_gid_max = " << level_set_gid_max.x() << ", "
                  << level_set_gid_max.y() << std::endl;
        m_zone_indicators_mpm(mpm_cell_flat_idx) = ZT_HYBRID;
        Vector4s zone;
        zone << cell_min.x(), cell_min.y(), cell_max.x(), cell_max.y();
        m_RZones.push_back(zone);
        m_RZoneIndices.push_back(Vector2u{unsigned(i), unsigned(j)});
        nHybrid_zones_mpm++;
      } else if (nDiscrete > 0) {
        if (m_allow_direct_transitions_between_discrete_and_continuum ||
            (m_prev_zone_indicators_mpm(mpm_cell_flat_idx) == ZT_DISCRETE) ||
            (m_prev_zone_indicators_mpm(mpm_cell_flat_idx) == ZT_HYBRID)) {
          m_zone_indicators_mpm(mpm_cell_flat_idx) = ZT_DISCRETE;
          nDiscrete_zones_mpm++;
        } else {
          m_zone_indicators_mpm(mpm_cell_flat_idx) = ZT_HYBRID;
          Vector4s zone;
          zone << cell_min.x(), cell_min.y(), cell_max.x(), cell_max.y();
          m_RZones.push_back(zone);
          m_RZoneIndices.push_back(Vector2u{unsigned(i), unsigned(j)});
          nHybrid_zones_mpm++;
        }
      } else if (nContinuum > 0) {
        if (m_allow_direct_transitions_between_discrete_and_continuum) {
          m_zone_indicators_mpm(mpm_cell_flat_idx) = ZT_CONTINUUM;
          nContinuum_zones_mpm++;
          if (m_prev_zone_indicators_mpm(mpm_cell_flat_idx) == ZT_DISCRETE) {
            Vector4s zone;
            zone << cell_min.x(), cell_min.y(), cell_max.x(), cell_max.y();
            m_DirectTransitionZonesFromDiscreteToContinuum.push_back(zone);
            m_DirectTransitionZonesFromDiscreteToContinuumIndices.push_back(
                Vector2u{unsigned(i), unsigned(j)});
          }
        } else {
          if ((m_prev_zone_indicators_mpm(mpm_cell_flat_idx) == ZT_CONTINUUM) ||
              (m_prev_zone_indicators_mpm(mpm_cell_flat_idx) == ZT_HYBRID)) {
            m_zone_indicators_mpm(mpm_cell_flat_idx) = ZT_CONTINUUM;
            nContinuum_zones_mpm++;
          } else {
            m_zone_indicators_mpm(mpm_cell_flat_idx) = ZT_HYBRID;
            Vector4s zone;
            zone << cell_min.x(), cell_min.y(), cell_max.x(), cell_max.y();
            m_RZones.push_back(zone);
            m_RZoneIndices.push_back(Vector2u{unsigned(i), unsigned(j)});
            nHybrid_zones_mpm++;
          }
        }
      } else {
        std::cerr << "This case should not happen..." << std::endl;
      }
    }
  }

  // for debug:
  // count how many zones changed their indicators
  int count = 0;
  for (int i = 0; i < int(m_mpm_cell_count.x() * m_mpm_cell_count.y()); i++) {
    if (m_prev_zone_indicators_mpm(i) != m_zone_indicators_mpm(i)) {
      count++;
    }
  }
  std::cout << "#diff zone indicators = " << count << std::endl;

  std::cout << "#hybrid zones = " << nHybrid_zones_mpm << std::endl;
  std::cout << "#Discrete zones = " << nDiscrete_zones_mpm << std::endl;
  std::cout << "#Continuum zones = " << nContinuum_zones_mpm << std::endl;
}

void ZoneTools::computeRZones(scalar rzone_level_set,
                              scalar rzone_half_thickness) {
  estimateLevelSetLevelZoneIndicators(rzone_level_set, rzone_half_thickness);
  checkIndicatorConsistency(m_zone_indicators_level_set, m_dist_grid.N);
  sampleMPMLevelZoneIndicators();
  checkIndicatorConsistency(m_zone_indicators_mpm, m_mpm_cell_count);
}

std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &
ZoneTools::getRZones() {
  return m_RZones;
}

std::vector<Vector2u> &ZoneTools::getRZoneIndices() { return m_RZoneIndices; }

std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &
ZoneTools::getDirectTransitionZonesFromDiscreteToContinuum() {
  return m_DirectTransitionZonesFromDiscreteToContinuum;
}

std::vector<Vector2u> &
ZoneTools::getDirectTransitionZonesFromDiscreteToContinuumIndices() {
  return m_DirectTransitionZonesFromDiscreteToContinuumIndices;
}

void ZoneTools::determineMPMParticlesToDelete(
    SimulationState &continuum_state, std::vector<unsigned> &id_list_to_del) {
  id_list_to_del.clear();

  // Delete material points that are in the full discrete zone
  const unsigned npoints_init{
      unsigned(continuum_state.material_points.npoints)};

  for (unsigned pnt_idx = 0; pnt_idx < npoints_init; ++pnt_idx) {
    const Eigen::Vector2d p = continuum_state.material_points.q.col(pnt_idx);

    const int idx_i =
        int(floor((p.x() - m_mpm_grid_min.x()) / m_mpm_cell_width));
    const int idx_j =
        int(floor((p.y() - m_mpm_grid_min.y()) / m_mpm_cell_width));

    if ((idx_i < 0) || (idx_i >= int(m_mpm_cell_count.x())) || (idx_j < 0) ||
        (idx_j >= int(m_mpm_cell_count.y()))) {
      id_list_to_del.push_back(pnt_idx);
    } else {
      const int flat_idx = idx_j * m_mpm_cell_count.x() + idx_i;
      if (m_zone_indicators_mpm(flat_idx) == ZT_DISCRETE) {
        id_list_to_del.push_back(pnt_idx);
      }
    }
  }
}

void ZoneTools::determineDEMGrainsToDelete(
    RigidBody2DState &discrete_state, std::vector<unsigned> &id_list_to_del) {
  id_list_to_del.clear();

  const unsigned nbodies_init{discrete_state.numBodies()};
  for (unsigned bdy_idx = 0; bdy_idx < nbodies_init; bdy_idx++) {
    // do not delete any fixed body even if it is outside of the mpm grid.
    if (discrete_state.fixed(bdy_idx))
      continue;

    const Eigen::Vector2d p = discrete_state.q().segment<2>(3 * bdy_idx);

    const int idx_i =
        int(floor((p.x() - m_mpm_grid_min.x()) / m_mpm_cell_width));
    const int idx_j =
        int(floor((p.y() - m_mpm_grid_min.y()) / m_mpm_cell_width));

    if ((idx_i < 0) || (idx_i >= int(m_mpm_cell_count.x())) || (idx_j < 0) ||
        (idx_j >= int(m_mpm_cell_count.y()))) {
      id_list_to_del.push_back(bdy_idx);
    } else {
      const int flat_idx = idx_j * m_mpm_cell_count.x() + idx_i;
      if (m_zone_indicators_mpm(flat_idx) == ZT_CONTINUUM) {
        id_list_to_del.push_back(bdy_idx);
      }
    }
  }
}

void ZoneTools::identifyInnerContinuumRegion(
    std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
    std::vector<Vector2u> &zone_indices) {
  zones.clear();
  zone_indices.clear();

  for (int j = 0; j < int(m_mpm_cell_count.y()); j++) {
    for (int i = 0; i < int(m_mpm_cell_count.x()); i++) {
      const int flat_idx = j * m_mpm_cell_count.x() + i;
      if (m_zone_indicators_mpm(flat_idx) != ZT_CONTINUUM)
        continue;

      const Eigen::Vector2d cell_minc =
          Eigen::Vector2d(m_mpm_grid_min.x() + i * m_mpm_cell_width,
                          m_mpm_grid_min.y() + j * m_mpm_cell_width);
      const Eigen::Vector2d cell_maxc =
          Eigen::Vector2d(m_mpm_grid_min.x() + (i + 1) * m_mpm_cell_width,
                          m_mpm_grid_min.y() + (j + 1) * m_mpm_cell_width);

      const Vector4s zone =
          Vector4s(cell_minc.x(), cell_minc.y(), cell_maxc.x(), cell_maxc.y());
      zones.push_back(zone);
      zone_indices.push_back(Vector2u{i, j});
    }
  }
}

void ZoneTools::identifyInnerContinuumRegionExcludingDirectTransitionZones(
    std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &zones,
    std::vector<Vector2u> &zone_indices) {
  zones.clear();

  for (int j = 0; j < int(m_mpm_cell_count.y()); j++) {
    for (int i = 0; i < int(m_mpm_cell_count.x()); i++) {
      const int flat_idx = j * m_mpm_cell_count.x() + i;
      if (m_zone_indicators_mpm(flat_idx) != ZT_CONTINUUM)
        continue;
      if (m_prev_zone_indicators_mpm(flat_idx) == ZT_DISCRETE)
        continue;

      const Eigen::Vector2d cell_minc =
          Eigen::Vector2d(m_mpm_grid_min.x() + i * m_mpm_cell_width,
                          m_mpm_grid_min.y() + j * m_mpm_cell_width);
      const Eigen::Vector2d cell_maxc =
          Eigen::Vector2d(m_mpm_grid_min.x() + (i + 1) * m_mpm_cell_width,
                          m_mpm_grid_min.y() + (j + 1) * m_mpm_cell_width);

      const Vector4s zone =
          Vector4s(cell_minc.x(), cell_minc.y(), cell_maxc.x(), cell_maxc.y());
      zones.push_back(zone);
      zone_indices.push_back(Vector2u{i, j});
    }
  }
}

CUniformGridH2D *
ZoneTools::setupUniformGridElementsAsPoints(RigidBody2DState &discrete_state) {
  CUniformGridH2D *ret = new CUniformGridH2D(
      Vector2s{m_mpm_grid_min}, m_mpm_cell_count, m_mpm_cell_width);

  for (int i = 0; i < int(discrete_state.nbodies()); i++) {
    const Vector2s x = discrete_state.q().segment<2>(3 * i);
    ret->registerPointData(x, i);
  }

  return ret;
}

CUniformGridH2D *
ZoneTools::setupUniformGridElementsAsPoints(SimulationState &continuum_state) {
  CUniformGridH2D *ret = new CUniformGridH2D(
      Vector2s{m_mpm_grid_min}, m_mpm_cell_count, m_mpm_cell_width);

  for (int i = 0; i < continuum_state.material_points.q.cols(); i++) {
    const Vector2s x = continuum_state.material_points.q.col(i);
    ret->registerPointData(x, i);
  }

  return ret;
}

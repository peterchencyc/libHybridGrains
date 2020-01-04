#include "Grid.h"

#include "mpmgrains2d/MaterialPoints.h"

static scalar computeCellWidth(const Vector2s &min, const Vector2s &max,
                               const Vector2i &cell_counts) {
  const Vector2s h{(max - min).array() / cell_counts.cast<scalar>().array()};
  assert(fabs(h(0) - h(1)) <= 1.0e-9);
  return h(0);
}

Grid::Grid(const Vector2s &min, const Vector2s &max,
           const Vector2i &cell_counts, const unsigned val_length)
    : m_N(cell_counts), m_min(min), m_max(max),
      m_h(computeCellWidth(min, max, cell_counts)),
      m_vals(MatrixXXsc::Zero(val_length,
                              (cell_counts(0) + 1) * (cell_counts(1) + 1))),
      m_m(VectorXs::Zero(m_vals.cols())),
      m_cell_occupied(cell_counts(0) * cell_counts(1), false),
      m_val_length(val_length) {}

const VectorXs &Grid::m() const { return m_m; }

scalar &Grid::m(const int node_idx) {
  assert(node_idx >= 0);
  assert(node_idx < m_m.size());
  assert(m_m(node_idx) >= 0.0);
  return m_m(node_idx);
}

const scalar &Grid::m(const int node_idx) const {
  assert(node_idx >= 0);
  assert(node_idx < m_m.size());
  assert(m_m(node_idx) >= 0.0);
  return m_m(node_idx);
}

scalar Grid::totalMass() const { return m_m.sum(); }

MatrixXXsc &Grid::vals() { return m_vals; }

std::vector<bool> &Grid::cellOccupied() { return m_cell_occupied; }

const std::vector<bool> &Grid::cellOccupied() const { return m_cell_occupied; }

VectorXs Grid::interpolateValue(const std::unique_ptr<BasisFunctions> &bf,
                                const MaterialPoints &points,
                                const Vector2s &x) const {
  assert(bf != nullptr);

  // TODO: Which h to use here?
  assert((points.hl.array() == points.hl(0)).all());
  const scalar h1{points.hl(0)};

  const std::pair<Array2u, Array2u> stencil{
      bf->computeStencil(x, h1, start(), cellWidth())};

  VectorXs val{VectorXs::Zero(m_vals.rows())};
#ifndef NDEBUG
  scalar w_sum{0.0};
#endif
  for (unsigned y_idx = stencil.first.y(); y_idx <= stencil.second.y();
       y_idx++) {
    for (unsigned x_idx = stencil.first.x(); x_idx <= stencil.second.x();
         x_idx++) {
      assert(nodeIndicesValid(x_idx, y_idx));
      const scalar w{bf->weight(x, h1, {x_idx, y_idx}, start(), cellWidth())};
#ifndef NDEBUG
      w_sum += w;
#endif
      const int grid_index{nodeIndex(x_idx, y_idx)};
      val += w * m_vals.col(grid_index);
    }
  }
  assert(fabs(w_sum - 1.0) <= 1.0e-9);
  return val;
}

MatrixXXsc Grid::interpolateGradient(const std::unique_ptr<BasisFunctions> &bf,
                                     const MaterialPoints &points,
                                     const Vector2s &x) const {
  assert(bf != nullptr);

  // TODO: Which h to use here?
  assert((points.hl.array() == points.hl(0)).all());
  const scalar h1{points.hl(0)};

  const std::pair<Array2u, Array2u> stencil{
      bf->computeStencil(x, h1, start(), cellWidth())};

  MatrixXXsc val{MatrixXXsc::Zero(m_vals.rows(), 2)};
#ifndef NDEBUG
  scalar w_sum{0.0};
  Vector2s wgrad_sum{0.0, 0.0};
#endif
  for (unsigned y_idx = stencil.first.y(); y_idx <= stencil.second.y();
       y_idx++) {
    for (unsigned x_idx = stencil.first.x(); x_idx <= stencil.second.x();
         x_idx++) {
      assert(nodeIndicesValid(x_idx, y_idx));
      const Vector2s gradw{
          bf->weightGrad(x, h1, {x_idx, y_idx}, start(), cellWidth())};
#ifndef NDEBUG
      w_sum += bf->weight(x, h1, {x_idx, y_idx}, start(), cellWidth());
      wgrad_sum += gradw;
#endif
      const int grid_index{nodeIndex(x_idx, y_idx)};
      val += m_vals.col(grid_index) * gradw.transpose();
    }
  }
  assert(fabs(w_sum - 1.0) <= 1.0e-9);
  assert((wgrad_sum.array() <= 1.0e-9).all());
  return val;
}

const Vector2s &Grid::start() const { return m_min; }

const scalar &Grid::cellWidth() const { return m_h; }

Vector2i Grid::containingCell(const Vector2s &x) const {
  // TODO: Floor operation redundant? Cast to int == floor because the value
  // shoudl be positive.
  return ((x - m_min) / m_h)
      .unaryExpr([](const scalar &y) { return floor(y); })
      .cast<int>();
}

int Grid::nodeIndex(const int i_idx, const int j_idx) const {
  assert((i_idx >= 0) && (j_idx >= 0) && (i_idx <= m_N.x()) &&
         (j_idx <= m_N.y()));
  return j_idx * (m_N.x() + 1) + i_idx;
}

int Grid::cellIndex(const int i_idx, const int j_idx) const {
  assert((i_idx >= 0) && (j_idx >= 0) && (i_idx < m_N.x()) &&
         (j_idx < m_N.y()));
  return j_idx * m_N.x() + i_idx;
}

int Grid::cellIndex(const Vector2s &x) const {
  const Vector2i containing_cell{containingCell(x)};
  return cellIndex(containing_cell[0], containing_cell[1]);
}

unsigned Grid::numNodes() const { return (m_N.x() + 1) * (m_N.y() + 1); }

const Vector2i &Grid::cellCount() const { return m_N; }

void Grid::rasterizeMass(const MaterialPoints &points, const int npoints,
                         const std::unique_ptr<BasisFunctions> &basis_funcs) {
  m_m.setZero();

  // TODO: Which h to use here?
  assert((points.hl.array() == points.hl(0)).all());
  const scalar h1{points.hl(0)};

  // For each material point
  for (unsigned pnt_idx = 0; pnt_idx < unsigned(npoints); pnt_idx++) {
    const Vector2s px{points.q.col(pnt_idx)};

    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(px, h1, m_min, m_h)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        // Switching to the new interface, for the compatibility with uGIMP
        // bases. const scalar weight{ basis_funcs->weight( px, { x_idx, y_idx
        // }, *this ) };
        const scalar weight{
            basis_funcs->weight(px, h1, {x_idx, y_idx}, m_min, m_h)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const int flat_node_idx{nodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < int(numNodes()));
        m_m(flat_node_idx) += weight * points.m(pnt_idx);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }
}

void Grid::rasterizeMassWeightedValuesVector(
    const MaterialPoints &points, const int npoints,
    const std::unique_ptr<BasisFunctions> &basis_funcs,
    const VectorXs &values) {
  m_vals.setZero();

  // TODO: Which h to use here?
  assert((points.hl.array() == points.hl(0)).all());
  const scalar h1{points.hl(0)};

  // For each material point
  for (unsigned pnt_idx = 0; pnt_idx < unsigned(npoints); pnt_idx++) {
    const Vector2s px{points.q.col(pnt_idx)};

    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(px, h1, m_min, m_h)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        // Switching to the new interface, for the compatibility with uGIMP
        // bases. const scalar weight{ basis_funcs->weight( px, { x_idx, y_idx
        // }, *this ) };
        const scalar weight{
            basis_funcs->weight(px, h1, {x_idx, y_idx}, m_min, m_h)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const int flat_node_idx{nodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < int(numNodes()));
        m_vals.col(flat_node_idx) +=
            weight * points.m(pnt_idx) *
            values.segment(pnt_idx * m_val_length, m_val_length);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }
}

void Grid::rasterizeMassWeightedValuesMatrix(
    const MaterialPoints &points, const int npoints,
    const std::unique_ptr<BasisFunctions> &basis_funcs,
    const MatrixXXsc &values) {
  m_vals.setZero();

  // TODO: Which h to use here?
  assert((points.hl.array() == points.hl(0)).all());
  const scalar h1{points.hl(0)};

  // For each material point
  for (unsigned pnt_idx = 0; pnt_idx < unsigned(npoints); pnt_idx++) {
    const Vector2s px{points.q.col(pnt_idx)};

    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(px, h1, m_min, m_h)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        // Switching to the new interface, for the compatibility with uGIMP
        // bases. const scalar weight{ basis_funcs->weight( px, { x_idx, y_idx
        // }, *this ) };
        const scalar weight{
            basis_funcs->weight(px, h1, {x_idx, y_idx}, m_min, m_h)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const int flat_node_idx{nodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < int(numNodes()));
        m_vals.col(flat_node_idx) +=
            weight * points.m(pnt_idx) * values.col(pnt_idx);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }
}

void Grid::rasterizeMassWeightedValues(
    const MaterialPoints &points, const int npoints,
    const std::unique_ptr<BasisFunctions> &basis_funcs,
    const std::vector<Matrix22sc, Eigen::aligned_allocator<Matrix22sc>>
        &values) {
  if (values.size() == 0)
    return;
  int cols = values[0].cols();
  int rows = values[0].rows();
  assert(m_vals.rows() == cols * rows);

  // TODO: Which h to use here?
  assert((points.hl.array() == points.hl(0)).all());
  const scalar h1{points.hl(0)};

  m_vals.setZero();

  // For each material point
  for (unsigned pnt_idx = 0; pnt_idx < unsigned(npoints); pnt_idx++) {
    const Vector2s px{points.q.col(pnt_idx)};

    const std::pair<Array2u, Array2u> stencil{
        basis_funcs->computeStencil(px, h1, m_min, m_h)};
    assert((stencil.first < stencil.second).all());

#ifndef NDEBUG
    scalar weight_sum{0.0};
#endif
    for (unsigned y_idx = stencil.first(1); y_idx <= stencil.second(1);
         y_idx++) {
      for (unsigned x_idx = stencil.first(0); x_idx <= stencil.second(0);
           x_idx++) {
        assert(nodeIndicesValid(x_idx, y_idx));
        // Switching to the new interface, for the compatibility with uGIMP
        // bases. const scalar weight{ basis_funcs->weight( px, { x_idx, y_idx
        // }, *this ) };
        const scalar weight{
            basis_funcs->weight(px, h1, {x_idx, y_idx}, m_min, m_h)};
#ifndef NDEBUG
        weight_sum += weight;
#endif
        const int flat_node_idx{nodeIndex(x_idx, y_idx)};
        assert(flat_node_idx < int(numNodes()));

        const Eigen::Map<const MatrixXXsc> data(values[pnt_idx].data(),
                                                cols * rows, 1);

        m_vals.col(flat_node_idx) += weight * points.m(pnt_idx) * data.col(0);
      }
    }
    assert(fabs(weight_sum - 1.0) <= 1.0e-6);
  }
}

void Grid::normalizationByMass() {
  for (unsigned node_idx = 0; node_idx < numNodes(); node_idx++) {
    if (m_m(node_idx) > 0.0) {
      m_vals.col(node_idx) /= m_m(node_idx);
    } else {
      m_vals.col(node_idx).setZero();
    }
  }
}

void Grid::normalizationByMass(const VectorXs &default_values_for_zero_mass) {
  for (unsigned node_idx = 0; node_idx < numNodes(); node_idx++) {
    if (m_m(node_idx) > 0.0) {
      m_vals.col(node_idx) /= m_m(node_idx);
    } else {
      m_vals.col(node_idx) =
          default_values_for_zero_mass.segment(0, m_val_length);
    }
  }
}

#ifndef NDEBUG
bool Grid::nodeIndicesValid(const int xidx, const int yidx) const {
  if (xidx < 0) {
    return false;
  }
  if (xidx > m_N(0)) {
    return false;
  }
  if (yidx < 0) {
    return false;
  }
  if (yidx > m_N(1)) {
    return false;
  }
  return true;
}
#endif

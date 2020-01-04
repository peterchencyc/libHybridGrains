#include "uniformgrid.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
using namespace std;

CUniformGridH2D::CUniformGridH2D(const Vector2s &in_min_coords,
                                 const Vector2u &in_res,
                                 const scalar in_cell_width)
    : res(in_res), min_coords(in_min_coords), cell_width(in_cell_width) {
  const int nCells = res.x() * res.y();
  m_nData = (int *)malloc(sizeof(int) * nCells);
  for (int i = 0; i < nCells; i++)
    m_nData[i] = 0;

  m_nAlloc = (int *)malloc(sizeof(int) * nCells);
  for (int i = 0; i < nCells; i++)
    m_nAlloc[i] = 0;

  m_IDs = (int **)malloc(sizeof(int *) * nCells);
  for (int i = 0; i < nCells; i++)
    m_IDs[i] = nullptr;
}

CUniformGridH2D::~CUniformGridH2D() {
  const int nCells = res.x() * res.y();
  for (int i = 0; i < nCells; i++)
    free(m_IDs[i]);
  free(m_IDs);
  free(m_nData);
  free(m_nAlloc);
}

void CUniformGridH2D::clearData() {
  const int nCells = res.x() * res.y();
  for (int i = 0; i < nCells; i++)
    m_nData[i] = 0;
}

void CUniformGridH2D::registerData(const Vector2s &in_min_coords,
                                   const Vector2s &in_max_coords, int id) {
  Vector2u idm;
  getGridID(in_min_coords, idm);
  Vector2u idM;
  getGridID(in_max_coords, idM);

  for (unsigned j = idm.y(); j <= idM.y(); j++) {
    for (unsigned i = idm.x(); i <= idM.x(); i++) {
      const unsigned flat_idx = j * res.x() + i;
      if (m_nAlloc[flat_idx] <= m_nData[flat_idx]) {
        int newSize;
        if (m_nAlloc[flat_idx] < 10)
          newSize = 16;
        else if (m_nAlloc[flat_idx] < 2000)
          newSize = m_nAlloc[flat_idx] * 2;
        else
          newSize = m_nAlloc[flat_idx] + 1024;

        m_nAlloc[flat_idx] = newSize;
        m_IDs[flat_idx] =
            (int *)realloc(m_IDs[flat_idx], sizeof(int) * m_nAlloc[flat_idx]);
      }

      m_IDs[flat_idx][m_nData[flat_idx]] = id;
      m_nData[flat_idx]++;
    }
  }
}

void CUniformGridH2D::registerPointData(const Vector2s &coords, int id) {
  Vector2u idm;
  getGridID(coords, idm);

  const unsigned flat_idx = idm.y() * res.x() + idm.x();
  if (m_nAlloc[flat_idx] <= m_nData[flat_idx]) {
    int newSize;
    if (m_nAlloc[flat_idx] < 10)
      newSize = 16;
    else if (m_nAlloc[flat_idx] < 2000)
      newSize = m_nAlloc[flat_idx] * 2;
    else
      newSize = m_nAlloc[flat_idx] + 1024;

    m_nAlloc[flat_idx] = newSize;
    m_IDs[flat_idx] =
        (int *)realloc(m_IDs[flat_idx], sizeof(int) * m_nAlloc[flat_idx]);
  }

  m_IDs[flat_idx][m_nData[flat_idx]] = id;
  m_nData[flat_idx]++;
}

void CUniformGridH2D::getGridID(const Vector2s &coords,
                                Vector2u &grid_idx) const {
  const Vector2s _idx = (coords - min_coords) / cell_width;
  grid_idx(0) =
      std::max<int>(0, std::min<int>(res.x() - 1, int(std::floor(_idx(0)))));
  grid_idx(1) =
      std::max<int>(0, std::min<int>(res.y() - 1, int(std::floor(_idx(1)))));
}

void CUniformGridH2D::getGridIDRange(const Vector2s &in_min_coords,
                                     const Vector2s &in_max_coords,
                                     Vector2u &min_grid_idx,
                                     Vector2u &max_grid_idx) const {
  getGridID(in_min_coords, min_grid_idx);
  getGridID(in_max_coords, max_grid_idx);
}

void CUniformGridH2D::getIDs(const Vector2u &grid_idx, int *out_nData,
                             const int **out_IDs) const {
  if ((grid_idx(0) >= res.x()) || (grid_idx(1) >= res.y())) {
    *out_nData = 0;
    *out_IDs = NULL;
  } else {
    const unsigned flat_idx = grid_idx(1) * res.x() + grid_idx(0);
    *out_nData = m_nData[flat_idx];
    *out_IDs = m_IDs[flat_idx];
  }
}

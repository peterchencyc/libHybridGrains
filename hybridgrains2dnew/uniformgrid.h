#ifndef UNIFORM_GRID_H2D_H
#define UNIFORM_GRID_H2D_H

#include <scisim/Math/MathDefines.h>

class CUniformGridH2D {
  CUniformGridH2D();

public:
  CUniformGridH2D(const Vector2s &in_min_coords, const Vector2u &in_res,
                  const scalar in_cell_width);
  ~CUniformGridH2D();

  void clearData();
  void registerData(const Vector2s &in_min_coords,
                    const Vector2s &in_max_coords, int id);
  void registerPointData(const Vector2s &coords, int id);
  void getGridID(const Vector2s &coords, Vector2u &grid_idx) const;
  void getGridIDRange(const Vector2s &in_min_coords,
                      const Vector2s &in_max_coords, Vector2u &min_grid_idx,
                      Vector2u &max_grid_idx) const;
  void getIDs(const Vector2u &grid_idx, int *out_nData,
              const int **out_IDs) const;

protected:
  Vector2u res;
  Vector2s min_coords;
  scalar cell_width;

  int *m_nData;
  int *m_nAlloc;
  int **m_IDs;
};

#endif

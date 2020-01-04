//
//  AvoidAVoid.h
//
//
//  Created by Yonghao Yue on 4/18/16.
//
//

#ifndef AvoidAVoid_h
#define AvoidAVoid_h

#include "quadtree.h"
#include "uniformgrid.h"
#include <functional>
#include <memory>

struct SPriorityData final {
  scalar value;
  Vector2s min_coords;
  Vector2s max_coords;
};

struct SQuadTreeNodePriorityData final {
  int nElem;
  std::vector<SPriorityData> data;
};

class CAvoidAVoid final {
public:
  CAvoidAVoid();
  ~CAvoidAVoid() = default;
  CAvoidAVoid(const CAvoidAVoid &other);
  CAvoidAVoid &operator=(const CAvoidAVoid &other);
  CAvoidAVoid(CAvoidAVoid &&) = default;
  CAvoidAVoid &operator=(CAvoidAVoid &&) = default;

  void setCheckConflictAgainstRegionBoundary(const bool do_check);

  // void setTargetRegion(const Vector2s& region_min, const Vector2s&
  // region_max);
  void setTargetRegions(
      const std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &regions);
  void setOldBodies(const Matrix2Xsc &pos, const VectorXs &radii,
                    const scalar &mean_radius);
  void avoidAVoid(std::function<double()> sampler);
  void getBodies(Matrix2Xsc &pos, VectorXs &radii);
  const std::vector<int> &getRZoneIDForNewlyInsertedElems();

  void serialize(std::ostream &output_stream) const;
  void deserialize(std::istream &input_stream);

  int getNumBodies();
  void getPos(const int i, Vector2s &pos);
  void getRadius(const int i, double &r);
  void getReduceRate(double &rate);

  int getTestID() const;
  void incrementTestID();

  CUniformGridH2D *getUniformGrid();
  std::vector<int> *getMailBox();

protected:
  void continueFillCell();
  void tryFillCell();
  void updateQuadTreeNodesActivity();
  int buildPriorityFromQuadTree();
  void drawUniformSample(int id, Vector2s &out_Pos);
  bool conflictCheck(const Vector2s &pos, const double &radius);
  void resizePoints(int nPnts);

  void randomFillRegion(const Vector2s &region_min, const Vector2s &region_max,
                        int numTrial);

private:
  Matrix2Xsc m_pos;
  VectorXs m_radius;
  scalar m_mean_radius;
  std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> m_regions;
  // Vector2s m_region_min;
  // Vector2s m_region_max;

  Vector2s m_current_region_min;
  Vector2s m_current_region_max;

  int m_testID;

  int m_nBodies;
  std::function<double()> m_RadiusSampler;
  double m_ReduceRate;

  bool m_CheckConflictAgainstRegionBoundary;

  std::unique_ptr<CQuadTree> m_qTree;
  double m_TolQTree;
  std::unique_ptr<SQuadTreeNodePriorityData> m_Priority;

  // TODO: Replace this with a unique pointer so it gets freed!
  std::unique_ptr<CUniformGridH2D> m_UniformGrid;
  std::vector<int> m_MailBox;

  std::vector<int> m_RZoneIDForNewlyInsertedElems;
};

#endif /* AvoidAVoid_h */

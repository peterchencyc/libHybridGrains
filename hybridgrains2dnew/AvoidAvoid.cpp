//
//  AvoidAVoid.cpp
//
//
//  Created by Yonghao Yue on 4/18/16.
//
//

#include "AvoidAvoid.h"
#include "random.h"
#include <stdio.h>

static void buildPriorityCallBack(SQuadTreeNode *in_Node, SNodeRange &in_Range,
                                  void *io_Data);
static void updateActivityCallBack(SQuadTreeNode *in_Node, SNodeRange &in_Range,
                                   void *io_Data);

inline bool operator<(const SPriorityData &lhs, const SPriorityData &rhs) {
  return (lhs.value < rhs.value) ||
         ((lhs.value == rhs.value) &&
          (std::tie(lhs.min_coords(0), lhs.min_coords(1), lhs.max_coords(0),
                    lhs.max_coords(1)) <
           std::tie(rhs.min_coords(0), rhs.min_coords(1), rhs.max_coords(0),
                    rhs.max_coords(1))));
}

static double defaultRadiusSampler() { return 0.1; }

CAvoidAVoid::CAvoidAVoid()
    : m_pos(), m_radius(), m_mean_radius(1.0), m_regions(), m_testID(100),
      m_nBodies(0), m_RadiusSampler(&defaultRadiusSampler), m_ReduceRate(0.97),
      m_CheckConflictAgainstRegionBoundary(false), m_qTree(nullptr),
      m_TolQTree(1.0e-4), m_Priority(nullptr), m_UniformGrid(nullptr) {}

CAvoidAVoid::CAvoidAVoid(const CAvoidAVoid &other) {
  std::cerr << "Code up copy constructor for CAvoidAVoid and children."
            << std::endl;
  // std::exit( EXIT_FAILURE );
}

CAvoidAVoid &CAvoidAVoid::operator=(const CAvoidAVoid &other) {
  CAvoidAVoid copy{other};
  using std::swap;
  swap(*this, copy);
  return *this;
}

void CAvoidAVoid::setCheckConflictAgainstRegionBoundary(const bool do_check) {
  m_CheckConflictAgainstRegionBoundary = do_check;
}

void CAvoidAVoid::setTargetRegions(
    const std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>> &regions) {
  // m_regions.clear();
  // for(std::vector<Vector4s>::size_type i=0; i<regions.size(); i++)
  //  m_regions.push_back( regions[i] );
  m_regions = regions;
}

void CAvoidAVoid::setOldBodies(const Matrix2Xsc &pos, const VectorXs &radii,
                               const scalar &mean_radius) {
  m_pos = pos;
  m_radius = radii;
  m_mean_radius = mean_radius;
  m_nBodies = int(m_pos.cols());
}

void CAvoidAVoid::getBodies(Matrix2Xsc &pos, VectorXs &radii) {
  pos = m_pos;
  radii = m_radius;
}

void CAvoidAVoid::serialize(std::ostream &output_stream) const {
  std::cerr << "Code up CAvoidAVoid::serialize." << std::endl;
  // std::exit( EXIT_FAILURE );
}

void CAvoidAVoid::deserialize(std::istream &input_stream) {
  std::cerr << "Code up CAvoidAVoid::deserialize." << std::endl;
  // std::exit( EXIT_FAILURE );
}

int CAvoidAVoid::getNumBodies() { return m_nBodies; }

void CAvoidAVoid::getPos(const int i, Vector2s &pos) { pos = m_pos.col(i); }

void CAvoidAVoid::getRadius(const int i, double &r) { r = m_radius(i); }

void CAvoidAVoid::getReduceRate(double &rate) { rate = m_ReduceRate; }

int CAvoidAVoid::getTestID() const { return m_testID; }

void CAvoidAVoid::incrementTestID() {
  m_testID++;
  if (m_testID == 0)
    m_testID = 1;
}

CUniformGridH2D *CAvoidAVoid::getUniformGrid() { return m_UniformGrid.get(); }

std::vector<int> *CAvoidAVoid::getMailBox() { return &m_MailBox; }

void CAvoidAVoid::avoidAVoid(std::function<double()> sampler) {
  m_RZoneIDForNewlyInsertedElems.clear();

  // assert( !m_regions.empty() );
  if (m_regions.empty()) {
    return;
  }
  Vector2s min_region = m_regions[0].segment<2>(0);
  Vector2s max_region = m_regions[0].segment<2>(2);
  for (std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>>::size_type i =
           1;
       i < m_regions.size(); i++) {
    min_region = min_region.cwiseMin(m_regions[i].segment<2>(0));
    max_region = max_region.cwiseMax(m_regions[i].segment<2>(2));
  }

  const Vector2s grid_min(min_region.x() - 6.0 * m_mean_radius,
                          min_region.y() - 6.0 * m_mean_radius);
  const Vector2s grid_max(max_region.x() + 6.0 * m_mean_radius,
                          max_region.y() + 6.0 * m_mean_radius);
  const scalar cell_width = 4.0 * m_mean_radius;
  const Vector2u grid_res(
      unsigned(ceil((grid_max.x() - grid_min.x()) / cell_width)),
      unsigned(ceil((grid_max.y() - grid_min.y()) / cell_width)));

  m_UniformGrid.reset(new CUniformGridH2D(grid_min, grid_res, cell_width));

  m_MailBox.clear();

  for (int i = 0; i < m_nBodies; i++) {
    const Vector2s min_c = m_pos.col(i).array() - m_radius(i);
    const Vector2s max_c = m_pos.col(i).array() + m_radius(i);
    m_UniformGrid->registerData(min_c, max_c, i);
    m_MailBox.push_back(0);
  }

  m_RadiusSampler = sampler;

  int count = m_nBodies;

  for (std::vector<Vector4s, Eigen::aligned_allocator<Vector4s>>::size_type i =
           0;
       i < m_regions.size(); i++) {
    m_current_region_min = m_regions[i].segment<2>(0);
    m_current_region_max = m_regions[i].segment<2>(2);

    randomFillRegion(m_current_region_min, m_current_region_max, 16);

    int new_count = m_nBodies;
    for (int j = count; j < new_count; j++) {
      m_RZoneIDForNewlyInsertedElems.push_back(i);
    }
    count = new_count;
  }
}

void CAvoidAVoid::randomFillRegion(const Vector2s &region_min,
                                   const Vector2s &region_max, int numTrial) {
  const Vector2s diff = region_max - region_min;

  for (int i = 0; i < numTrial; i++) {
    Vector2s pos{static_cast<scalar>(randomMT()) * diff(0) + region_min(0),
                 static_cast<scalar>(randomMT()) * diff(1) + region_min(1)};

    const double radius = m_RadiusSampler();

    const bool conflict = conflictCheck(pos, radius);
    if (!conflict) {
      const int id = m_nBodies;
      resizePoints(m_nBodies + 1);
      m_nBodies++;
      m_pos.col(id) = pos;
      m_radius(id) = radius;

      const Vector2s min_c = pos.array() - radius;
      const Vector2s max_c = pos.array() + radius;
      m_UniformGrid->registerData(min_c, max_c, id);
      m_MailBox.push_back(0);
    }
  }
}

const std::vector<int> &CAvoidAVoid::getRZoneIDForNewlyInsertedElems() {
  return m_RZoneIDForNewlyInsertedElems;
}

void CAvoidAVoid::continueFillCell() {
  buildPriorityFromQuadTree();

  const int nSample =
      std::min<int>(65536, static_cast<int>(m_Priority->data.size()));
  for (int i = 0; i < nSample; i++) {
    Vector2s pos;
    drawUniformSample(i, pos);
    const double radius = m_RadiusSampler();

    SQuadTreeNode *ptr = m_qTree->locateNode(pos);
    if (ptr->active == AS_INACTIVE) {
      continue;
    }
    m_qTree->subdivNode(ptr);

    const bool conflict = conflictCheck(pos, radius);
    if (!conflict) {
      const int id = m_nBodies;
      resizePoints(m_nBodies + 1);
      m_nBodies++;
      m_pos.col(id) = pos;
      m_radius(id) = radius;

      const Vector2s min_c = pos.array() - radius;
      const Vector2s max_c = pos.array() + radius;
      m_UniformGrid->registerData(min_c, max_c, id);
      m_MailBox.push_back(0);
    }
  }

  updateQuadTreeNodesActivity();
  m_qTree->mergeInactiveNodes();
}

void CAvoidAVoid::tryFillCell() {
  m_qTree->resetTree();
  for (int i = 0; i < 6; i++) {
    continueFillCell();
    const int nNonInactiveNodes = countNonInactiveNodes(m_qTree.get());
    if (nNonInactiveNodes == 0) {
      break;
    }
  }
}

void CAvoidAVoid::updateQuadTreeNodesActivity() {
  m_qTree->processAllNonInactiveLeafNodes(updateActivityCallBack, this);
}

int CAvoidAVoid::buildPriorityFromQuadTree() {
  int nNonInactiveNodes = countNonInactiveNodes(m_qTree.get());
  m_Priority->nElem = nNonInactiveNodes;
  m_Priority->data.resize(nNonInactiveNodes);
  m_Priority->data.clear();
  m_qTree->processAllNonInactiveLeafNodes(
      buildPriorityCallBack, static_cast<void *>(m_Priority.get()));

  std::sort(m_Priority->data.begin(), m_Priority->data.end());

  return nNonInactiveNodes;
}

void CAvoidAVoid::drawUniformSample(int id, Vector2s &out_Pos) {
  const Vector2s min_C = m_Priority->data[id].min_coords;
  const Vector2s diff =
      m_Priority->data[id].max_coords - m_Priority->data[id].min_coords;

  out_Pos(0) = static_cast<scalar>(randomMT()) * diff(0) + min_C(0);
  out_Pos(1) = static_cast<scalar>(randomMT()) * diff(1) + min_C(1);
}

bool CAvoidAVoid::conflictCheck(const Vector2s &pos, const double &radius) {
  const int testID = getTestID();

  Vector2s min_c = pos.array() - radius;
  Vector2s max_c = pos.array() + radius;

  Vector2u idm, idM;
  m_UniformGrid->getGridIDRange(min_c, max_c, idm, idM);

  for (unsigned v = idm.y(); v <= idM.y(); v++) {
    for (unsigned u = idm.x(); u <= idM.x(); u++) {
      Vector2u grid_idx = Vector2u{u, v};
      int nData;
      const int *ids;
      m_UniformGrid->getIDs(grid_idx, &nData, &ids);
      for (int k = 0; k < nData; k++) {
        if (m_MailBox[ids[k]] == testID)
          continue;
        m_MailBox[ids[k]] = testID;

        const double dist = (m_pos.col(ids[k]) - pos).norm();
        if (dist < (m_radius(ids[k]) + radius) * m_ReduceRate) {
          incrementTestID();
          return true;
        }
      }
    }
  }

  incrementTestID();

  if (m_CheckConflictAgainstRegionBoundary) {
    if (pos(0) < m_current_region_min(0) + radius)
      return true;
    if (pos(0) > m_current_region_max(0) - radius)
      return true;
    if (pos(1) < m_current_region_min(1) + radius)
      return true;
    if (pos(1) > m_current_region_max(1) - radius)
      return true;
  }

  return false;
}

void CAvoidAVoid::resizePoints(int nPnts) {
  m_pos.conservativeResize(2, nPnts);
  m_radius.conservativeResize(nPnts);
}

void buildPriorityCallBack(SQuadTreeNode *in_Node, SNodeRange &in_Range,
                           void *io_Data) {
  SQuadTreeNodePriorityData *dat = (SQuadTreeNodePriorityData *)io_Data;
  scalar volume =
      (in_Range.maxX - in_Range.minX) * (in_Range.maxY - in_Range.minY);

  SPriorityData the_Data;
  the_Data.min_coords(0) = in_Range.minX;
  the_Data.max_coords(0) = in_Range.maxX;
  the_Data.min_coords(1) = in_Range.minY;
  the_Data.max_coords(1) = in_Range.maxY;
  the_Data.value = volume;
  dat->data.push_back(the_Data);
}

void updateActivityCallBack(SQuadTreeNode *in_Node, SNodeRange &in_Range,
                            void *io_Data) {
  CAvoidAVoid *aav = (CAvoidAVoid *)io_Data;
  const Vector2s min_coords(in_Range.minX, in_Range.minY);
  const Vector2s max_coords(in_Range.maxX, in_Range.maxY);
  const double half_diag = (max_coords - min_coords).norm() * 0.5;
  const Vector2s cR = (min_coords + max_coords) * 0.5;

  Vector2u idm, idM;
  aav->getUniformGrid()->getGridIDRange(min_coords, max_coords, idm, idM);

  const int testID = aav->getTestID();
  std::vector<int> *mailBoxPtr = aav->getMailBox();

  for (unsigned v = idm.y(); v <= idM.y(); v++) {
    for (unsigned u = idm.x(); u <= idM.x(); u++) {
      Vector2u grid_idx = Vector2u{u, v};
      int nData;
      const int *ids;
      aav->getUniformGrid()->getIDs(grid_idx, &nData, &ids);
      for (int k = 0; k < nData; k++) {
        if ((*mailBoxPtr)[ids[k]] == testID)
          continue;
        (*mailBoxPtr)[ids[k]] = testID;

        double ri;
        aav->getRadius(ids[k], ri);
        double reduceRate;
        aav->getReduceRate(reduceRate);
        const double r = ri * reduceRate;

        Vector2s pos_i;
        aav->getPos(ids[k], pos_i);
        const double d_C_CR = (pos_i - cR).norm();
        const double dm = std::max<double>(0.0, d_C_CR - (r + half_diag));
        const double dM = std::max<double>(0.0, d_C_CR - (r - half_diag));

        if (dm >= r)
          in_Node->active = AS_ACTIVE;
        else if (dM <= r)
          in_Node->active = AS_INACTIVE;
        else
          in_Node->active = AS_PARTIAL_ACTIVE;
      }
    }
  }

  aav->incrementTestID();
}

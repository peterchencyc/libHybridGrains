#ifndef QUADTREE_CALLBACK_H
#define QUADTREE_CALLBACK_H

#include "scisim/Math/MathDefines.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

enum EActiveState {
  AS_ACTIVE,
  AS_INACTIVE,
  AS_PARTIAL_ACTIVE,
};

struct SQuadTreeNode {
  SQuadTreeNode *child[4];
  EActiveState active;
};

struct SNodeRange {
  double minX;
  double maxX;
  double minY;
  double maxY;
};

struct SQuadTreeStack {
  int ptr;
  int stackSize;
  SQuadTreeNode **nodePtrs;
  double *minXs;
  double *maxXs;
  double *minYs;
  double *maxYs;
};

inline void initStack(SQuadTreeStack *io_Stack) { io_Stack->ptr = 0; }

inline void pushStack(SQuadTreeStack *io_Stack, SQuadTreeNode *in_Node,
                      double in_MinX, double in_MaxX, double in_MinY,
                      double in_MaxY) {
  io_Stack->nodePtrs[io_Stack->ptr] = in_Node;
  io_Stack->minXs[io_Stack->ptr] = in_MinX;
  io_Stack->minYs[io_Stack->ptr] = in_MinY;
  io_Stack->maxXs[io_Stack->ptr] = in_MaxX;
  io_Stack->maxYs[io_Stack->ptr] = in_MaxY;
  io_Stack->ptr++;
}

inline void popStack(SQuadTreeStack *io_Stack, SQuadTreeNode **out_Node,
                     double *out_MinX, double *out_MaxX, double *out_MinY,
                     double *out_MaxY) {
  io_Stack->ptr--;
  *out_Node = io_Stack->nodePtrs[io_Stack->ptr];
  *out_MinX = io_Stack->minXs[io_Stack->ptr];
  *out_MinY = io_Stack->minYs[io_Stack->ptr];
  *out_MaxX = io_Stack->maxXs[io_Stack->ptr];
  *out_MaxY = io_Stack->maxYs[io_Stack->ptr];
}

inline bool isStackEmpty(const SQuadTreeStack *in_Stack) {
  return in_Stack->ptr == 0;
}

typedef void (*quadNodeCallBackFunc)(SQuadTreeNode *, SNodeRange &, void *);

class CQuadTree {
  CQuadTree();

public:
  CQuadTree(const Vector2s &in_Min, const Vector2s &in_Max) {
    m_RootNode = (SQuadTreeNode *)malloc(sizeof(SQuadTreeNode));
    m_RootNode->child[0] = m_RootNode->child[1] = m_RootNode->child[2] =
        m_RootNode->child[3] = NULL;
    m_RootNode->active = AS_ACTIVE;

    m_MinCoords = in_Min;
    m_MaxCoords = in_Max;

    m_Stack.stackSize = 1024 * 256 * 4;
    m_Stack.nodePtrs =
        (SQuadTreeNode **)malloc(sizeof(SQuadTreeNode *) * m_Stack.stackSize);
    m_Stack.minXs = (double *)malloc(sizeof(double) * m_Stack.stackSize);
    m_Stack.minYs = (double *)malloc(sizeof(double) * m_Stack.stackSize);
    m_Stack.maxXs = (double *)malloc(sizeof(double) * m_Stack.stackSize);
    m_Stack.maxYs = (double *)malloc(sizeof(double) * m_Stack.stackSize);
    initStack(&m_Stack);
  }

  ~CQuadTree() {
    freeNode(m_RootNode);
    free(m_RootNode);
    m_RootNode = NULL;

    free(m_Stack.nodePtrs);
    m_Stack.nodePtrs = NULL;
    free(m_Stack.minXs);
    m_Stack.minXs = NULL;
    free(m_Stack.minYs);
    m_Stack.minYs = NULL;
    free(m_Stack.maxXs);
    m_Stack.maxXs = NULL;
    free(m_Stack.maxYs);
    m_Stack.maxYs = NULL;
    m_Stack.stackSize = 0;
    m_Stack.ptr = 0;
  }

  const SQuadTreeNode *getRoot() const { return m_RootNode; }

  void resetTree() {
    freeNode(m_RootNode);
    m_RootNode->active = AS_ACTIVE;
  }

  SQuadTreeNode *locateNode(const Vector2s &p) {
    SQuadTreeNode *the_Node = m_RootNode;
    double xmin = m_MinCoords(0);
    double xmax = m_MaxCoords(0);
    double ymin = m_MinCoords(1);
    double ymax = m_MaxCoords(1);

    while (1) {
      if (the_Node->child[0] == NULL)
        return the_Node;

      double cx = (xmin + xmax) * 0.5;
      double cy = (ymin + ymax) * 0.5;

      int b0 = p(0) < cx ? 0 : 1;
      int b1 = p(1) < cy ? 0 : 1;
      int cid = (b1 << 1) | b0;
      the_Node = the_Node->child[cid];

      if (b0 == 0)
        xmax = cx;
      else
        xmin = cx;
      if (b1 == 0)
        ymax = cy;
      else
        ymin = cy;
    }
  }

  void subdivNode(SQuadTreeNode *io_Node) {
    for (int j = 0; j < 4; j++) {
      io_Node->child[j] = (SQuadTreeNode *)malloc(sizeof(SQuadTreeNode));
      for (int i = 0; i < 4; i++)
        io_Node->child[j]->child[i] = NULL;
      io_Node->child[j]->active = AS_ACTIVE;
    }
  }

  void processAllNonInactiveLeafNodes(quadNodeCallBackFunc in_Func,
                                      void *io_UserData) {
    initStack(&m_Stack);
    pushStack(&m_Stack, m_RootNode, m_MinCoords(0), m_MaxCoords(0),
              m_MinCoords(1), m_MaxCoords(1));

    while (!isStackEmpty(&m_Stack)) {
      SQuadTreeNode *the_Node;
      double minX, maxX, minY, maxY;
      popStack(&m_Stack, &the_Node, &minX, &maxX, &minY, &maxY);
      if (the_Node->active == AS_INACTIVE)
        continue;

      if (the_Node->child[0] != NULL) {
        double cx = (minX + maxX) * 0.5;
        double cy = (minY + maxY) * 0.5;
        pushStack(&m_Stack, the_Node->child[0], minX, cx, minY, cy);
        pushStack(&m_Stack, the_Node->child[1], cx, maxX, minY, cy);
        pushStack(&m_Stack, the_Node->child[2], minX, cx, cy, maxY);
        pushStack(&m_Stack, the_Node->child[3], cx, maxX, cy, maxY);
      } else {
        SNodeRange the_Range;
        the_Range.minX = minX;
        the_Range.maxX = maxX;
        the_Range.minY = minY;
        the_Range.maxY = maxY;
        in_Func(the_Node, the_Range, io_UserData);
      }
    }
  }

  void processAllLeafNodes(quadNodeCallBackFunc in_Func, void *io_UserData) {
    initStack(&m_Stack);
    pushStack(&m_Stack, m_RootNode, m_MinCoords(0), m_MaxCoords(0),
              m_MinCoords(1), m_MaxCoords(1));

    while (!isStackEmpty(&m_Stack)) {
      SQuadTreeNode *the_Node;
      double minX, maxX, minY, maxY;
      popStack(&m_Stack, &the_Node, &minX, &maxX, &minY, &maxY);

      if (the_Node->child[0] != NULL) {
        double cx = (minX + maxX) * 0.5;
        double cy = (minY + maxY) * 0.5;
        pushStack(&m_Stack, the_Node->child[0], minX, cx, minY, cy);
        pushStack(&m_Stack, the_Node->child[1], cx, maxX, minY, cy);
        pushStack(&m_Stack, the_Node->child[2], minX, cx, cy, maxY);
        pushStack(&m_Stack, the_Node->child[3], cx, maxX, cy, maxY);
      } else {
        SNodeRange the_Range;
        the_Range.minX = minX;
        the_Range.maxX = maxX;
        the_Range.minY = minY;
        the_Range.maxY = maxY;
        in_Func(the_Node, the_Range, io_UserData);
      }
    }
  }

  void mergeInactiveNodes() { mergeInactiveNodes(m_RootNode); }

protected:
  void freeNode(SQuadTreeNode *io_Node) {
    if (io_Node->child[0] == NULL)
      return;

    for (int i = 0; i < 4; i++) {
      freeNode(io_Node->child[i]);
      free(io_Node->child[i]);
      io_Node->child[i] = NULL;
    }
  }

  void mergeInactiveNodes(SQuadTreeNode *io_Node) {
    if (io_Node->child[0] == NULL)
      return;
    for (int i = 0; i < 4; i++)
      mergeInactiveNodes(io_Node->child[i]);

    if ((io_Node->child[0]->active == AS_INACTIVE) &&
        (io_Node->child[1]->active == AS_INACTIVE) &&
        (io_Node->child[2]->active == AS_INACTIVE) &&
        (io_Node->child[3]->active == AS_INACTIVE)) {
      io_Node->active = AS_INACTIVE;
      for (int i = 0; i < 4; i++) {
        free(io_Node->child[i]);
        io_Node->child[i] = NULL;
      }
    }
  }

  SQuadTreeNode *m_RootNode;
  Vector2s m_MinCoords;
  Vector2s m_MaxCoords;

  SQuadTreeStack m_Stack;
};

inline void countCallBack(SQuadTreeNode *in_Node, SNodeRange &in_Range,
                          void *io_Data) {
  int *c = (int *)io_Data;
  (*c)++;
}

inline int countNonInactiveNodes(CQuadTree *in_QuadTree) {
  int count = 0;
  in_QuadTree->processAllNonInactiveLeafNodes(countCallBack, (void *)&count);
  return count;
}

inline int countNodes(CQuadTree *in_QuadTree) {
  int count = 0;
  in_QuadTree->processAllLeafNodes(countCallBack, (void *)&count);
  return count;
}

#endif

#pragma once

#include <stdlib.h>

#include "bbox.h"
#include "error.h"
#include "mat.h"
#include "perm.h"
#include "ptr_array.h"
#include "tree.h"

typedef struct BfOctree BfOctree;
typedef struct BfOctreeNode BfOctreeNode;

typedef void (*BfOctreeFunc)(BfOctree *, BfOctreeNode *, void *);

enum BfOctreeNodeFlags {
  BF_OCTREE_NODE_FLAG_NONE = 0,
  BF_OCTREE_NODE_FLAG_CHILD_0 = 1 << 0,
  BF_OCTREE_NODE_FLAG_CHILD_1 = 1 << 1,
  BF_OCTREE_NODE_FLAG_CHILD_2 = 1 << 2,
  BF_OCTREE_NODE_FLAG_CHILD_3 = 1 << 3,
  BF_OCTREE_NODE_FLAG_CHILD_4 = 1 << 4,
  BF_OCTREE_NODE_FLAG_CHILD_5 = 1 << 5,
  BF_OCTREE_NODE_FLAG_CHILD_6 = 1 << 6,
  BF_OCTREE_NODE_FLAG_CHILD_7 = 1 << 7,
  BF_OCTREE_NODE_FLAG_ROOT = 1 << 8,

  /* This is a dirty bit which can be used to signal, e.g., that a
   * node has been explored in a tree traversal. The assumption is
   * that this bit is private and will be unset upon entering and
   * exiting any public API function. */
  BF_OCTREE_NODE_FLAG_DIRTY = 1 << 9,
};

struct BfOctree {
  BfTree base;

  /* The root node of the octree. */
  BfOctreeNode *root;

  /* The permutation from the original ordering to the octree
   * ordering. To undo this permutation, `bfPermGetReversePerm` can be
   * used to construct the inverse permutation from the octree
   * ordering back to the original ordering. */
  BfPerm perm;

  /* A pointer to the point set upon which the octree has been
   * built. The points in this array in the original order, and have
   * not been permuted into the octree order. */
  BfPoints3 const *points;

  /* Pointer to an array of unit normals, each associated with a point
   * in `points`. */
  BfVectors3 const *unitNormals;
};

void bfInitOctreeFromPoints(BfOctree *tree, BfPoints3 const *points, BfVectors3 const *unitNormals);
void bfOctreeSaveBoxesToTextFile(BfOctree const *tree, char const *path);
void bfFreeOctree(BfOctree *tree);
void bfGetOctreeNode(BfOctree const *tree, BfSize l, BfSize i, BfOctreeNode **node);
void bfMapOctree(BfOctree *tree, enum BfTreeTraversals traversal, BfOctreeFunc func, void *arg);
void bfMapOctreeNodes(BfOctree *tree, BfOctreeNode *node, enum BfTreeTraversals traversal, BfOctreeFunc func, void *arg);
BfPerm const *bfOctreeGetPerm(BfOctree const *octree);

struct BfOctreeNode {
  enum BfOctreeNodeFlags flags;

  /* The parent of this node. If this is the root node, then this is a
   * pointer to the containing `BfOctree`. Otherwise, it points to
   * the parent node, which will be a `BfOctreeNode`. */
  void *parent;

  /* Index offsets for this node. These are relative to this node's
   * parent. To find the absolute indices, we follow the parent
   * pointer back to the root. The first and last entries hold
   * sentinel values: offset[0] == 0 and offset[4] == size (the number
   * of points contained in this node).
   */
  BfSize offset[5];

  /* Pointers to this node's children. If child[i] == NULL, then no
   * points are contained in the corresponding octree box, and
   * offset[i + 1] - offset[i] == 0 should hold. */
  BfOctreeNode *child[4];

  /* The bounding box for this octree node. All of the points
   * indexed by this node must lie in the rectangle determined by the
   * product of `[bbox.min[0], bbox.max[0]]` and `[bbox.min[1],
   * bbox.max[1]]`. */
  BfBbox3 bbox;

  /* The split for this octree node. If a `point` is in `child[0]`,
   * it satisfies `point[0] <= split[0]` and `point[1] <=
   * split[1]`. If it's in `child[1]`, then `point[0] <= split[0]` and
   * `point[1] > split[1]`. Likewise, if it's in `child[2]`, then
   * `point[0] > split[0]` and `point[1] <= split[1]`; and if it's in
   * `child[3]`, then `point[0] > split[0]` and `point[1]` >
   * `split[1]`. */
  BfPoint3 split;

  /* The depth of this node. */
  BfSize depth;
};

BfSphere bfGetOctreeNodeBoundingSphere(BfOctreeNode const *node);
BfSize bfOctreeNodeNumChildren(BfOctreeNode const *node);
bool bfOctreeNodeIsLeaf(BfOctreeNode const *node);
BfSize bfOctreeNodeDepth(BfOctreeNode const *node);
BfSize bfGetMaxDepthBelowOctreeNode(BfOctreeNode const *node);
BfSize bfOctreeNodeNumPoints(BfOctreeNode const *node);
BfOctree *bfGetOctreeFromNode(BfOctreeNode const *node);

void bfGetOctreeNodePoints(BfOctree const *tree, BfOctreeNode const *node, BfPoints3 *points);
void bfOctreeNodeGetUnitNormals(BfOctree const *tree, BfOctreeNode const *node, BfVectors3 *vectors);
bool bfOctreeNodesAreSeparated(BfOctreeNode const *node1, BfOctreeNode const *node2);

typedef struct BfOctreeLevelIter {
  enum BfTreeTraversals traversal;
  BfPtrArray nodes;
  BfPtrArray levelNodes;
  void *aux;
} BfOctreeLevelIter;

BfOctreeLevelIter bfInitOctreeLevelIter(enum BfTreeTraversals traversal, BfOctreeNode *node);
BfSize bfOctreeLevelIterCurrentDepth(BfOctreeLevelIter const *iter);
bool bfOctreeLevelIterIsDone(BfOctreeLevelIter const *iter);
void bfOctreeLevelIterNext(BfOctreeLevelIter *iter);
void bfFreeOctreeLevelIter(BfOctreeLevelIter *iter);
BfSize bfOctreeLevelIterNumPoints(BfOctreeLevelIter const *iter);

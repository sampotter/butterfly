#pragma once

#include <stdlib.h>

#include "bbox.h"
#include "error.h"
#include "mat.h"
#include "perm.h"
#include "ptr_array.h"
#include "tree.h"

typedef struct BfQuadtree BfQuadtree;
typedef struct BfQuadtreeNode BfQuadtreeNode;

enum BfQuadtreeNodeFlags {
  BF_QUADTREE_NODE_FLAG_NONE = 0,
  BF_QUADTREE_NODE_FLAG_CHILD_0 = 1 << 0,
  BF_QUADTREE_NODE_FLAG_CHILD_1 = 1 << 1,
  BF_QUADTREE_NODE_FLAG_CHILD_2 = 1 << 2,
  BF_QUADTREE_NODE_FLAG_CHILD_3 = 1 << 3,
  BF_QUADTREE_NODE_FLAG_ROOT = 1 << 4,

  /* This is a dirty bit which can be used to signal, e.g., that a
   * node has been explored in a tree traversal. The assumption is
   * that this bit is private and will be unset upon entering and
   * exiting any public API function. */
  BF_QUADTREE_NODE_FLAG_DIRTY = 1 << 5,
};

struct BfQuadtreeNode {
  enum BfQuadtreeNodeFlags flags;

  /* The parent of this node. If this is the root node, then this is a
   * pointer to the containing `BfQuadtree`. Otherwise, it points to
   * the parent node, which will be a `BfQuadtreeNode`. */
  void *parent;

  /* Index offsets for this node. These are relative to this node's
   * parent. To find the absolute indices, we follow the parent
   * pointer back to the root. The first and last entries hold
   * sentinel values: offset[0] == 0 and offset[4] == size (the number
   * of points contained in this node).
   */
  BfSize offset[5];

  /* Pointers to this node's children. If child[i] == NULL, then no
   * points are contained in the corresponding quadtree box, and
   * offset[i + 1] - offset[i] == 0 should hold. */
  BfQuadtreeNode *child[4];

  /* The bounding box for this quadtree node. All of the points
   * indexed by this node must lie in the rectangle determined by the
   * product of `[bbox.min[0], bbox.max[0]]` and `[bbox.min[1],
   * bbox.max[1]]`. */
  BfBbox2 bbox;

  /* The split for this quadtree node. If a `point` is in `child[0]`,
   * it satisfies `point[0] <= split[0]` and `point[1] <=
   * split[1]`. If it's in `child[1]`, then `point[0] <= split[0]` and
   * `point[1] > split[1]`. Likewise, if it's in `child[2]`, then
   * `point[0] > split[0]` and `point[1] <= split[1]`; and if it's in
   * `child[3]`, then `point[0] > split[0]` and `point[1]` >
   * `split[1]`. */
  BfPoint2 split;

  /* The depth of this node. */
  BfSize depth;
};

typedef void (*BfQuadtreeFunc)(BfQuadtree *, BfQuadtreeNode *, void *);

struct BfQuadtree {
  /* A pointer to the point set upon which the quadtree has been
   * built. The points in this array in the original order, and have
   * not been permuted into the quadtree order. */
  BfPoints2 const *points;

  /* The permutation from the original ordering to the quadtree
   * ordering. For example, to put the rows of a matrix `X` which have
   * the same ordering as `points`, call `bfMatPermuteRows(X,
   * &perm)`. To undo this permutation, `bfPermGetReversePerm` can be
   * used to construct the inverse permutation from the quadtree
   * ordering back to the original ordering. */
  BfPerm perm;

  BfQuadtreeNode *root;
};

void bfInitQuadtreeFromPoints(BfQuadtree *tree, BfPoints2 const *points);

void bfInitQuadtreeFromMat(BfQuadtree *tree, BfMat const *mat);

void bfQuadtreeSaveBoxesToTextFile(BfQuadtree const *tree, char const *path);

void bfFreeQuadtree(BfQuadtree *tree);

void
bfGetQuadtreeNode(BfQuadtree const *tree, BfSize l, BfSize i,
                  BfQuadtreeNode **node);

BfCircle2 bfGetQuadtreeNodeBoundingCircle(BfQuadtreeNode const *node);

BfSize bfQuadtreeNodeNumChildren(BfQuadtreeNode const *node);

bool bfQuadtreeNodeIsLeaf(BfQuadtreeNode const *node);

BfSize bfQuadtreeNodeDepth(BfQuadtreeNode const *node);

BfSize bfGetMaxDepthBelowQuadtreeNode(BfQuadtreeNode const *node);

BfSize bfQuadtreeNodeNumPoints(BfQuadtreeNode const *node);

BfQuadtree *bfGetQuadtreeFromNode(BfQuadtreeNode const *node);

void
bfGetQuadtreeNodePoints(BfQuadtree const *tree, BfQuadtreeNode const *node,
                        BfPoints2 *points);

bool bfQuadtreeNodesAreSeparated(BfQuadtreeNode const *node1,
                                 BfQuadtreeNode const *node2);

void
bfMapQuadtree(BfQuadtree *tree, enum BfTreeTraversals traversal,
              BfQuadtreeFunc func, void *arg);

void
bfMapQuadtreeNodes(BfQuadtree *tree, BfQuadtreeNode *node,
                   enum BfTreeTraversals traversal,
                   BfQuadtreeFunc func, void *arg);

typedef struct BfQuadtreeLevelIter {
  enum BfTreeTraversals traversal;
  BfPtrArray nodes;
  BfPtrArray levelNodes;
  void *aux;
} BfQuadtreeLevelIter;

BfQuadtreeLevelIter bfInitQuadtreeLevelIter(enum BfTreeTraversals traversal, BfQuadtreeNode *node);
BfSize bfQuadtreeLevelIterCurrentDepth(BfQuadtreeLevelIter const *iter);
bool bfQuadtreeLevelIterIsDone(BfQuadtreeLevelIter const *iter);
void bfQuadtreeLevelIterNext(BfQuadtreeLevelIter *iter);
void bfFreeQuadtreeLevelIter(BfQuadtreeLevelIter *iter);
BfSize bfQuadtreeLevelIterNumPoints(BfQuadtreeLevelIter const *iter);

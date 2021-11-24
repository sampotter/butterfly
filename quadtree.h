#pragma once

#include <stdlib.h>

#include "def.h"
#include "error.h"
#include "geom.h"

typedef struct BfQuadtree BfQuadtree;
typedef struct BfQuadtreeNode BfQuadtreeNode;

enum BfQuadtreeNodeFlags {
  BF_QUADTREE_NODE_FLAG_ROOT = 0,
  BF_QUADTREE_NODE_FLAG_CHILD_0 = 1 << 0,
  BF_QUADTREE_NODE_FLAG_CHILD_1 = 1 << 1,
  BF_QUADTREE_NODE_FLAG_CHILD_2 = 1 << 2,
  BF_QUADTREE_NODE_FLAG_CHILD_3 = 1 << 3
};

struct BfQuadtreeNode {
  enum BfQuadtreeNodeFlags flags;

  /* The parent of this node. If this is the root node, then this is a
   * pointer to the containing `BfQuadtree`. Otherwise, it points to
   * the parent node, which will be a `BfQuadtreeNode`. */
  void *parent;

  size_t offset[5]; // offset[0] == 0 and offset[4] == size are
                    // sentinel values
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
};

struct BfQuadtree {
  size_t num_points;
  BfPoint2 const *points;

  size_t *perm;

  BfQuadtreeNode *root;

  size_t depth;
};

enum BfError
bfInitQuadtreeFromPoints(BfQuadtree *tree,
                         size_t num_points, BfPoint2 const *points);

enum BfError
bfGetQuadtreeNode(BfQuadtree const *tree, size_t l, size_t i,
                  BfQuadtreeNode **node);

enum BfError
bfGetQuadtreeNodeIndices(BfQuadtreeNode const *node,
                         size_t *num_indices, size_t **indices);

enum BfError
bfGetQuadtreeNodeBoundingBox(BfQuadtreeNode const *node, BfBbox2 *bbox);

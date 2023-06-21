#include <bf/octree_node.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/vectors.h>

#include "macros.h"

#define NUM_CHILDREN_ 8
static BfSize const NUM_CHILDREN = NUM_CHILDREN_;

/** Interface(TreeNode, OctreeNode) */

static BfTreeNodeVtable TreeNodeVtable = {
  .GetType = (__typeof__(&bfTreeNodeGetType))bfOctreeNodeGetType,
  .Delete = (__typeof__(&bfTreeNodeDelete))bfOctreeNodeDelete
};

BfType bfOctreeNodeGetType(BfTreeNode const *treeNode) {
  (void)treeNode;
  return BF_TYPE_OCTREE_NODE;
}

void bfOctreeNodeDelete(BfOctreeNode **octreeNode) {
  bfOctreeNodeDeinit(*octreeNode);
  bfOctreeNodeDealloc(octreeNode);
}

/** Upcasting: OctreeNode -> TreeNode */

BfTreeNode *bfOctreeNodeToTreeNode(BfOctreeNode *node) {
  return &node->super;
}

BfTreeNode const *bfOctreeNodeConstToTreeNodeConst(BfOctreeNode const *node) {
  return &node->super;
}

/** Downcasting: TreeNode -> OctreeNode */

BfOctreeNode *bfTreeNodeToOctreeNode(BfTreeNode *node) {
  if (!bfTreeNodeInstanceOf(node, BF_TYPE_OCTREE_NODE)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfOctreeNode *)node;
  }
}

BfOctreeNode const *bfTreeNodeConstToOctreeNodeConst(BfTreeNode const *node) {
  if (!bfTreeNodeInstanceOf(node, BF_TYPE_OCTREE_NODE)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfOctreeNode const *)node;
  }
}

/** Implementation: OctreeNode */

BfOctreeNode *bfOctreeNodeNew() {
  BF_ERROR_BEGIN();

  BfOctreeNode *node = bfMemAlloc(1, sizeof(BfOctreeNode));
  if (node == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {}

  return node;
}

static void sift(BfSize q, BfSize *offset, BfPoint3 const *point,
                 BfReal const *split, BfSize *perm,
                 bool (*inOctant)(BfPoint3 const, BfReal const *)) {
  BfSize const lastOffset = offset[NUM_CHILDREN];

  BfSize i = offset[q];

  /* Skip over any points at the start of this block indices that are
   * already in the correct octant. This moves `i` to the position of
   * the first point which isn't in the current octant.  */
  while (i < lastOffset && inOctant(point[perm[i]], split))
    ++i;

  BfSize j = (i == lastOffset) ? i : i + 1;

  /* Scan over the remaining points, swapping points in the current
   * octant into position. */
  while (j < lastOffset) {
    if (inOctant(point[perm[j]], split) && !inOctant(point[perm[i]], split)) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }

  offset[q + 1] = i;
}

static bool inOctant1(BfPoint3 const x, BfReal const *split) {
  return x[0] <= split[0] && x[1] <= split[1] && x[2] <= split[2];
}

static bool inOctant2(BfPoint3 const x, BfReal const *split) {
  return x[0] <= split[0] && x[1] <= split[1] && x[2] > split[2];
}

static bool inOctant3(BfPoint3 const x, BfReal const *split) {
  return x[0] <= split[0] && x[1] > split[1] && x[2] <= split[2];
}

static bool inOctant4(BfPoint3 const x, BfReal const *split) {
  return x[0] <= split[0] && x[1] > split[1] && x[2] > split[2];
}

static bool inOctant5(BfPoint3 const x, BfReal const *split) {
  return x[0] > split[0] && x[1] <= split[1] && x[2] <= split[2];
}

static bool inOctant6(BfPoint3 const x, BfReal const *split) {
  return x[0] > split[0] && x[1] <= split[1] && x[2] > split[2];
}

static bool inOctant7(BfPoint3 const x, BfReal const *split) {
  return x[0] > split[0] && x[1] > split[1] && x[2] <= split[2];
}

static bool inOctant8(BfPoint3 const x, BfReal const *split) {
  return x[0] > split[0] && x[1] > split[1] && x[2] > split[2];
}

static __typeof__(&inOctant1) inOctant[NUM_CHILDREN_] = {
  inOctant1, inOctant2, inOctant3, inOctant4,
  inOctant5, inOctant6, inOctant7, inOctant8
};

/* Initialize a octree node, and allocate and initialize each node
 * beneath this node. Calling this function for the octree's root
 * node has the effect of building the complete octree.
 *
 * `node`: the root node at this level of the recursion.
 * `points`: the entire point set.
 * `boundingBox`: `node`'s bounding box.
 * `i0`, `i1`: `[perm[i0], ..., perm[i1])` indexes the points
 *   contained by `node`. After calling `recInitOctreeNode`,
 *   these entries of `perm` will be put into Z order.
 * `perm`: an array of `points->size` indices indexing `points`.
 * `currentDepth`: `node`'s depth */
static void octreeNodeInitRecursive(BfOctreeNode *node,
                                    BfPoints3 const *points, BfSize maxLeafSize,
                                    BfBoundingBox3 boundingBox,
                                    BfSize i0, BfSize i1, BfSize *perm,
                                    BfSize currentDepth) {
  BF_ERROR_BEGIN();

  BF_ASSERT(i0 <= i1);

  BfPoint3 const *point = (BfPoint3 const *)points->data;
  BfReal const *split = node->split;
  BfTreeNode **child = &node->super.child[0];
  BfSize *offset = &node->super.offset[0];

  /* Sanity check: make sure all points are in bounding box for node. */
#if BF_DEBUG
  for (BfSize i = i0; i < i1; ++i)
    BF_ASSERT(bfBoundingBox3ContainsPoint(&boundingBox, point[perm[i]]));
#endif

  /* Sift the points into the correct octants at this level. Note that
   * because of the way we do the sifting, it's unnecessary to call
   * `sift` for the last child node. */
  offset[0] = i0;
  offset[NUM_CHILDREN] = i1;
  for (BfSize q = 0; q < NUM_CHILDREN - 1; ++q)
    sift(q, offset, point, split, perm, inOctant[q]);

  /* Check that each point is in the correct octant now. */
#ifdef BF_DEBUG
  for (BfSize q = 0; q < NUM_CHILDREN; ++q)
    for (BfSize i = offset[q]; i < offset[q + 1]; ++i)
      BF_ASSERT(inOctant[q](point[perm[i]], split));
#endif

  /* Compute the bounding boxes each child node */
  BfBoundingBox3 childBoundingBox[NUM_CHILDREN_] = {
    [0] = {
      .min = {boundingBox.min[0], boundingBox.min[1], boundingBox.min[2]},
      .max = {          split[0],           split[1],           split[2]}
    },
    [1] = {
      .min = {boundingBox.min[0], boundingBox.min[1],           split[2]},
      .max = {          split[0],           split[1], boundingBox.max[2]}
    },
    [2] = {
      .min = {boundingBox.min[0],           split[1], boundingBox.min[2]},
      .max = {          split[0], boundingBox.max[1],           split[2]}
    },
    [3] = {
      .min = {boundingBox.min[0],           split[1],           split[2]},
      .max = {          split[0], boundingBox.max[1], boundingBox.max[2]}
    },
    [4] = {
      .min = {          split[0], boundingBox.min[1], boundingBox.min[2]},
      .max = {boundingBox.max[0],           split[1],           split[2]}
    },
    [5] = {
      .min = {          split[0], boundingBox.min[1],           split[2]},
      .max = {boundingBox.max[0],           split[1], boundingBox.max[2]}
    },
    [6] = {
      .min = {          split[0],           split[1], boundingBox.min[2]},
      .max = {boundingBox.max[0], boundingBox.max[1],           split[2]}
    },
    [7] = {
      .min = {          split[0],           split[1],           split[2]},
      .max = {boundingBox.max[0], boundingBox.max[1], boundingBox.max[2]}
    },
  };

  for (BfSize q = 0; q < NUM_CHILDREN; ++q) {
    BfSize numChildPoints = offset[q + 1] - offset[q];
    if (numChildPoints == 0)
      continue;

    BF_ASSERT(!bfBoundingBox3IsEmpty(&childBoundingBox[q]));
    BF_ASSERT(bfBoundingBox3GetVolume(&childBoundingBox[q]) >= 1e-15);

    /* Create and initialize a new octree node */
    BfOctreeNode *newChild = bfOctreeNodeNew();
    bfTreeNodeInit(&newChild->super, &TreeNodeVtable, false, (void *)node,
                   NUM_CHILDREN, q, currentDepth + 1);

    /* Compute bounding box and split for the `q`th child node */
    newChild->boundingBox = childBoundingBox[q];
    bfBoundingBox3GetCenter(&childBoundingBox[q], newChild->split);

    /* If the node has few enough points, it's a leaf and no more
     * initialization needs to be done. Otherwise, we continue
     * building the octree recursively. */
    if (numChildPoints > maxLeafSize) {
      octreeNodeInitRecursive(
        newChild, points, maxLeafSize,
        childBoundingBox[q], offset[q], offset[q + 1],
        perm, currentDepth + 1);
      HANDLE_ERROR();
    }

    /* Set the `q`th child to `newChild` for the current node */
    child[q] = bfOctreeNodeToTreeNode(newChild);
  }

  /* Sanity check: a child is `NULL` exactly when there are *no*
   * points in the corresponding index range. */
#if BF_DEBUG
  for (BfSize q = 0; q < NUM_CHILDREN; ++q)
    BF_ASSERT((offset[q] == offset[q + 1] && child[q] == NULL) ||
              (offset[q] < offset[q + 1] && child[q] != NULL));
#endif

  BF_ERROR_END() {
    for (BfSize q = 0; q < NUM_CHILDREN; ++q) {
      bfMemFree(child[q]);
    }
  }
}

void bfOctreeNodeInitRoot(BfOctreeNode *node, BfOctree const *tree, BfSize maxLeafSize) {
  BF_ERROR_BEGIN();

  bfTreeNodeInit(&node->super, &TreeNodeVtable,
                 true, (void *)tree, NUM_CHILDREN, BF_SIZE_BAD_VALUE, 0);
  HANDLE_ERROR();

  /* Compute scaled bounding box for the entire octree */
  BfBoundingBox3 boundingBox = bfPoints3GetBoundingBox(tree->points);
  bfBoundingBox3RescaleToCube(&boundingBox);
  node->boundingBox = boundingBox;

  /* Perturb the sides of the bounding box slightly to make sure all
   * points are inside: */
  for (BfSize i = 0; i < 3; ++i) boundingBox.min[i] -= 1e2*BF_EPS;
  for (BfSize i = 0; i < 3; ++i) boundingBox.max[i] += 1e2*BF_EPS;

  /* Compute the split for the node */
  bfBoundingBox3GetCenter(&boundingBox, node->split);

  octreeNodeInitRecursive(node, tree->points, maxLeafSize,
                          boundingBox, 0, tree->points->size,
                          tree->super.perm.index, 0);
  HANDLE_ERROR();

  BF_ERROR_END() {}
}

void bfOctreeNodeDeinit(BfOctreeNode *octreeNode) {
  bfTreeNodeDeinit(&octreeNode->super);
}

void bfOctreeNodeDealloc(BfOctreeNode **octreeNode) {
  bfMemFree(*octreeNode);
  *octreeNode = NULL;
}

// BfOctreeNode *bfOctreeNodeGetChild(BfOctreeNode *node, BfSize i);
// BfOctreeNode const *bfOctreeNodeGetChildConst(BfOctreeNode const *node, BfSize i);

BfSphere bfOctreeNodeGetBoundingSphere(BfOctreeNode const *node) {
  return bfBoundingBox3GetBoundingSphere(&node->boundingBox);
}

/* Fill `points` with the points contained in `node`. They will be
 * added in octree order.
 *
 * If `tree == NULL`, then this function will retrieve the containing
 * `BfOctree` from `node`, which takes `O(log N)` time. */
BfPoints3 bfOctreeNodeGetPoints(BfOctreeNode const *octreeNode, BfOctree const *octree) {
  BF_ERROR_BEGIN();

  BfPoints3 points;

  BfTreeNode const *treeNode = bfOctreeNodeConstToTreeNodeConst(octreeNode);

  BfTree const *tree = octree == NULL ?
    bfTreeNodeGetTreeConst(treeNode) :
    bfOctreeConstToTreeConst(octree);

  /* determine the number of points contained by `node` and find the
   * offset into `tree->perm` */
  BfSize numInds = bfTreeNodeGetNumPoints(treeNode);
  BfSize const *inds = bfTreeNodeGetIndexPtrConst(treeNode, tree);

  bfPoints3GetByIndex(octree->points, numInds, inds, &points);
  HANDLE_ERROR();

  BF_ERROR_END() {}

  return points;
}

BfVectors3 bfOctreeNodeGetUnitNormals(BfOctreeNode const *octreeNode, BfOctree const *octree) {
  BF_ERROR_BEGIN();

  if (octree->unitNormals == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfTreeNode const *treeNode = bfOctreeNodeConstToTreeNodeConst(octreeNode);

  BfTree const *tree = octree == NULL ?
    bfTreeNodeGetTreeConst(treeNode) :
    bfOctreeConstToTreeConst(octree);

  /* determine the number of points contained by `node` and find the
   * offset into `tree->perm` */
  BfSize numInds = bfTreeNodeGetNumPoints(treeNode);
  BfSize const *inds = bfTreeNodeGetIndexPtrConst(treeNode, tree);

  BfVectors3 unitNormals;
  bfVectors3GetByIndex(octree->unitNormals, numInds, inds, &unitNormals);
  HANDLE_ERROR();

  BF_ERROR_END() {}

  return unitNormals;
}

bool bfOctreeNodesAreSeparated(BfOctreeNode const *node1, BfOctreeNode const *node2) {
  BfSphere sphere1 = bfOctreeNodeGetBoundingSphere(node1);
  BfSphere sphere2 = bfOctreeNodeGetBoundingSphere(node2);

  BfReal R = bfPoint3Dist(sphere1.center, sphere2.center);

  return R > sphere1.r + sphere2.r + 1e1*BF_EPS_MACH;
}

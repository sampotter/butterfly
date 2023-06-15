#include <bf/quadtree_node.h>

#include <math.h>

#include <bf/assert.h>
#include <bf/circle.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/points.h>
#include <bf/vectors.h>

#include "macros.h"

static BfSize const NUM_CHILDREN = 4;

static BfSize const LEAF_SIZE_THRESHOLD = 1;

/** Interface(TreeNode, QuadtreeNode) */

static BfTreeNodeVtable TreeNodeVtable = {
  .GetType = (__typeof__(&bfTreeNodeGetType))bfQuadtreeNodeGetType,
  .Delete = (__typeof__(&bfTreeNodeDelete))bfQuadtreeNodeDelete,
};

BfType bfQuadtreeNodeGetType(BfTreeNode const *treeNode) {
  (void)treeNode;
  return BF_TYPE_QUADTREE_NODE;
}

void bfQuadtreeNodeDelete(BfQuadtreeNode **quadtreeNode) {
  bfQuadtreeNodeDeinit(*quadtreeNode);
  bfQuadtreeNodeDealloc(quadtreeNode);
}

/** Upcasting: QuadtreeNode -> TreeNode */

BfTreeNode *bfQuadtreeNodeToTreeNode(BfQuadtreeNode *node) {
  return &node->super;
}

BfTreeNode const *bfQuadtreeNodeConstToTreeNodeConst(BfQuadtreeNode const *node) {
  return &node->super;
}

/** Downcasting: TreeNode -> QuadtreeNode */

BfQuadtreeNode *bfTreeNodeToQuadtreeNode(BfTreeNode *node) {
  if (!bfTreeNodeInstanceOf(node, BF_TYPE_QUADTREE_NODE)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfQuadtreeNode *)node;
  }
}

BfQuadtreeNode const *bfTreeNodeConstToQuadtreeNodeConst(BfTreeNode const *node) {
  if (!bfTreeNodeInstanceOf(node, BF_TYPE_QUADTREE_NODE)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfQuadtreeNode const *)node;
  }
}

/** Implementation: QuadtreeNode */

BfQuadtreeNode *bfQuadtreeNodeNew() {
  BF_ERROR_BEGIN();

  BfQuadtreeNode *node = bfMemAlloc(1, sizeof(BfQuadtreeNode));
  HANDLE_ERROR();

  BF_ERROR_END() {}

  return node;
}

/* Initialize a quadtree node, and allocate and initialize each node
 * beneath this node. Calling this function for the quadtree's root
 * node has the effect of building the complete quadtree.
 *
 * `node`: the root node at this level of the recursion.
 * `points`: the entire point set.
 * `bbox`: `node`'s bounding box.
 * `i0`, `i1`: `[perm[i0], ..., perm[i1])` indexes the points
 *   contained by `node`. After calling `recInitQuadtreeNode`,
 *   these entries of `perm` will be put into Z order.
 * `perm`: an array of `points->size` indices indexing `points`.
 * `currentDepth`: `node`'s depth */
static void quadtreeNodeInitRecursive(BfQuadtreeNode *node,
                                      BfPoints2 const *points, BfBbox2 bbox,
                                      BfSize i0, BfSize i1, BfSize *perm,
                                      BfSize currentDepth) {
  BF_ERROR_BEGIN();

  BF_ASSERT(i0 <= i1);

  BfPoint2 *point = points->data;
  BfReal const *split = node->split;

  BfTreeNode **child = &node->super.child[0];
  BfSize *offset = &node->super.offset[0];

  /* Set sentinel values */
  offset[0] = i0;
#if BF_DEBUG
  offset[1] = BF_SIZE_BAD_VALUE;
  offset[2] = BF_SIZE_BAD_VALUE;
  offset[3] = BF_SIZE_BAD_VALUE;
#endif
  offset[4] = i1;

  BfSize i, j;

  /* Sift permutation indices for child[0] into place... */
  i = offset[0];
  while (i < offset[4] &&
         (point[perm[i]][0] <= split[0] && point[perm[i]][1] <= split[1]))
    ++i;
  j = i == offset[4] ? i : i + 1;
  while (j < offset[4]) {
    if (point[perm[j]][0] <= split[0] && point[perm[j]][1] <= split[1] &&
        (point[perm[i]][0] > split[0] || point[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  offset[1] = i;

  /* ... For child[1]... */
  i = offset[1]; // TODO: unnecessary
  while (i < offset[4] &&
         (point[perm[i]][0] <= split[0] && point[perm[i]][1] > split[1]))
    ++i;
  j = i == offset[4] ? i : i + 1;
  while (j < offset[4]) {
    if (point[perm[j]][0] <= split[0] && point[perm[j]][1] > split[1] &&
        (point[perm[i]][0] > split[0] || point[perm[i]][1] <= split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  offset[2] = i;

  /* ... For child[2]... */
  i = offset[2]; // TODO: unnecessary
  while (i < offset[4] &&
         (point[perm[i]][0] > split[0] && point[perm[i]][1] <= split[1]))
    ++i;
  j = i == offset[4] ? i : i + 1;
  while (j < offset[4]) {
    if (point[perm[j]][0] > split[0] && point[perm[j]][1] <= split[1] &&
        (point[perm[i]][0] <= split[0] || point[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  offset[3] = i;

  /* we don't have to do any sifting for the last child! ... but
   * should verify anyway */

  /* Compute the bounding boxes each child node */
  BfBbox2 childBbox[4] = {
    [0] = {
      .min = {bbox.min[0], bbox.min[1]},
      .max = {split[0], split[1]}
    },
    [1] = {
      .min = {bbox.min[0], split[1]},
      .max = {split[0], bbox.max[1]}
    },
    [2] = {
      .min = {split[0], bbox.min[1]},
      .max = {bbox.max[0], split[1]}
    },
    [3] = {
      .min = {split[0], split[1]},
      .max = {bbox.max[0], bbox.max[1]}
    }
  };

#ifdef BF_DEBUG
  /* Sanity checks follow: */

  BF_ASSERT(j == offset[4]);

  /* child[0]'s indices are correctly sifted */
  for (BfSize k = offset[0]; k < offset[1]; ++k)
    BF_ASSERT(point[perm[k]][0] <= split[0] && point[perm[k]][1] <= split[1]);

  /* child[1]'s indices are correctly sifted */
  for (BfSize k = offset[1]; k < offset[2]; ++k)
    BF_ASSERT(point[perm[k]][0] <= split[0] && point[perm[k]][1] > split[1]);

  /* child[2]'s indices are correctly sifted */
  for (BfSize k = offset[2]; k < offset[3]; ++k)
    BF_ASSERT(point[perm[k]][0] > split[0] && point[perm[k]][1] <= split[1]);

  /* child[3]'s indices are correctly sifted */
  for (BfSize k = offset[3]; k < offset[4]; ++k)
    BF_ASSERT(point[perm[k]][0] > split[0] && point[perm[k]][1] > split[1]);
#endif

  for (BfSize q = 0; q < NUM_CHILDREN; ++q) {
    BfSize numChildPoints = offset[q + 1] - offset[q];
    if (numChildPoints == 0)
      continue;

    BfQuadtreeNode *newChild = bfQuadtreeNodeNew();
    HANDLE_ERROR();

    bfTreeNodeInit(&newChild->super, &TreeNodeVtable, false, (void *)node,
                   NUM_CHILDREN, q, currentDepth + 1);
    HANDLE_ERROR();

    /* Compute bounding box and split for the `q`th child node */
    newChild->bbox = childBbox[q];
    bfBbox2GetCenter(&childBbox[q], newChild->split);

    /* If the node has few enough points, it's a leaf and no more
     * initialization needs to be done. Otherwise, we continue
     * building the quadtree recursively. */
    if (numChildPoints > LEAF_SIZE_THRESHOLD) {
      quadtreeNodeInitRecursive(
        newChild, points, childBbox[q], offset[q], offset[q + 1],
        perm, currentDepth + 1);
      HANDLE_ERROR();
    }

    /* Set the `q`th child to `newChild` for the current node */
    child[q] = bfQuadtreeNodeToTreeNode(newChild);
  }

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

void bfQuadtreeNodeInitRoot(BfQuadtreeNode *node, BfQuadtree const *tree) {
  BF_ERROR_BEGIN();

  bfTreeNodeInit(&node->super, &TreeNodeVtable,
                 true, (void *)tree, NUM_CHILDREN, BF_SIZE_BAD_VALUE, 0);
  HANDLE_ERROR();

  /* Compute scaled bounding box for the entire quadtree */
  BfBbox2 bbox = bfPoints2GetBoundingBox(tree->points);
  bfBbox2RescaleToSquare(&bbox);
  node->bbox = bbox;

  /* Compute the split for the node node */
  bfBbox2GetCenter(&bbox, node->split);

  quadtreeNodeInitRecursive(node, tree->points, bbox, 0, tree->points->size,
                            tree->super.perm.index, 0);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfQuadtreeNodeDeinit(BfQuadtreeNode *quadtreeNode) {
  bfTreeNodeDeinit(&quadtreeNode->super);
}

void bfQuadtreeNodeDealloc(BfQuadtreeNode **quadtreeNode) {
  bfMemFree(*quadtreeNode);
  *quadtreeNode = NULL;
}

// void bfQuadtreeNodeDeinitAndDealloc(BfQuadtreeNode **node);
// BfQuadtreeNode *bfQuadtreeNodeGetChild(BfQuadtreeNode *node, BfSize i);
// BfQuadtreeNode const *bfQuadtreeNodeGetChildConst(BfQuadtreeNode const *node, BfSize i);

BfCircle bfQuadtreeNodeGetBoundingCircle(BfQuadtreeNode const *node) {
  BfReal const *min = node->bbox.min, *max = node->bbox.max;

  BfCircle circ;
  circ.r = hypot(max[0] - min[0], max[1] - min[1])/2;
  circ.center[0] = (min[0] + max[0])/2;
  circ.center[1] = (min[1] + max[1])/2;

  return circ;
}

/* Fill `points` with the points contained in `node`. They will be
 * added in quadtree order.
 *
 * If `tree == NULL`, then this function will retrieve the containing
 * `BfQuadtree` from `node`, which takes `O(log N)` time. */
BfPoints2 bfQuadtreeNodeGetPoints(BfQuadtreeNode const *quadtreeNode,
                                  BfQuadtree const *quadtree) {
  BF_ERROR_BEGIN();

  BfPoints2 points;

  BfTreeNode const *treeNode = bfQuadtreeNodeConstToTreeNodeConst(quadtreeNode);

  BfTree const *tree = quadtree == NULL ?
    bfTreeNodeGetTreeConst(treeNode) :
    bfQuadtreeConstToTreeConst(quadtree);

  /* determine the number of points contained by `node` and find the
   * offset into `tree->perm` */
  BfSize numInds = bfTreeNodeGetNumPoints(treeNode);
  BfSize const *inds = bfTreeNodeGetIndexPtrConst(treeNode, tree);

  bfGetPointsByIndex(quadtree->points, numInds, inds, &points);
  HANDLE_ERROR();

  BF_ERROR_END() {}

  return points;
}

BfVectors2 *bfQuadtreeNodeGetUnitNormals(BfQuadtreeNode const *quadtreeNode,
                                         BfQuadtree const *quadtree) {
  BF_ERROR_BEGIN();

  if (quadtree->unitNormals == NULL)
    return NULL;

  BfVectors2 *unitNormals = bfMemAlloc(1, sizeof(BfVectors2));
  HANDLE_ERROR();

  BfTreeNode const *treeNode = bfQuadtreeNodeConstToTreeNodeConst(quadtreeNode);

  BfTree const *tree = quadtree == NULL ?
    bfTreeNodeGetTreeConst(treeNode) :
    bfQuadtreeConstToTreeConst(quadtree);

  /* determine the number of points containined by `node` and find the
   * offset into `tree->perm` */
  BfSize numInds = bfTreeNodeGetNumPoints(treeNode);
  BfSize const *inds = bfTreeNodeGetIndexPtrConst(treeNode, tree);

  bfGetVectorsByIndex(quadtree->unitNormals, numInds, inds, unitNormals);
  HANDLE_ERROR();

  BF_ERROR_END() {}

  return unitNormals;
}

bool bfQuadtreeNodesAreSeparated(BfQuadtreeNode const *node1,
                                 BfQuadtreeNode const *node2) {
  BfCircle circ1 = bfQuadtreeNodeGetBoundingCircle(node1);
  BfCircle circ2 = bfQuadtreeNodeGetBoundingCircle(node2);

  BfReal R = bfPoint2Dist(circ1.center, circ2.center);

  return R > circ1.r + circ2.r + 1e1*BF_EPS_MACH;
}

#include <bf/octree_node.h>

#include <bf/bbox.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/points.h>

#include "macros.h"

static BfSize const NUM_CHILDREN = 8;

/** Interface: TreeNode -> OctreeNode */

static BfTreeNodeVtable TreeNodeVtable = {
  .GetType = (__typeof__(&bfTreeNodeGetType))bfOctreeNodeGetType
};

BfType bfOctreeNodeGetType(BfTreeNode const *treeNode) {
  (void)treeNode;
  return BF_TYPE_OCTREE_NODE;
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

BfOctreeNode *bfOctreeNodeNew() {
  BEGIN_ERROR_HANDLING();

  BfOctreeNode *node = malloc(sizeof(BfOctreeNode));
  if (node == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    node = NULL;

  return node;
}

#define IN_CHILD_0(x) (x[0] <= split[0] && x[1] <= split[1] && x[2] <= split[2])
#define IN_CHILD_1(x) (x[0] <= split[0] && x[1] <= split[1] && x[2] >  split[2])
#define IN_CHILD_2(x) (x[0] <= split[0] && x[1] >  split[1] && x[2] <= split[2])
#define IN_CHILD_3(x) (x[0] <= split[0] && x[1] >  split[1] && x[2] >  split[2])
#define IN_CHILD_4(x) (x[0] >  split[0] && x[1] <= split[1] && x[2] <= split[2])
#define IN_CHILD_5(x) (x[0] >  split[0] && x[1] <= split[1] && x[2] >  split[2])
#define IN_CHILD_6(x) (x[0] >  split[0] && x[1] >  split[1] && x[2] <= split[2])
#define IN_CHILD_7(x) (x[0] >  split[0] && x[1] >  split[1] && x[2] >  split[2])

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
static void octreeNodeInitRecursive(
  BfOctreeNode *node, BfPoints3 const *points, BfBoundingBox3 const *boundingBox,
  BfSize i0, BfSize i1, BfSize *perm, BfSize currentDepth)
{
  BEGIN_ERROR_HANDLING();

  assert(i0 <= i1);

  BfPoint3 *point = points->data;
  BfReal const *split = node->split;
  node->super.depth = currentDepth;

  BfTreeNode **child = &node->super.child[0];
  BfSize *offset = &node->super.offset[0];

  /* Set sentinel values */
  offset[0] = i0;
#if BF_DEBUG
  for (BfSize q = 1; q < NUM_CHILDREN; ++q)
    offset[q] = BF_SIZE_BAD_VALUE;
#endif
  offset[8] = i1;

#define SIFT_POINTS(LAST_OFFSET, CURRENT_OFFSET) do {               \
    BfSize i = offset[CURRENT_OFFSET];                        \
    while (i < offset[LAST_OFFSET] &&                         \
           IN_CHILD(point[perm[i]]))                                \
      ++i;                                                          \
    BfSize j = i;                                                   \
    assert(j <= offset[LAST_OFFSET]);                         \
    if (j < offset[LAST_OFFSET])                              \
      ++j;                                                          \
    while (j < offset[LAST_OFFSET]) {                         \
      if (IN_CHILD(point[perm[j]]) && !IN_CHILD(point[perm[i]])) {  \
        SWAP(perm[j], perm[i]);                                     \
        ++i;                                                        \
      }                                                             \
      ++j;                                                          \
    }                                                               \
    offset[CURRENT_OFFSET + 1] = i;                           \
  } while (0)

  /* Sift all children */
#define IN_CHILD(x) IN_CHILD_0(x)
  SIFT_POINTS(8, 0);
#undef IN_CHILD
#define IN_CHILD(x) IN_CHILD_1(x)
  SIFT_POINTS(8, 1);
#undef IN_CHILD
#define IN_CHILD(x) IN_CHILD_2(x)
  SIFT_POINTS(8, 2);
#undef IN_CHILD
#define IN_CHILD(x) IN_CHILD_3(x)
  SIFT_POINTS(8, 3);
#undef IN_CHILD
#define IN_CHILD(x) IN_CHILD_4(x)
  SIFT_POINTS(8, 4);
#undef IN_CHILD
#define IN_CHILD(x) IN_CHILD_5(x)
  SIFT_POINTS(8, 5);
#undef IN_CHILD
#define IN_CHILD(x) IN_CHILD_6(x)
  SIFT_POINTS(8, 6);
#undef IN_CHILD
  /* (... last child doesn't need to be sifted!) */

  /* Compute the bounding boxes each child node */
  BfBoundingBox3 childBoundingBox[8] = {
    [0] = {
      .min = {boundingBox->min[0], boundingBox->min[1], boundingBox->min[2]},
      .max = {split[0],            split[1],            split[2]}
    },
    [1] = {
      .min = {boundingBox->min[0], boundingBox->min[1], split[2]},
      .max = {split[0],            split[1],            boundingBox->max[2]}
    },
    [2] = {
      .min = {boundingBox->min[0], split[1],            boundingBox->min[2]},
      .max = {split[0],            boundingBox->max[1], split[2]}
    },
    [3] = {
      .min = {boundingBox->min[0], split[1],            split[2]},
      .max = {split[0],            boundingBox->max[1], boundingBox->max[2]}
    },
    [4] = {
      .min = {split[0],            boundingBox->min[1], boundingBox->min[2]},
      .max = {boundingBox->max[0], split[1],            split[2]}
    },
    [5] = {
      .min = {split[0],            boundingBox->min[1], split[2]},
      .max = {boundingBox->max[0], split[1],            boundingBox->max[2]}
    },
    [6] = {
      .min = {split[0],            split[1],            boundingBox->min[2]},
      .max = {boundingBox->max[0], boundingBox->max[1], split[2]}
    },
    [7] = {
      .min = {split[0],            split[1],            split[2]},
      .max = {boundingBox->max[0], boundingBox->max[1], boundingBox->max[2]}
    }
  };

#ifdef BF_DEBUG
  for (BfSize k = offset[0]; k < offset[1]; ++k)
    assert(IN_CHILD_0(point[perm[k]]));
  for (BfSize k = offset[1]; k < offset[2]; ++k)
    assert(IN_CHILD_1(point[perm[k]]));
  for (BfSize k = offset[2]; k < offset[3]; ++k)
    assert(IN_CHILD_2(point[perm[k]]));
  for (BfSize k = offset[3]; k < offset[4]; ++k)
    assert(IN_CHILD_3(point[perm[k]]));
  for (BfSize k = offset[4]; k < offset[5]; ++k)
    assert(IN_CHILD_4(point[perm[k]]));
  for (BfSize k = offset[5]; k < offset[6]; ++k)
    assert(IN_CHILD_5(point[perm[k]]));
  for (BfSize k = offset[6]; k < offset[7]; ++k)
    assert(IN_CHILD_6(point[perm[k]]));
  for (BfSize k = offset[7]; k < offset[8]; ++k)
    assert(IN_CHILD_7(point[perm[k]]));
#endif

  for (BfSize q = 0; q < 8; ++q) {
    BfSize childPermSize = offset[q + 1] - offset[q];

    /* TODO: right now, building octree to maximum depth possible */
    if (childPermSize <= 1) {
      child[q] = NULL;
      continue;
    }

    BfOctreeNode *newChild = bfOctreeNodeNew();

    bfTreeNodeInit(&newChild->super, &TreeNodeVtable, false, (void *)node,
                   NUM_CHILDREN, q, currentDepth + 1);

    newChild->boundingBox = childBoundingBox[q];

    /* Compute the split for the `q`th child node */
    newChild->split[0] = (childBoundingBox[q].min[0] + childBoundingBox[q].max[0])/2;
    newChild->split[1] = (childBoundingBox[q].min[1] + childBoundingBox[q].max[1])/2;
    newChild->split[2] = (childBoundingBox[q].min[2] + childBoundingBox[q].max[2])/2;

    octreeNodeInitRecursive(newChild, points, &newChild->boundingBox,
                            offset[q], offset[q + 1], perm, currentDepth + 1);
    HANDLE_ERROR();

    child[q] = bfOctreeNodeToTreeNode(newChild);
  }

  END_ERROR_HANDLING() {
    for (BfSize q = 0; q < NUM_CHILDREN; ++q)
      free(child[q]);
  }
}

void bfOctreeNodeInitRoot(BfOctreeNode *node, BfOctree const *tree) {
  BEGIN_ERROR_HANDLING();

  bfTreeNodeInit(&node->super, &TreeNodeVtable, true, (void *)tree,
                 NUM_CHILDREN, BF_SIZE_BAD_VALUE, 0);
  HANDLE_ERROR();

  /* Compute the bounding box for the entire octree */
  node->boundingBox = bfPoints3GetBoundingBox(tree->points);
  assert(!bfBoundingBox3IsEmpty(&node->boundingBox));

  /* Rescale the bounding box so that it's square */
  bfBoundingBox3RescaleToCube(&node->boundingBox);

  /* Compute the split for the root node */
  bfBoundingBox3GetCenter(&node->boundingBox, node->split);

  /* Recursively initialize octree starting from root */
  octreeNodeInitRecursive(node, tree->points, &node->boundingBox, 0,
                          tree->points->size, tree->super.perm.index, 0);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

void bfOctreeNodeDeinit(BfOctreeNode *node) {
  (void)node;
  assert(false);
}

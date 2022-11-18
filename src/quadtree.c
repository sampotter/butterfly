#include <bf/quadtree.h>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <bf/circle.h>
#include <bf/error_macros.h>
#include <bf/points.h>
#include <bf/ptr_array.h>
#include <bf/vectors.h>

#define SWAP(x, y) do {                         \
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = tmp;                                    \
  } while (0)

enum BfQuadtreeNodeFlags childFlags[4] = {
  BF_QUADTREE_NODE_FLAG_CHILD_0,
  BF_QUADTREE_NODE_FLAG_CHILD_1,
  BF_QUADTREE_NODE_FLAG_CHILD_2,
  BF_QUADTREE_NODE_FLAG_CHILD_3
};

enum BfQuadtreeNodeFlags child_mask =
  BF_QUADTREE_NODE_FLAG_CHILD_0 |
  BF_QUADTREE_NODE_FLAG_CHILD_1 |
  BF_QUADTREE_NODE_FLAG_CHILD_2 |
  BF_QUADTREE_NODE_FLAG_CHILD_3;

BfSize child_flag_to_index[] = {
  [BF_QUADTREE_NODE_FLAG_CHILD_0] = 0,
  [BF_QUADTREE_NODE_FLAG_CHILD_1] = 1,
  [BF_QUADTREE_NODE_FLAG_CHILD_2] = 2,
  [BF_QUADTREE_NODE_FLAG_CHILD_3] = 3
};

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
static
void
recInitQuadtreeNode(BfQuadtreeNode *node,
                    BfPoints2 const *points, BfBbox2 bbox,
                    BfSize i0, BfSize i1, BfSize *perm,
                    BfSize currentDepth)
{
  BEGIN_ERROR_HANDLING();

  assert(i0 <= i1);

  BfPoint2 *point = points->data;
  BfReal const *split = node->split;
  node->depth = currentDepth;

  /* Set sentinel values */
  node->offset[0] = i0;
#if BF_DEBUG
  node->offset[1] = BF_SIZE_BAD_VALUE;
  node->offset[2] = BF_SIZE_BAD_VALUE;
  node->offset[3] = BF_SIZE_BAD_VALUE;
#endif
  node->offset[4] = i1;

  BfSize i, j;

  /* Sift permutation indices for child[0] into place... */
  i = node->offset[0];
  while (i < node->offset[4] &&
         (point[perm[i]][0] <= split[0] && point[perm[i]][1] <= split[1]))
    ++i;
  j = i == node->offset[4] ? i : i + 1;
  while (j < node->offset[4]) {
    if (point[perm[j]][0] <= split[0] && point[perm[j]][1] <= split[1] &&
        (point[perm[i]][0] > split[0] || point[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[1] = i;

  /* ... For child[1]... */
  i = node->offset[1]; // TODO: unnecessary
  while (i < node->offset[4] &&
         (point[perm[i]][0] <= split[0] && point[perm[i]][1] > split[1]))
    ++i;
  j = i == node->offset[4] ? i : i + 1;
  while (j < node->offset[4]) {
    if (point[perm[j]][0] <= split[0] && point[perm[j]][1] > split[1] &&
        (point[perm[i]][0] > split[0] || point[perm[i]][1] <= split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[2] = i;

  /* ... For child[2]... */
  i = node->offset[2]; // TODO: unnecessary
  while (i < node->offset[4] &&
         (point[perm[i]][0] > split[0] && point[perm[i]][1] <= split[1]))
    ++i;
  j = i == node->offset[4] ? i : i + 1;
  while (j < node->offset[4]) {
    if (point[perm[j]][0] > split[0] && point[perm[j]][1] <= split[1] &&
        (point[perm[i]][0] <= split[0] || point[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[3] = i;

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

  assert(j == node->offset[4]);

  /* child[0]'s indices are correctly sifted */
  for (BfSize k = node->offset[0]; k < node->offset[1]; ++k)
    assert(point[perm[k]][0] <= split[0] && point[perm[k]][1] <= split[1]);

  /* child[1]'s indices are correctly sifted */
  for (BfSize k = node->offset[1]; k < node->offset[2]; ++k)
    assert(point[perm[k]][0] <= split[0] && point[perm[k]][1] > split[1]);

  /* child[2]'s indices are correctly sifted */
  for (BfSize k = node->offset[2]; k < node->offset[3]; ++k)
    assert(point[perm[k]][0] > split[0] && point[perm[k]][1] <= split[1]);

  /* child[3]'s indices are correctly sifted */
  for (BfSize k = node->offset[3]; k < node->offset[4]; ++k)
    assert(point[perm[k]][0] > split[0] && point[perm[k]][1] > split[1]);
#endif

  for (BfSize q = 0; q < 4; ++q) {
    BfSize childPermSize = node->offset[q + 1] - node->offset[q];

    /* TODO: right now, building quadtree to maximum depth possible */
    if (childPermSize <= 1) {
      node->child[q] = NULL;
      continue;
    }

    node->child[q] = malloc(sizeof(BfQuadtreeNode));
    node->child[q]->flags = childFlags[q];
    node->child[q]->parent = node;
    node->child[q]->bbox = childBbox[q];

    /* Compute the split for the `q`th child node */
    node->child[q]->split[0] = (childBbox[q].min[0] + childBbox[q].max[0])/2;
    node->child[q]->split[1] = (childBbox[q].min[1] + childBbox[q].max[1])/2;

    recInitQuadtreeNode(
      node->child[q],
      points, childBbox[q],
      node->offset[q], node->offset[q + 1], perm,
      currentDepth + 1);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {
    for (BfSize q = 0; q < 4; ++q)
      free(node->child[q]);
  }
}

void bfInitQuadtreeFromPoints(BfQuadtree *tree, BfPoints2 const *points, BfVectors2 const *unitNormals) {
  BEGIN_ERROR_HANDLING();

  tree->points = points;
  tree->unitNormals = unitNormals;

  BfSize numPoints = points->size;
  if (numPoints == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Initialize permutation to identity */
  tree->perm = bfPermIdentity(numPoints);
  HANDLE_ERROR();

  tree->root = malloc(sizeof(BfQuadtreeNode));
  if (tree->root == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  tree->root->flags = BF_QUADTREE_NODE_FLAG_ROOT;
  tree->root->parent = tree;

  /* Compute the bounding box for the entire quadtree */
  BfBbox2 bbox = bfGetPoints2BoundingBox(points);
  assert(!bfBbox2IsEmpty(&bbox));

  /* Rescale the bounding box so that it's square */
  BfReal w = bbox.max[0] - bbox.min[0];
  BfReal h = bbox.max[1] - bbox.min[1];
  if (w > h) {
    BfReal c = (bbox.min[1] + bbox.max[1])/2;
    bbox.min[1] = w*(bbox.min[1] - c)/h + c;
    bbox.max[1] = w*(bbox.max[1] - c)/h + c;
  } else {
    BfReal c = (bbox.min[0] + bbox.max[1])/2;
    bbox.min[0] = h*(bbox.min[0] - c)/h + c;
    bbox.max[0] = h*(bbox.max[0] - c)/h + c;
  }

  tree->root->bbox = bbox;

  /* Compute the split for the root node */
  tree->root->split[0] = (bbox.min[0] + bbox.max[0])/2;
  tree->root->split[1] = (bbox.min[1] + bbox.max[1])/2;

  recInitQuadtreeNode(tree->root, tree->points, bbox,
                      0, numPoints, tree->perm.index, 0);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFreeQuadtree(tree);
  }
}

void saveBoxesToTextFileRec(BfQuadtreeNode const *node, FILE *fp) {
  BfReal xmin = node->bbox.min[0];
  BfReal xmax = node->bbox.max[0];
  BfReal ymin = node->bbox.min[1];
  BfReal ymax = node->bbox.max[1];

  fprintf(fp, "%g %g %g %g\n", xmin, xmax, ymin, ymax);

  for (BfSize i = 0; i < 4; ++i)
    if (node->child[i])
      saveBoxesToTextFileRec(node->child[i], fp);
}

void bfQuadtreeSaveBoxesToTextFile(BfQuadtree const *tree, char const *path) {
  FILE *fp = fopen(path, "w");

  saveBoxesToTextFileRec(tree->root, fp);

  fclose(fp);
}

static void recFreeQuadtreeNode(BfQuadtreeNode *node) {
  for (BfSize i = 0; i < 4; ++i)
    if (node->child[i])
      recFreeQuadtreeNode(node->child[i]);

  free(node);
}

void bfFreeQuadtree(BfQuadtree *tree) {
  recFreeQuadtreeNode(tree->root);

  bfPermDeinit(&tree->perm);
}

void
bfGetQuadtreeNode(BfQuadtree const *tree, BfSize depth, BfSize nodeIndex,
                  BfQuadtreeNode **node)
{
  BEGIN_ERROR_HANDLING();

  BfSize nodesAtDepth = pow(4.0, depth);
  if (nodeIndex >= nodesAtDepth)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize i, r = nodeIndex;
  BfQuadtreeNode *currentNode = tree->root;
  do {
    nodesAtDepth /= 4;
    i = r/nodesAtDepth;
    r = r % nodesAtDepth;
    currentNode = currentNode->child[i];
    if (!currentNode)
      RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  } while (--depth > 0);
  *node = currentNode;

  END_ERROR_HANDLING() {}
}

BfCircle bfGetQuadtreeNodeBoundingCircle(BfQuadtreeNode const *node)
{
  BfReal const *min = node->bbox.min, *max = node->bbox.max;

  BfCircle circ;
  circ.r = hypot(max[0] - min[0], max[1] - min[1])/2;
  circ.center[0] = (min[0] + max[0])/2;
  circ.center[1] = (min[1] + max[1])/2;

  return circ;
}

BfSize bfQuadtreeNodeNumChildren(BfQuadtreeNode const *node) {
  BfSize numChildren = 0;
  for (BfSize i = 0; i < 4; ++i)
    numChildren += node->child[i] != NULL;
  return numChildren;
}

bool bfQuadtreeNodeIsLeaf(BfQuadtreeNode const *node) {
  return node->child[0] == NULL
      && node->child[1] == NULL
      && node->child[2] == NULL
      && node->child[3] == NULL;
}

BfSize bfQuadtreeNodeDepth(BfQuadtreeNode const *node) {
  return node->depth;
}

void
getMaxDepthBelowQuadtreeFunc(BfQuadtree *tree, BfQuadtreeNode *node,
                             BfSize *maxDepth)
{
  (void)tree;

  *maxDepth = node->depth > *maxDepth ? node->depth : *maxDepth;
}

BfSize bfGetMaxDepthBelowQuadtreeNode(BfQuadtreeNode const *node) {
  BfSize maxDepth = node->depth;

  bfMapQuadtreeNodes(
    bfGetQuadtreeFromNode(node),
    (BfQuadtreeNode *)node,
    BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    (BfQuadtreeFunc)getMaxDepthBelowQuadtreeFunc,
    &maxDepth);

  return maxDepth;
}

BfSize bfQuadtreeNodeNumPoints(BfQuadtreeNode const *node) {
  return node->offset[4] - node->offset[0];
}

BfQuadtree *bfGetQuadtreeFromNode(BfQuadtreeNode const *node) {
  while (node->flags & child_mask)
    node = node->parent;
  return node->parent;
}

/* Fill `points` with the points contained in `node`. They will be
 * added in quadtree order.
 *
 * If `tree == NULL`, then this function will retrieve the containing
 * `BfQuadtree` from `node`, which takes `O(log N)` time. */
void
bfGetQuadtreeNodePoints(BfQuadtree const *tree, BfQuadtreeNode const *node,
                        BfPoints2 *points)
{
  BEGIN_ERROR_HANDLING();

  if (tree == NULL)
    tree = bfGetQuadtreeFromNode(node);

  /* determine the number of points containined by `node` and find the
   * offset into `tree->perm` */
  BfSize numInds = node->offset[4] - node->offset[0];
  BfSize const *inds = &tree->perm.index[node->offset[0]];

  bfGetPointsByIndex(tree->points, numInds, inds, points);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

void bfQuadtreeNodeGetUnitNormals(BfQuadtree const *tree, BfQuadtreeNode const *node, BfVectors2 *unitNormals) {
  BEGIN_ERROR_HANDLING();

  if (tree == NULL)
    tree = bfGetQuadtreeFromNode(node);

  if (tree->unitNormals == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* determine the number of points containined by `node` and find the
   * offset into `tree->perm` */
  BfSize numInds = node->offset[4] - node->offset[0];
  BfSize const *inds = &tree->perm.index[node->offset[0]];

  bfGetVectorsByIndex(tree->unitNormals, numInds, inds, unitNormals);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

bool bfQuadtreeNodesAreSeparated(BfQuadtreeNode const *node1,
                                 BfQuadtreeNode const *node2) {
  BfCircle circ1 = bfGetQuadtreeNodeBoundingCircle(node1);
  BfCircle circ2 = bfGetQuadtreeNodeBoundingCircle(node2);

  BfReal R = bfPoint2Dist(circ1.center, circ2.center);

  return R > circ1.r + circ2.r + 1e1*BF_EPS_MACH;
}

static void clearDirtyBit(BfQuadtreeNode *node, void *arg) {
  (void)arg;
  node->flags &= ~BF_QUADTREE_NODE_FLAG_DIRTY;
}

static
void
fillWithLrLevelOrderNodePtrs(BfPtrArray *nodes, BfQuadtreeNode *current)
{
  BEGIN_ERROR_HANDLING();

  /* insert root node into array */
  current->flags |= BF_QUADTREE_NODE_FLAG_DIRTY;
  bfPtrArrayAppend(nodes, current);
  HANDLE_ERROR();

  /* insert all nodes beneath the initial node into nodes */
  BfSize i = 0;
  while (i < bfPtrArraySize(nodes)) {
    current = bfPtrArrayGet(nodes, i);
    for (BfSize j = 0; j < 4; ++j) {
      if (current->child[j] == NULL)
        continue;
      current->flags |= BF_QUADTREE_NODE_FLAG_DIRTY;
      bfPtrArrayAppend(nodes, current->child[j]);
      HANDLE_ERROR();
    }
    ++i;
  }

  END_ERROR_HANDLING() {}
}

static
void
findLevelOrderOffsets(BfPtrArray *nodes, BfSize *numLevels, BfSize **offsets)
{
  BEGIN_ERROR_HANDLING();

  /* start by computing the number of levels */

  BfQuadtreeNode const *prev = NULL;
  bfPtrArrayGetFirst(nodes, (BfPtr *)&prev);
  HANDLE_ERROR();
  BfSize minDepth = prev->depth;

  BfQuadtreeNode const *node = NULL;
  bfPtrArrayGetLast(nodes, (BfPtr *)&node);
  HANDLE_ERROR();
  BfSize maxDepth = node->depth;

  assert(maxDepth >= minDepth);
  *numLevels = maxDepth - minDepth + 1;

  /* allocate space for offsets */
  *offsets = malloc((*numLevels + 1)*sizeof(BfSize));
  if (*offsets == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* set sentinel values */
  (*offsets)[0] = 0;
  (*offsets)[*numLevels] = bfPtrArraySize(nodes);

  /* find each index in the level order where the depth increases */
  BfSize i = 1;
  for (BfSize j = 1; j < bfPtrArraySize(nodes); ++j) {
    node = bfPtrArrayGet(nodes, j);
    HANDLE_ERROR();
    if (node->depth != prev->depth)
      (*offsets)[i++] = j;
    prev = node;
  }

  END_ERROR_HANDLING() {
    free(*offsets);
  }
}

typedef struct {
  BfQuadtreeFunc func;
  BfQuadtree *tree;
  void *arg;
} WrappedArgs;

static
void
wrappedFunc(BfQuadtreeNode *node, WrappedArgs *args) {
  args->func(args->tree, node, args->arg);
}

static
void
mapQuadtreeNodesLrLevelOrder(BfQuadtree *tree, BfQuadtreeNode *node,
                             BfQuadtreeFunc func, void *arg)
{
  BEGIN_ERROR_HANDLING();

  /* initialize a queue of node pointers for the BFS */
  BfPtrArray queue;
  bfInitPtrArrayWithDefaultCapacity(&queue);
  HANDLE_ERROR();

  fillWithLrLevelOrderNodePtrs(&queue, node);
  HANDLE_ERROR();

  WrappedArgs wrappedArgs = {.func = func, .tree = tree, .arg = arg};
  bfMapPtrArray(&queue, (BfPtrFunc)wrappedFunc, &wrappedArgs);

  END_ERROR_HANDLING() {}

  bfMapPtrArray(&queue, (BfPtrFunc)clearDirtyBit, arg);
  bfFreePtrArray(&queue);
}

static
void
mapQuadtreeNodesLrReverseLevelOrder(BfQuadtree *tree, BfQuadtreeNode *node,
                                    BfQuadtreeFunc func, void *arg)
{
  BEGIN_ERROR_HANDLING();

  BfPtrArray queue;
  BfSize numLevels, *offsets = NULL;

  /* initialize a queue of node pointers for the BFS */
  bfInitPtrArrayWithDefaultCapacity(&queue);
  HANDLE_ERROR();

  /* fill the queue with the level order of the nodes */
  fillWithLrLevelOrderNodePtrs(&queue, node);
  HANDLE_ERROR();

  /* get the offsets to each level */
  findLevelOrderOffsets(&queue, &numLevels, &offsets);
  HANDLE_ERROR();

  /* now map each node in reverse level order, with each level
   * enumerated from left to right */
  for (BfSize j = numLevels; j > 0; --j) {
    for (BfSize k = offsets[j - 1]; k < offsets[j]; ++k) {
      node = bfPtrArrayGet(&queue, k);
      HANDLE_ERROR();
      func(tree, node, arg);
      HANDLE_ERROR();
    }
  }

  END_ERROR_HANDLING() {}

  bfMapPtrArray(&queue, (BfPtrFunc)clearDirtyBit, arg);
  bfFreePtrArray(&queue);
  free(offsets);
}

void
bfMapQuadtree(BfQuadtree *tree, enum BfTreeTraversals traversal,
              BfQuadtreeFunc func, void *arg)
{
  bfMapQuadtreeNodes(tree, tree->root, traversal, func, arg);
}

void
bfMapQuadtreeNodes(BfQuadtree *tree, BfQuadtreeNode *node,
                   enum BfTreeTraversals traversal,
                   BfQuadtreeFunc func, void *arg)
{
  BEGIN_ERROR_HANDLING();

  if (traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER) {
    mapQuadtreeNodesLrLevelOrder(tree, node, func, arg);
    HANDLE_ERROR();
    return;
  }

  if (traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    mapQuadtreeNodesLrReverseLevelOrder(tree, node, func, arg);
    HANDLE_ERROR();
    return;
  }

  RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  END_ERROR_HANDLING() {}
}

typedef struct LrLevelOrderInfo {
  BfSize currentLevel, numLevels, *offsets;
} LrLevelOrderInfo;

static
void
initLrLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                                  BfQuadtreeNode *node)
{
  BEGIN_ERROR_HANDLING();

  LrLevelOrderInfo *info = NULL;

  bfInitPtrArrayWithDefaultCapacity(&iter->nodes);
  HANDLE_ERROR();

  fillWithLrLevelOrderNodePtrs(&iter->nodes, node);
  HANDLE_ERROR();

  info = malloc(sizeof(LrLevelOrderInfo));
  if (info == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  iter->aux = info;

  findLevelOrderOffsets(&iter->nodes,&info->numLevels,&info->offsets);
  HANDLE_ERROR();

  info->currentLevel = 0;

  bfPtrArrayGetRangeView(
    &iter->nodes,
    info->offsets[info->currentLevel],
    info->offsets[info->currentLevel + 1],
    &iter->levelNodes);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFreePtrArray(&iter->nodes);
    free(info);
  }
}

static
void
initLrReverseLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                                         BfQuadtreeNode *node)
{
  BEGIN_ERROR_HANDLING();

  LrLevelOrderInfo *info = NULL;

  bfInitPtrArrayWithDefaultCapacity(&iter->nodes);
  HANDLE_ERROR();

  fillWithLrLevelOrderNodePtrs(&iter->nodes, node);
  HANDLE_ERROR();

  info = malloc(sizeof(LrLevelOrderInfo));
  if (info == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  iter->aux = info;

  findLevelOrderOffsets(&iter->nodes,&info->numLevels,&info->offsets);
  HANDLE_ERROR();

  info->currentLevel = info->numLevels - 1;

  bfPtrArrayGetRangeView(
    &iter->nodes,
    info->offsets[info->currentLevel],
    info->offsets[info->currentLevel + 1],
    &iter->levelNodes);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    bfFreePtrArray(&iter->nodes);
    free(info);
  }
}

BfQuadtreeLevelIter bfGetInvalidQuadtreeLevelIter() {
  return (BfQuadtreeLevelIter) {
    .traversal = BF_TREE_TRAVERSAL_UNKNOWN,
    .nodes = bfGetUninitializedPtrArray(),
    .levelNodes = bfGetUninitializedPtrArray(),
    .aux = NULL
  };
}

BfQuadtreeLevelIter
bfInitQuadtreeLevelIter(enum BfTreeTraversals traversal, BfQuadtreeNode *node)
{
  BEGIN_ERROR_HANDLING();

  BfQuadtreeLevelIter iter = {.traversal = traversal};

  if (traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER) {
    initLrLevelOrderQuadtreeLevelIter(&iter, node);
    HANDLE_ERROR();
    return iter;
  }

  if (traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    initLrReverseLevelOrderQuadtreeLevelIter(&iter, node);
    HANDLE_ERROR();
    return iter;
  }

  RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  END_ERROR_HANDLING() {}

  return bfGetInvalidQuadtreeLevelIter();
}

static
void
levelOrderQuadtreeLevelIterCurrentDepth(BfQuadtreeLevelIter const *iter,
                                        BfSize *depth)
{
  BEGIN_ERROR_HANDLING();

  LrLevelOrderInfo *info = iter->aux;
  BfSize i = info->offsets[info->currentLevel];
  BfQuadtreeNode *node = NULL;
  node = bfPtrArrayGet(&iter->nodes, i);
  HANDLE_ERROR();

  *depth = node->depth;

  END_ERROR_HANDLING() {}
}

BfSize
bfQuadtreeLevelIterCurrentDepth(BfQuadtreeLevelIter const *iter)
{
  BEGIN_ERROR_HANDLING();

  BfSize depth = BF_SIZE_BAD_VALUE;

  if (iter->traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER ||
      iter->traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    levelOrderQuadtreeLevelIterCurrentDepth(iter, &depth);
    HANDLE_ERROR();
  } else {
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  END_ERROR_HANDLING() {}

  return depth;
}

void
lrLevelOrderQuadtreeLevelIterNext(BfQuadtreeLevelIter *iter)
{
  LrLevelOrderInfo *info = iter->aux;

  if (info->currentLevel == info->numLevels)
    return;

  ++info->currentLevel;

  bfMakeEmptyPtrArrayView(&iter->levelNodes);
  bfPtrArrayGetRangeView(
    &iter->nodes,
    info->offsets[info->currentLevel],
    info->offsets[info->currentLevel + 1],
    &iter->levelNodes);
}

void
lrReverseLevelOrderQuadtreeLevelIterNext(BfQuadtreeLevelIter *iter)
{
  LrLevelOrderInfo *info = iter->aux;

  if (info->currentLevel == 0)
    return;

  --info->currentLevel;

  bfMakeEmptyPtrArrayView(&iter->levelNodes);
  bfPtrArrayGetRangeView(
    &iter->nodes,
    info->offsets[info->currentLevel],
    info->offsets[info->currentLevel + 1],
    &iter->levelNodes);
}

bool
bfQuadtreeLevelIterIsDone(BfQuadtreeLevelIter const *iter)
{
  return bfPtrArrayIsEmpty(&iter->levelNodes);
}

void
bfQuadtreeLevelIterNext(BfQuadtreeLevelIter *iter)
{
  if (iter->traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER) {
    lrLevelOrderQuadtreeLevelIterNext(iter);
    return;
  }

  if (iter->traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    lrReverseLevelOrderQuadtreeLevelIterNext(iter);
    return;
  }

  bfSetError(BF_ERROR_INVALID_ARGUMENTS);
}

static
void
freeLrLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter) {
  LrLevelOrderInfo *info = iter->aux;
  free(info->offsets);
  free(info);

  iter->aux = NULL;
}

void
bfFreeQuadtreeLevelIter(BfQuadtreeLevelIter *iter)
{
  bfFreePtrArray(&iter->nodes);

  if (iter->traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER ||
      iter->traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    freeLrLevelOrderQuadtreeLevelIter(iter);
    return;
  }

  bfSetError(BF_ERROR_INVALID_ARGUMENTS);
}

BfSize bfQuadtreeLevelIterNumPoints(BfQuadtreeLevelIter const *iter) {
  BfPtrArray const *levelNodes = &iter->levelNodes;
  BfSize numPoints = 0;
  for (BfSize i = 0; i < bfPtrArraySize(levelNodes); ++i) {
    BfQuadtreeNode *node = bfPtrArrayGet(levelNodes, i);
    numPoints += bfQuadtreeNodeNumPoints(node);
  }
  return numPoints;
}

#include "quadtree.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "error_macros.h"
#include "ptr_array.h"

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
  enum BfError error;
  bool erred = false;

  assert(i0 <= i1);

  BfPoint2 *point = points->data;
  BfReal const *split = node->split;
  node->depth = currentDepth;

  // set sentinel values
  node->offset[0] = i0;
#if BF_DEBUG
  node->offset[1] = BF_SIZE_BAD_VALUE;
  node->offset[2] = BF_SIZE_BAD_VALUE;
  node->offset[3] = BF_SIZE_BAD_VALUE;
#endif
  node->offset[4] = i1;

  BfSize i, j;

  // sift permutation indices for child[0] into place...
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

  // ... for child[1]...
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

  // ... for child[2]...
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

  // we don't have to do any sifting for the last child! nice.
  // ... but let's verify that things are as they should be.

  // compute the bounding boxes each child node
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
  /* sanity checks */
  assert(j == node->offset[4]);

  /* verify that child[0]'s indices are correctly sifted */
  for (BfSize k = node->offset[0]; k < node->offset[1]; ++k)
    assert(point[perm[k]][0] <= split[0] && point[perm[k]][1] <= split[1]);

  /* verify that child[1]'s indices are correctly sifted */
  for (BfSize k = node->offset[1]; k < node->offset[2]; ++k)
    assert(point[perm[k]][0] <= split[0] && point[perm[k]][1] > split[1]);

  /* verify that child[2]'s indices are correctly sifted */
  for (BfSize k = node->offset[2]; k < node->offset[3]; ++k)
    assert(point[perm[k]][0] > split[0] && point[perm[k]][1] <= split[1]);

  /* verify that child[3]'s indices are correctly sifted */
  for (BfSize k = node->offset[3]; k < node->offset[4]; ++k)
    assert(point[perm[k]][0] > split[0] && point[perm[k]][1] > split[1]);
#endif

  for (BfSize q = 0; q < 4; ++q) {
    BfSize childPermSize = node->offset[q + 1] - node->offset[q];

    // TODO: for now, we build the quadtree out to the maximum depth possible
    if (childPermSize <= 1) {
      node->child[q] = NULL;
      continue;
    }

    node->child[q] = malloc(sizeof(BfQuadtreeNode));
    node->child[q]->flags = childFlags[q];
    node->child[q]->parent = node;
    node->child[q]->bbox = childBbox[q];

    // compute the split for the `q`th child node
    node->child[q]->split[0] = (childBbox[q].min[0] + childBbox[q].max[0])/2;
    node->child[q]->split[1] = (childBbox[q].min[1] + childBbox[q].max[1])/2;

    recInitQuadtreeNode(
      node->child[q],
      points, childBbox[q],
      node->offset[q], node->offset[q + 1], perm,
      currentDepth + 1);
    HANDLE_ERROR();
  }

cleanup:
  if (erred)
    for (BfSize q = 0; q < 4; ++q)
      free(node->child[q]);
}

void
bfInitQuadtreeFromPoints(BfQuadtree *tree, BfPoints2 const *points)
{
  enum BfError error;
  bool erred = false;

  tree->points = points;

  BfSize numPoints = points->size;
  if (numPoints == 0)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  // initialize permutation to identity
  tree->perm = malloc(numPoints*sizeof(BfSize));
  if (tree->perm == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  for (BfSize i = 0; i < numPoints; ++i)
    tree->perm[i] = i;

  tree->root = malloc(sizeof(BfQuadtreeNode));
  if (tree->root == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);
  tree->root->flags = BF_QUADTREE_NODE_FLAG_ROOT;
  tree->root->parent = tree;

  // compute the bounding box for the entire quadtree
  BfBbox2 bbox = bfGetPoints2BoundingBox(points);
  assert(!bfBbox2IsEmpty(&bbox));

  // rescale the bounding box so that it's square
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

  // compute the split for the root node
  tree->root->split[0] = (bbox.min[0] + bbox.max[0])/2;
  tree->root->split[1] = (bbox.min[1] + bbox.max[1])/2;

  recInitQuadtreeNode(tree->root, tree->points, bbox,
                      0, numPoints, tree->perm, 0);
  HANDLE_ERROR();

cleanup:
  if (erred)
    bfFreeQuadtree(tree);
}

static void recFreeQuadtreeNode(BfQuadtreeNode *node) {
  for (BfSize i = 0; i < 4; ++i)
    if (node->child[i])
      recFreeQuadtreeNode(node->child[i]);

  free(node);
}

void bfFreeQuadtree(BfQuadtree *tree) {
  recFreeQuadtreeNode(tree->root);

  free(tree->perm);
}

void
bfGetQuadtreeNode(BfQuadtree const *tree, BfSize depth, BfSize nodeIndex,
                  BfQuadtreeNode **node)
{
  BfSize nodesAtDepth = pow(4.0, depth);
  if (nodeIndex >= nodesAtDepth) {
    bfSetError(BF_ERROR_INVALID_ARGUMENTS);
    return;
  }

  BfSize i, r = nodeIndex;
  BfQuadtreeNode *currentNode = tree->root;
  do {
    nodesAtDepth /= 4;
    i = r/nodesAtDepth;
    r = r % nodesAtDepth;
    currentNode = currentNode->child[i];
    if (!currentNode) {
      bfSetError(BF_ERROR_INVALID_ARGUMENTS);
      return;
    }
  } while (--depth > 0);
  *node = currentNode;
}

BfCircle2 bfGetQuadtreeNodeBoundingCircle(BfQuadtreeNode const *node)
{
  BfReal const *min = node->bbox.min, *max = node->bbox.max;

  BfCircle2 circ;
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
  enum BfError error;
  bool erred = false;

  if (tree == NULL)
    tree = bfGetQuadtreeFromNode(node);

  /* determine the number of points containined by `node` and find the
   * offset into `tree->perm` */
  BfSize numInds = node->offset[4] - node->offset[0];
  BfSize const *inds = &tree->perm[node->offset[0]];

  bfGetPointsByIndex(tree->points, numInds, inds, points);
  HANDLE_ERROR();

cleanup:
  if (erred)
    bfFreePoints2(points);
}

static void clearDirtyBit(BfQuadtreeNode *node, void *arg) {
  (void)arg;
  node->flags &= ~BF_QUADTREE_NODE_FLAG_DIRTY;
}

static
void
fillWithLrLevelOrderNodePtrs(BfPtrArray *nodes, BfQuadtreeNode *current)
{
  enum BfError error;

  /* insert root node into array */
  current->flags |= BF_QUADTREE_NODE_FLAG_DIRTY;
  bfPtrArrayAppend(nodes, current);
  error = bfGetError();
  if (error) {
    bfSetError(error);
    return;
  }

  /* insert all nodes beneath the initial node into nodes */
  BfSize i = 0;
  while (i < bfPtrArraySize(nodes)) {
    bfPtrArrayGet(nodes, i, (BfPtr *)&current);
    for (BfSize j = 0; j < 4; ++j) {
      if (current->child[j] == NULL)
        continue;
      current->flags |= BF_QUADTREE_NODE_FLAG_DIRTY;
      bfPtrArrayAppend(nodes, current->child[j]);
      error = bfGetError();
      if (error) {
        bfSetError(error);
        return;
      }
    }
    ++i;
  }
}

static
void
findLevelOrderOffsets(BfPtrArray *nodes, BfSize *numLevels, BfSize **offsets)
{
  enum BfError error;
  bool erred = false;

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
    bfPtrArrayGet(nodes, j, (BfPtr *)&node);
    HANDLE_ERROR();
    if (node->depth != prev->depth)
      (*offsets)[i++] = j;
    prev = node;
  }

cleanup:
  if (erred)
    free(*offsets);
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
  enum BfError error;
  bool erred = false;

  (void)erred;

  /* initialize a queue of node pointers for the BFS */
  BfPtrArray queue;
  bfInitPtrArrayWithDefaultCapacity(&queue);
  HANDLE_ERROR();

  fillWithLrLevelOrderNodePtrs(&queue, node);
  HANDLE_ERROR();

  WrappedArgs wrappedArgs = {.func = func, .tree = tree, .arg = arg};
  bfMapPtrArray(&queue, (BfPtrFunc)wrappedFunc, &wrappedArgs);

cleanup:
  bfMapPtrArray(&queue, (BfPtrFunc)clearDirtyBit, arg);
  bfFreePtrArray(&queue);
}

static
void
mapQuadtreeNodesLrReverseLevelOrder(BfQuadtree *tree, BfQuadtreeNode *node,
                                    BfQuadtreeFunc func, void *arg)
{
  enum BfError error;
  bool erred = false;

  (void)erred;

  /* initialize a queue of node pointers for the BFS */
  BfPtrArray queue;
  bfInitPtrArrayWithDefaultCapacity(&queue);
  HANDLE_ERROR();

  /* fill the queue with the level order of the nodes */
  fillWithLrLevelOrderNodePtrs(&queue, node);
  HANDLE_ERROR();

  /* get the offsets to each level */
  BfSize numLevels, *offsets = NULL;
  findLevelOrderOffsets(&queue, &numLevels, &offsets);
  HANDLE_ERROR();

  /* now map each node in reverse level order, with each level
   * enumerated from left to right */
  for (BfSize j = numLevels; j > 0; --j) {
    for (BfSize k = offsets[j - 1]; k < offsets[j]; ++k) {
      bfPtrArrayGet(&queue, k, (BfPtr *)&node);
      HANDLE_ERROR();
      func(tree, node, arg);
      HANDLE_ERROR();
    }
  }

cleanup:
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
  if (traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER) {
    mapQuadtreeNodesLrLevelOrder(tree, node, func, arg);
    return;
  }

  if (traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    mapQuadtreeNodesLrReverseLevelOrder(tree, node, func, arg);
    return;
  }

  bfSetError(BF_ERROR_INVALID_ARGUMENTS);
}

typedef struct LrLevelOrderInfo {
  BfSize currentLevel, numLevels, *offsets;
} LrLevelOrderInfo;

static
void
initLrLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                                  BfQuadtreeNode *node)
{
  enum BfError error;
  bool erred = false;

  bfInitPtrArrayWithDefaultCapacity(&iter->nodes);
  HANDLE_ERROR();

  fillWithLrLevelOrderNodePtrs(&iter->nodes, node);
  HANDLE_ERROR();

  LrLevelOrderInfo *info = malloc(sizeof(LrLevelOrderInfo));
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

cleanup:
  if (erred) {
    bfFreePtrArray(&iter->nodes);
    free(info);
  }
}

static
void
initLrReverseLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                                         BfQuadtreeNode *node)
{
  enum BfError error;
  bool erred = false;

  bfInitPtrArrayWithDefaultCapacity(&iter->nodes);
  HANDLE_ERROR();

  fillWithLrLevelOrderNodePtrs(&iter->nodes, node);
  HANDLE_ERROR();

  LrLevelOrderInfo *info = malloc(sizeof(LrLevelOrderInfo));
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

cleanup:
  if (erred) {
    bfFreePtrArray(&iter->nodes);
    free(info);
  }
}

void
bfInitQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                        enum BfTreeTraversals traversal, BfQuadtreeNode *node)
{
  enum BfError error;

  iter->traversal = traversal;

  if (traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER) {
    initLrLevelOrderQuadtreeLevelIter(iter, node);
    error = bfGetError();
    if (error)
      bfSetError(error);
    return;
  }

  if (traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    initLrReverseLevelOrderQuadtreeLevelIter(iter, node);
    error = bfGetError();
    if (error)
      bfSetError(error);
    return;
  }

  bfSetError(BF_ERROR_INVALID_ARGUMENTS);
}

static
void
levelOrderQuadtreeLevelIterCurrentDepth(BfQuadtreeLevelIter const *iter,
                                        BfSize *depth)
{
  LrLevelOrderInfo *info = iter->aux;

  BfSize i = info->offsets[info->currentLevel];

  BfQuadtreeNode *node = NULL;
  bfPtrArrayGet(&iter->nodes, i, (BfPtr *)&node);
  enum BfError error = bfGetError();
  if (error) {
    bfSetError(error);
    return;
  }

  *depth = node->depth;
}

BfSize
bfQuadtreeLevelIterCurrentDepth(BfQuadtreeLevelIter const *iter)
{
  enum BfError error;

  BfSize depth;

  if (iter->traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER ||
      iter->traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    levelOrderQuadtreeLevelIterCurrentDepth(iter, &depth);
    error = bfGetError();
    if (error) {
      bfSetError(error);
      return BF_SIZE_BAD_VALUE;
    }
    return depth;
  }

  bfSetError(BF_ERROR_INVALID_ARGUMENTS);

  return BF_SIZE_BAD_VALUE;
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

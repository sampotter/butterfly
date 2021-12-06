#include "quadtree.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#include "ptr_array.h"

#define SWAP(x, y) do {                         \
    __typeof__(x) tmp = x;                      \
    x = y;                                      \
    y = tmp;                                    \
  } while (0)

enum BfQuadtreeNodeFlags child_flags[4] = {
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

size_t child_flag_to_index[] = {
  [BF_QUADTREE_NODE_FLAG_CHILD_0] = 0,
  [BF_QUADTREE_NODE_FLAG_CHILD_1] = 1,
  [BF_QUADTREE_NODE_FLAG_CHILD_2] = 2,
  [BF_QUADTREE_NODE_FLAG_CHILD_3] = 3
};

static
enum BfError
recInitQuadtreeFromPoints(BfQuadtreeNode *node,
                          BfPoint2 const *points, BfBbox2 bbox,
                          size_t perm_size, size_t *perm,
                          BfSize current_depth)
{
  enum BfError error = BF_ERROR_NO_ERROR;
  size_t i, j;

  BfReal const *split = node->split;

  node->depth = current_depth;

  // set sentinel values
  node->offset[0] = 0;
  node->offset[4] = perm_size;

  // sift permutation indices for child[0] into place...
  i = node->offset[0];
  while (i < perm_size &&
         (points[perm[i]][0] <= split[0] && points[perm[i]][1] <= split[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] <= split[0] && points[perm[j]][1] <= split[1] &&
        (points[perm[i]][0] > split[0] || points[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[1] = i;

#ifdef BF_DEBUG
  // (verify that child[0]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] <= split[0] && points[perm[k]][1] <= split[1];
    assert(node->offset[0] <= k && k < node->offset[1] ?
           in_quadrant : !in_quadrant);
  }
#endif

  // ... for child[1]...
  i = node->offset[1];
  while (i < perm_size &&
         (points[perm[i]][0] <= split[0] && points[perm[i]][1] > split[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] <= split[0] && points[perm[j]][1] > split[1] &&
        (points[perm[i]][0] > split[0] || points[perm[i]][1] <= split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[2] = i;

#ifdef BF_DEBUG
  // (verify that child[1]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] <= split[0] && points[perm[k]][1] > split[1];
    assert(node->offset[1] <= k && k < node->offset[2] ?
           in_quadrant : !in_quadrant);
  }
#endif

  // ... for child[2]...
  i = node->offset[2];
  while (i < perm_size &&
         (points[perm[i]][0] > split[0] && points[perm[i]][1] <= split[1]))
    ++i;
  j = i == perm_size ? i : i + 1;
  while (j < perm_size) {
    if (points[perm[j]][0] > split[0] && points[perm[j]][1] <= split[1] &&
        (points[perm[i]][0] <= split[0] || points[perm[i]][1] > split[1])) {
      SWAP(perm[i], perm[j]);
      ++i;
    }
    ++j;
  }
  node->offset[3] = i;

#ifdef BF_DEBUG
  // (verify that child[2]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] > split[0] && points[perm[k]][1] <= split[1];
    assert(node->offset[2] <= k && k < node->offset[3] ?
           in_quadrant : !in_quadrant);
  }
#endif

  // we don't have to do any sifting for the last child! nice.
  // ... but let's verify that things are as they should be.

#ifdef BF_DEBUG
  assert(j == perm_size);
  assert(j == node->offset[4]);
#endif

#ifdef BF_DEBUG
  // (verify that child[2]'s indices are correctly sifted)
  for (size_t k = 0; k < perm_size; ++k) {
    bool in_quadrant =
      points[perm[k]][0] > split[0] && points[perm[k]][1] > split[1];
    assert(node->offset[3] <= k && k < node->offset[4] ?
           in_quadrant : !in_quadrant);
  }
#endif

#ifdef BF_DEBUG
  // (verify that child[3]'s indices are correctly sifted)
  for (size_t i = node->offset[3], j; i < perm_size; ++i) {
    j = perm[i];
    assert(points[j][0] > split[0] && points[j][1] > split[1]);
  }
#endif

  // compute the bounding boxes each child node
  BfBbox2 child_bbox[4] = {
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

  for (size_t q = 0, perm_size; q < 4; ++q) {
    perm_size = node->offset[q + 1] - node->offset[q];

    // TODO: for now, we build the quadtree out to the maximum depth possible
    if (perm_size <= 1) {
      node->child[q] = NULL;
      continue;
    }

    node->child[q] = malloc(sizeof(BfQuadtreeNode));
    node->child[q]->flags = child_flags[q];
    node->child[q]->parent = node;
    node->child[q]->bbox = child_bbox[q];

    // compute the split for the `q`th child node
    node->child[q]->split[0] = (child_bbox[q].min[0] + child_bbox[q].max[0])/2;
    node->child[q]->split[1] = (child_bbox[q].min[1] + child_bbox[q].max[1])/2;

    error |= recInitQuadtreeFromPoints(
      node->child[q], points, child_bbox[q], perm_size, &perm[node->offset[q]],
      current_depth + 1);

    if (error)
      break;
  }

  if (error)
    for (size_t q = 0; q < 4; ++q)
      free(node->child[q]);

  return error;
}

enum BfError
bfInitQuadtreeFromPoints(BfQuadtree *tree,
                         size_t num_points, BfPoint2 const *points)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  tree->num_points = num_points;
  tree->points = points;

  // initialize permutation to identity
  tree->perm = malloc(num_points*sizeof(size_t));
  for (size_t i = 0; i < num_points; ++i)
    tree->perm[i] = i;

  tree->root = malloc(sizeof(BfQuadtreeNode));
  tree->root->flags = BF_QUADTREE_NODE_FLAG_ROOT;
  tree->root->parent = tree;

  // compute the bounding box for the entire quadtree
  BfBbox2 bbox = {
    .min = {INFINITY, INFINITY},
    .max = {-INFINITY, -INFINITY}
  };
  for (size_t i = 0; i < num_points; ++i) {
    bbox.min[0] = fmin(bbox.min[0], points[i][0]);
    bbox.max[0] = fmax(bbox.max[0], points[i][0]);
    bbox.min[1] = fmin(bbox.min[1], points[i][1]);
    bbox.max[1] = fmax(bbox.max[1], points[i][1]);
  }

  assert(bbox.min[0] < bbox.max[0]);
  assert(bbox.min[1] < bbox.max[1]);

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

  error |= recInitQuadtreeFromPoints(tree->root, tree->points, bbox,
                                     tree->num_points, tree->perm, 0);

  return error;
}

enum BfError
bfGetQuadtreeNode(BfQuadtree const *tree, size_t depth, size_t node_index,
                  BfQuadtreeNode **node)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  size_t nodes_at_depth = pow(4.0, depth);
  if (node_index >= nodes_at_depth)
    return error | BF_ERROR_INVALID_ARGUMENTS;

  size_t i, r = node_index;

  BfQuadtreeNode *current_node = tree->root;
  do {
    nodes_at_depth /= 4;
    i = r/nodes_at_depth;
    r = r % nodes_at_depth;
    current_node = current_node->child[i];
    if (!current_node)
      return error | BF_ERROR_INVALID_ARGUMENTS | BF_ERROR_RUNTIME_ERROR;
  } while (--depth > 0);

  *node = current_node;

  return error;
}

enum BfError
bfGetQuadtreeNodeIndices(BfQuadtreeNode const *node,
                         size_t *num_indices, size_t **indices)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  assert(node->offset[0] == 0); // I think?

  *num_indices = node->offset[4];

  BfQuadtreeNode const *parent_node;
  size_t q, offset = 0;

  while (node->flags & child_mask) {
    q = child_flag_to_index[node->flags];
    parent_node = node->parent;
    offset += parent_node->offset[q];
    node = parent_node;
  }

  BfQuadtree const *tree = node->parent;

  *indices = tree->perm + offset;

  return error;
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

bool bfQuadtreeNodeIsLeaf(BfQuadtreeNode const *node) {
  return node->child[0] == NULL
      && node->child[1] == NULL
      && node->child[2] == NULL
      && node->child[3] == NULL;
}

BfSize bfQuadtreeNodeDepth(BfQuadtreeNode const *node) {
  return node->depth;
}

static
enum BfError
findMaxDepth(BfQuadtreeNode const *node, void *arg) {
  BfSize *max_depth = arg;
  if (node->depth > *max_depth)
    *max_depth = node->depth;
  return BF_ERROR_NO_ERROR;
}

BfSize bfGetMaxDepthBelowQuadtreeNode(BfQuadtreeNode const *node) {
  BfSize max_depth = 0;
  bfMapQuadtreeNodes(
    (BfQuadtreeNode *)node, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
    findMaxDepth, &max_depth);
  return max_depth;
}

BfSize bfQuadtreeNodeNumPoints(BfQuadtreeNode const *node) {
  return node->offset[4];
}

BfQuadtree *bfGetQuadtreeFromNode(BfQuadtreeNode const *node) {
  while (node->flags & child_mask)
    node = node->parent;
  return node->parent;
}

/* Fill `X` with the points contained in `node`. They will be added in
 * "quadtree order", so that they are permuted to match the points
 * contained by `node`'s children. E.g., points `X[offset[0]]`, ...,
 * `X[offset[1] - 1]` will correspond to the points returned by
 * `bfGetQuadtreeFromNode(node->child[0])`, and in the same order. */
enum BfError
bfGetQuadtreeNodePoints(BfQuadtreeNode const *node, BfMat *X)
{
  size_t num_points;
  size_t *node_indices;
  bfGetQuadtreeNodeIndices(node, &num_points, &node_indices);

  if (X->shape[0] != num_points)
    return BF_ERROR_INVALID_ARGUMENTS;

  BfQuadtree const *tree = bfGetQuadtreeFromNode(node);

  for (BfSize i = 0; i < num_points; ++i) {
    enum BfError error = bfSetMatRow(X, i, tree->points[node_indices[i]]);
    if (error)
      return error;
  }

  return BF_ERROR_NO_ERROR;
}

static enum BfError clearDirtyBit(BfQuadtreeNode *node, void *arg) {
  (void)arg;
  node->flags &= ~BF_QUADTREE_NODE_FLAG_DIRTY;
  return BF_ERROR_NO_ERROR;
}

static
enum BfError
fillWithLrLevelOrderNodePtrs(BfPtrArray *nodes, BfQuadtreeNode *current)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  /* insert root node into array */
  current->flags |= BF_QUADTREE_NODE_FLAG_DIRTY;
  error = bfPtrArrayAppend(nodes, current);
  if (error)
    return error;

  /* insert all nodes beneath the initial node into nodes */
  BfSize i = 0;
  while (i < bfPtrArraySize(nodes)) {
    bfPtrArrayGet(nodes, i, (BfPtr *)&current);
    for (BfSize j = 0; j < 4; ++j) {
      if (current->child[j] == NULL)
        continue;
      current->flags |= BF_QUADTREE_NODE_FLAG_DIRTY;
      error = bfPtrArrayAppend(nodes, current->child[j]);
      if (error)
        return error;
    }
    ++i;
  }

  return error;
}

static
enum BfError
findLevelOrderOffsets(BfPtrArray *nodes, BfSize *num_levels, BfSize **offsets)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  BfQuadtreeNode *prev, *node;

  /* start by computing the number of levels */

  BfSize min_depth;
  bfPtrArrayGetFirst(nodes, (BfPtr *)&prev);
  min_depth = prev->depth;

  BfSize max_depth;
  bfPtrArrayGetLast(nodes, (BfPtr *)&node);
  max_depth = node->depth;

  assert(max_depth >= min_depth);
  *num_levels = max_depth - min_depth + 1;

  /* allocate space for offsets */
  *offsets = malloc((*num_levels + 1)*sizeof(BfSize));
  if (*offsets == NULL)
    return BF_ERROR_MEMORY_ERROR;

  /* set sentinel values */
  (*offsets)[0] = 0;
  (*offsets)[*num_levels] = bfPtrArraySize(nodes);

  /* find each index in the level order where the depth increases */
  BfSize i = 1;
  for (BfSize j = 1; j < bfPtrArraySize(nodes); ++j) {
    bfPtrArrayGet(nodes, j, (BfPtr *)&node);
    if (node->depth != prev->depth)
      (*offsets)[i++] = j;
    prev = node;
  }

  return error;
}

static
enum BfError
mapQuadtreeNodesLrLevelOrder(BfQuadtreeNode *node,
                             BfQuadtreeNodeFunc func, void *arg)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  /* initialize a queue of node pointers for the BFS */
  BfPtrArray queue;
  error = bfInitPtrArrayWithDefaultCapacity(&queue);
  if (error)
    goto cleanup;

  error = fillWithLrLevelOrderNodePtrs(&queue, node);
  if (error)
    goto cleanup;

  error = bfMapPtrArray(&queue, (BfPtrFunc)func, arg);

cleanup:
  bfMapPtrArray(&queue, (BfPtrFunc)clearDirtyBit, arg);
  bfFreePtrArray(&queue);

  return error;
}

static
enum BfError
mapQuadtreeNodesLrReverseLevelOrder(BfQuadtreeNode *node,
                                    BfQuadtreeNodeFunc func, void *arg)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  /* initialize a queue of node pointers for the BFS */
  BfPtrArray queue;
  error = bfInitPtrArrayWithDefaultCapacity(&queue);
  if (error)
    goto cleanup;

  /* fill the queue with the level order of the nodes */
  fillWithLrLevelOrderNodePtrs(&queue, node);

  /* get the offsets to each level */
  BfSize num_levels, *offsets;
  error = findLevelOrderOffsets(&queue, &num_levels, &offsets);
  if (error)
    goto cleanup;

  /* now map each node in reverse level order, with each level
   * enumerated from left to right */
  for (BfSize j = num_levels; j > 0; --j) {
    for (BfSize k = offsets[j - 1]; k < offsets[j]; ++k) {
      bfPtrArrayGet(&queue, k, (BfPtr *)&node);
      error = func(node, arg);
      if (error)
        goto cleanup;
    }
  }

cleanup:
  bfMapPtrArray(&queue, (BfPtrFunc)clearDirtyBit, arg);
  bfFreePtrArray(&queue);
  free(offsets);

  return error;
}

enum BfError
bfMapQuadtreeNodes(BfQuadtreeNode *node, enum BfTreeTraversals traversal,
                   BfQuadtreeNodeFunc func, void *arg)
{
  switch (traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
    return mapQuadtreeNodesLrLevelOrder(node, func, arg);
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER:
    return mapQuadtreeNodesLrReverseLevelOrder(node, func, arg);
  default:
    return BF_ERROR_INVALID_ARGUMENTS;
  }
}

enum BfError
bfMapQuadtreeNodeLeaves(BfQuadtreeNode const *node,
                        BfQuadtreeNodeFunc func, void *arg)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  if (bfQuadtreeNodeIsLeaf(node)) {
    error |= func(node, arg);
    return error;
  }

  for (BfSize i = 0; i < 4; ++i) {
    if (node->child[i] == NULL)
      continue;

    error = bfMapQuadtreeNodeLeaves(node->child[i], func, arg);
    if (error)
      return error;
  }

  return error;
}

typedef struct LrLevelOrderInfo {
  BfSize current_level, num_levels, *offsets;
} LrLevelOrderInfo;

static
enum BfError
initLrLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                                  BfQuadtreeNode *node)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  error = bfInitPtrArrayWithDefaultCapacity(&iter->nodes);
  if (error)
    return error;

  error = fillWithLrLevelOrderNodePtrs(&iter->nodes, node);
  if (error)
    return error;

  LrLevelOrderInfo *info = malloc(sizeof(LrLevelOrderInfo));
  iter->aux = info;

  error = findLevelOrderOffsets(&iter->nodes,&info->num_levels,&info->offsets);
  if (error)
    return error;

  info->current_level = 0;

  error = bfPtrArrayGetRangeView(
    &iter->nodes,
    info->offsets[info->current_level],
    info->offsets[info->current_level + 1],
    &iter->level_nodes);

  return error;
}

static
enum BfError
initLrReverseLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                                         BfQuadtreeNode *node)
{
  enum BfError error = BF_ERROR_NO_ERROR;

  error = bfInitPtrArrayWithDefaultCapacity(&iter->nodes);
  if (error)
    return error;

  error = fillWithLrLevelOrderNodePtrs(&iter->nodes, node);
  if (error)
    return error;

  LrLevelOrderInfo *info = malloc(sizeof(LrLevelOrderInfo));
  iter->aux = info;

  error = findLevelOrderOffsets(&iter->nodes,&info->num_levels,&info->offsets);
  if (error)
    return error;

  info->current_level = info->num_levels;

  error = bfPtrArrayGetRangeView(
    &iter->nodes,
    info->offsets[info->current_level - 1],
    info->offsets[info->current_level],
    &iter->level_nodes);

  return error;
}

enum BfError
bfInitQuadtreeLevelIter(BfQuadtreeLevelIter *iter,
                        enum BfTreeTraversals traversal, BfQuadtreeNode *node)
{
  iter->traversal = traversal;

  switch (iter->traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
    return initLrLevelOrderQuadtreeLevelIter(iter, node);
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER:
    return initLrReverseLevelOrderQuadtreeLevelIter(iter, node);
  default:
    return BF_ERROR_INVALID_ARGUMENTS;
  }
}

static
enum BfError
levelOrderQuadtreeLevelIterCurrentDepth(BfQuadtreeLevelIter const *iter,
                                        BfSize *depth)
{
  LrLevelOrderInfo *info = iter->aux;

  BfSize i = info->offsets[info->current_level];

  BfQuadtreeNode *node;
  bfPtrArrayGet(&iter->nodes, i, (BfPtr *)&node);

  *depth = node->depth;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfQuadtreeLevelIterCurrentDepth(BfQuadtreeLevelIter const *iter, BfSize *depth)
{
  switch (iter->traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER:
    return levelOrderQuadtreeLevelIterCurrentDepth(iter, depth);
  default:
    return BF_ERROR_INVALID_ARGUMENTS;
  }
}

enum BfError
lrLevelOrderQuadtreeLevelIterNext(BfQuadtreeLevelIter *iter)
{
  LrLevelOrderInfo *info = iter->aux;

  ++info->current_level;

  bfMakeEmptyPtrArrayView(&iter->level_nodes);

  return info->current_level == info->num_levels ?
    BF_ERROR_NO_ERROR :
    bfPtrArrayGetRangeView(
      &iter->nodes,
      info->offsets[info->current_level],
      info->offsets[info->current_level + 1],
      &iter->level_nodes);
}

enum BfError
lrReverseLevelOrderQuadtreeLevelIterNext(BfQuadtreeLevelIter *iter)
{
  LrLevelOrderInfo *info = iter->aux;

  bfMakeEmptyPtrArrayView(&iter->level_nodes);

  enum BfError error = info->current_level == 0 ?
    BF_ERROR_NO_ERROR :
    bfPtrArrayGetRangeView(
      &iter->nodes,
      info->offsets[info->current_level - 1],
      info->offsets[info->current_level],
      &iter->level_nodes);

  if (info->current_level > 0)
    --info->current_level;

  return error;
}

bool
bfQuadtreeLevelIterIsDone(BfQuadtreeLevelIter const *iter)
{
  return bfPtrArrayIsEmpty(&iter->level_nodes);
}

enum BfError
bfQuadtreeLevelIterNext(BfQuadtreeLevelIter *iter)
{
  switch (iter->traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
    return lrLevelOrderQuadtreeLevelIterNext(iter);
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER:
    return lrReverseLevelOrderQuadtreeLevelIterNext(iter);
  default:
    return BF_ERROR_INVALID_ARGUMENTS;
  }
}

static
enum BfError
freeLrLevelOrderQuadtreeLevelIter(BfQuadtreeLevelIter *iter) {
  free((LrLevelOrderInfo *)iter->aux);

  iter->aux = NULL;

  return BF_ERROR_NO_ERROR;
}

enum BfError
bfFreeQuadtreeLevelIter(BfQuadtreeLevelIter *iter)
{
  bfFreePtrArray(&iter->nodes);

  switch (iter->traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER:
    return freeLrLevelOrderQuadtreeLevelIter(iter);
  default:
    return BF_ERROR_INVALID_ARGUMENTS;
  }
}

#include <bf/tree_level_iter.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/tree_node.h>

static
void
fillWithLrLevelOrderNodePtrs(BfPtrArray *nodes, BfTreeNode *current)
{
  BEGIN_ERROR_HANDLING();

  /* Insert root node into array */
  bfPtrArrayAppend(nodes, current);
  HANDLE_ERROR();

  /* Insert all nodes beneath the initial node into nodes */
  BfSize i = 0;
  while (i < bfPtrArraySize(nodes)) {
    current = bfPtrArrayGet(nodes, i);
    for (BfSize j = 0; j < 4; ++j) {
      if (current->child[j] == NULL)
        continue;
      bfPtrArrayAppend(nodes, current->child[j]);
      HANDLE_ERROR();
    }
    ++i;
  }

#if BF_DEBUG
  BfTreeNode const *prevNode = bfPtrArrayGet(nodes, 0);
  for (BfSize i = 1; i < bfPtrArraySize(nodes); ++i) {
    BfTreeNode const *node = bfPtrArrayGet(nodes, i);
    assert(prevNode->depth <= node->depth);
    prevNode = node;
  }
#endif

  END_ERROR_HANDLING() {}
}

static
void
findLevelOrderOffsets(BfPtrArray *nodes, BfSize *numLevels, BfSize **offsets)
{
  BEGIN_ERROR_HANDLING();

  /* start by computing the number of levels */

  BfTreeNode const *prev = NULL;
  bfPtrArrayGetFirst(nodes, (BfPtr *)&prev);
  HANDLE_ERROR();
  BfSize minDepth = prev->depth;

  BfTreeNode const *node = NULL;
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

typedef struct LrLevelOrderInfo {
  BfSize currentLevel, numLevels, *offsets;
} LrLevelOrderInfo;

static void treeLevelIterInit_lrLevelOrder(BfTreeLevelIter *iter, BfTreeNode *node) {
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
    bfPtrArrayDeinit(&iter->nodes);
    free(info);
  }
}

static
void
treeLevelIterInit_lrReverseLevelOrder(BfTreeLevelIter *iter, BfTreeNode *node) {
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
    bfPtrArrayDeinit(&iter->nodes);
    free(info);
  }
}

/** Implementation: TreeLevelIter */

void bfTreeLevelIterInit(BfTreeLevelIter *iter, BfTreeTraversal traversal, BfTreeNode *node) {
  BEGIN_ERROR_HANDLING();

  iter->traversal = traversal;

  switch (traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
    treeLevelIterInit_lrLevelOrder(iter, node);
    HANDLE_ERROR();
    break;
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER:
    treeLevelIterInit_lrReverseLevelOrder(iter, node);
    HANDLE_ERROR();
    break;
  default:
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  END_ERROR_HANDLING() {}
}

static void currentDepth_levelOrder(BfTreeLevelIter const *iter, BfSize *depth) {
  BEGIN_ERROR_HANDLING();

  LrLevelOrderInfo *info = iter->aux;
  BfSize i = info->offsets[info->currentLevel];
  BfTreeNode *node = NULL;
  node = bfPtrArrayGet(&iter->nodes, i);
  HANDLE_ERROR();

  *depth = node->depth;

  END_ERROR_HANDLING() {}
}

BfSize bfTreeLevelIterCurrentDepth(BfTreeLevelIter const *iter) {
  BEGIN_ERROR_HANDLING();

  BfSize depth = BF_SIZE_BAD_VALUE;

  if (iter->traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER ||
      iter->traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    currentDepth_levelOrder(iter, &depth);
    HANDLE_ERROR();
  } else {
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  END_ERROR_HANDLING() {}

  return depth;
}

bool bfTreeLevelIterIsDone(BfTreeLevelIter const *iter) {
  LrLevelOrderInfo *info = iter->aux;
  assert(info->currentLevel <= info->numLevels);
  return info->currentLevel == info->numLevels;
}

static void next_lrLevelOrder(BfTreeLevelIter *iter) {
  LrLevelOrderInfo *info = iter->aux;

  if (++info->currentLevel == info->numLevels)
    return;

  assert(info->currentLevel < info->numLevels);

  bfMakeEmptyPtrArrayView(&iter->levelNodes);
  bfPtrArrayGetRangeView(
    &iter->nodes,
    info->offsets[info->currentLevel],
    info->offsets[info->currentLevel + 1],
    &iter->levelNodes);
}

static void next_lrReverseLevelOrder(BfTreeLevelIter *iter) {
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

void bfTreeLevelIterNext(BfTreeLevelIter *iter) {
  if (iter->traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER) {
    next_lrLevelOrder(iter);
    return;
  }

  if (iter->traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    next_lrReverseLevelOrder(iter);
    return;
  }

  bfSetError(BF_ERROR_INVALID_ARGUMENTS);
}

static void deinit_levelOrder(BfTreeLevelIter *iter) {
  LrLevelOrderInfo *info = iter->aux;
  free(info->offsets);
  free(info);

  iter->aux = NULL;
}

void bfTreeLevelIterDeinit(BfTreeLevelIter *iter) {
  bfPtrArrayDeinit(&iter->nodes);

  if (iter->traversal == BF_TREE_TRAVERSAL_LR_LEVEL_ORDER ||
      iter->traversal == BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER) {
    deinit_levelOrder(iter);
    return;
  }

  bfSetError(BF_ERROR_INVALID_ARGUMENTS);
}

BfSize bfTreeLevelIterGetNumPoints(BfTreeLevelIter const *iter) {
  BfPtrArray const *levelNodes = &iter->levelNodes;
  BfSize numPoints = 0;
  for (BfSize i = 0; i < bfPtrArraySize(levelNodes); ++i) {
    BfTreeNode const *node = bfPtrArrayGet(levelNodes, i);
    numPoints += bfTreeNodeGetNumPoints(node);
  }
  return numPoints;
}

bool bfTreeLevelIterCurrentLevelIsInternal(BfTreeLevelIter const *iter) {
  BfPtrArray const *levelNodes = &iter->levelNodes;
  for (BfSize i = 0; i < bfPtrArraySize(levelNodes); ++i) {
    BfTreeNode const *treeNode = bfPtrArrayGet(levelNodes, i);
    if (bfTreeNodeIsLeaf(treeNode))
      return false;
  }
  return true;
}

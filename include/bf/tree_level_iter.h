#pragma once

#include "ptr_array.h"
#include "tree_traversals.h"
#include "types.h"

struct BfTreeLevelIter {
  BfTreeTraversal traversal;
  BfPtrArray nodes;
  BfPtrArray *levelNodes;
  void *aux;
};

BfTreeLevelIter *bfTreeLevelIterNewFromTree(BfTreeTraversal traversal, BfTree *tree);
void bfTreeLevelIterInit(BfTreeLevelIter *iter, BfTreeTraversal traversal, BfTreeNode *node);
void bfTreeLevelIterDeinit(BfTreeLevelIter *iter);
BfSize bfTreeLevelIterCurrentDepth(BfTreeLevelIter const *iter);
bool bfTreeLevelIterIsDone(BfTreeLevelIter const *iter);
void bfTreeLevelIterNext(BfTreeLevelIter *iter);
BfSize bfTreeLevelIterGetNumPoints(BfTreeLevelIter const *iter);
bool bfTreeLevelIterCurrentLevelIsInternal(BfTreeLevelIter const *iter);
BfPtrArray *bfTreeLevelIterGetLevelNodes(BfTreeLevelIter *iter);

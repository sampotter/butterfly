#pragma once

#include "tree_traversals.h"

struct BfTreeLevelIter {
  BfTreeTraversal traversal;
  BfPtrArray nodes;
  BfPtrArray levelNodes;
  void *aux;
};

BfTreeLevelIter bfTreeLevelIterInit(BfTreeTraversal traversal, BfTreeNode *node);
void bfTreeLevelIterFree(BfTreeLevelIter *iter);
BfSize bfTreeLevelIterCurrentDepth(BfTreeLevelIter const *iter);
bool bfTreeLevelIterIsDone(BfTreeLevelIter const *iter);
void bfTreeLevelIterNext(BfTreeLevelIter *iter);
BfSize bfTreeLevelIterGetNumPoints(BfTreeLevelIter const *iter);

#pragma once

#include <bf/ptr_array.h>

typedef enum BfTreeTraversals {
  BF_TREE_TRAVERSAL_UNKNOWN,
  BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
  BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER,
  BF_TREE_TRAVERSAL_POST_ORDER
} BfTreeTraversal;

typedef struct BfTree {
} BfTree;

typedef struct BfTreeNode {
} BfTreeNode;

typedef struct BfTreeIter BfTreeIter;

typedef struct BfTreeLevelIter BfTreeLevelIter;

BfSize bfTreeGetNumPoints(BfTree const *tree);
void bfTreeMap(BfTree *tree, BfTreeTraversal traversal, BfMapFunc func, void *arg);
BfTreeLevelIter *bfTreeGetLevelIter(BfTree *tree, BfTreeTraversal traversal, BfSize skip);

BfSize bfTreeNodeGetNumPoints(BfTreeNode const *node);
BfSize bfTreeNodeNumChildren(BfTreeNode const *node);
bool bfTreeNodeIsLeaf(BfTreeNode const *node);

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

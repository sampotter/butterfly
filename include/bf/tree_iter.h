#pragma once

#include "tree_traversals.h"
#include "types.h"

typedef struct BfTreeIter {
  BfTreeTraversal traversal;
  BfTreeNode *currentNode;
} BfTreeIter;

bool bfTreeIterIsDone(BfTreeIter const *iter);

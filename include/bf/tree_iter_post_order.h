#pragma once

#include "ptr_array.h"
#include "tree_iter.h"

/** Interface: TreeIter */

BfType bfTreeIterPostOrderGetType(BfTreeIterPostOrder const *iter);
BfTreeNode *bfTreeIterPostOrderGetCurrentNode(BfTreeIterPostOrder *iter);
bool bfTreeIterPostOrderIsDone(BfTreeIterPostOrder const *iter);
void bfTreeIterPostOrderNext(BfTreeIterPostOrder *iter);

/** Implementation: TreeIterPostOrder */

struct BfTreeIterPostOrder {
  BfTreeIter super;

  /* Stack used internally by the travresal. */
  BfPtrArray stack;
};

BfTreeIter *bfTreeIterPostOrderToTreeIter(BfTreeIterPostOrder *iter);

BfTreeIterPostOrder *bfTreeIterPostOrderNew();
void bfTreeIterPostOrderInit(BfTreeIterPostOrder *iter, BfTree const *tree);
void bfTreeIterPostOrderDeinit(BfTreeIterPostOrder *iter);

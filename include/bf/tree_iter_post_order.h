#pragma once

#include "ptr_array.h"
#include "tree_iter.h"

/** Interface: TreeIter */

BfType bfTreeIterPostOrderGetType(BfTreeIterPostOrder const *iter);
void bfTreeIterPostOrderDelete(BfTreeIterPostOrder **iter);
BfTreeNode *bfTreeIterPostOrderGetCurrentNode(BfTreeIterPostOrder *iter);
bool bfTreeIterPostOrderIsDone(BfTreeIterPostOrder const *iter);
void bfTreeIterPostOrderNext(BfTreeIterPostOrder *iter);

/** Implementation: TreeIterPostOrder */

struct BfTreeIterPostOrder {
  BfTreeIter super;

  /* Stack used internally by the traversal. */
  BfPtrArray stack;
};

BfTreeIter *bfTreeIterPostOrderToTreeIter(BfTreeIterPostOrder *iter);

BfTreeIterPostOrder *bfTreeIterPostOrderNew(void);
void bfTreeIterPostOrderInit(BfTreeIterPostOrder *iter, BfTree const *tree);
void bfTreeIterPostOrderDeinit(BfTreeIterPostOrder *iter);
void bfTreeIterPostOrderDealloc(BfTreeIterPostOrder **iter);

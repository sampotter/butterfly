#pragma once

#include "tree_traversals.h"
#include "types.h"

/** Interface: TreeIter */

BfType bfTreeIterGetType(BfTreeIter const *iter);
void bfTreeIterDelete(BfTreeIter **treeIter);
BfTreeNode *bfTreeIterGetCurrentNode(BfTreeIter *iter);
bool bfTreeIterIsDone(BfTreeIter const *iter);
void bfTreeIterNext(BfTreeIter *iter);

typedef struct BfTreeIterVtable {
  __typeof__(&bfTreeIterGetType) GetType;
  __typeof__(&bfTreeIterDelete) Delete;
  __typeof__(&bfTreeIterGetCurrentNode) GetCurrentNode;
  __typeof__(&bfTreeIterIsDone) IsDone;
  __typeof__(&bfTreeIterNext) Next;
} BfTreeIterVtable;

/** Implementation: TreeIter */

struct BfTreeIter {
  BfTreeIterVtable *vtable;
  BfTree const *tree;
};

void bfTreeIterInit(BfTreeIter *iter, BfTreeIterVtable *vtable, BfTree const *tree);
void bfTreeIterDeinit(BfTreeIter *treeIter);

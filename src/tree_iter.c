#include <bf/tree_iter.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

/** Interface: TreeIter */

BfType bfTreeIterGetType(BfTreeIter const *iter) {
  return iter->vtable->GetType(iter);
}

BfTreeNode *bfTreeIterGetCurrentNode(BfTreeIter *iter) {
  return iter->vtable->GetCurrentNode(iter);
}

bool bfTreeIterIsDone(BfTreeIter const *iter) {
  return iter->vtable->IsDone(iter);
}

void bfTreeIterNext(BfTreeIter *iter) {
  return iter->vtable->Next(iter);
}

/** Implementation: TreeIter */

void bfTreeIterInit(BfTreeIter *iter, BfTreeIterVtable *vtable, BfTree const *tree) {
  iter->vtable = vtable;
  iter->tree = tree;
}

void bfTreeIterDeinit(BfTreeIter *iter) {
  iter->vtable = NULL;
  iter->tree = NULL;
}

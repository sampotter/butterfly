#include <bf/node_span.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/tree_node.h>

BfNodeSpan *bfNodeSpanNewFromPtrArray(BfPtrArray *ptrArray, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfNodeSpan *nodeSpan = bfMemAlloc(1, sizeof(BfNodeSpan));
  HANDLE_ERROR();

  switch (policy) {
  case BF_POLICY_VIEW:
    nodeSpan->nodes = bfPtrArrayGetView(ptrArray);
    break;
  case BF_POLICY_COPY:
    nodeSpan->nodes = bfPtrArrayCopy(ptrArray);
    break;
  case BF_POLICY_STEAL:
    nodeSpan->nodes = bfPtrArraySteal(ptrArray);
    break;
  }
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return nodeSpan;
}

BfSize bfNodeSpanGetSize(BfNodeSpan const *nodeSpan) {
  return bfPtrArraySize(nodeSpan->nodes);
}

bool bfNodeSpanIsEmpty(BfNodeSpan const *nodeSpan) {
  return bfPtrArrayIsEmpty(nodeSpan->nodes);
}

BfType bfNodeSpanGetNodeType(BfNodeSpan const *nodeSpan) {
  return bfNodeSpanIsEmpty(nodeSpan) ?
    BF_TYPE_NONE :
    bfTreeNodeGetType(bfNodeSpanGetFirstNode(nodeSpan));
}

BfTree *bfNodeSpanGetTree(BfNodeSpan const *nodeSpan) {
  return bfNodeSpanIsEmpty(nodeSpan) ?
    NULL :
    bfTreeNodeGetTree(bfNodeSpanGetNode(nodeSpan, 0));
}

bool bfNodeSpanContains(BfNodeSpan const *nodeSpan, BfTreeNode const *node) {
  /* TODO: we could speed this up by assuming that the nodes in
   * `nodeSpan` are sorted and then sorting using `node`'s i0 and i1
   * indices. */
  for (BfSize i = 0; i < bfNodeSpanGetSize(nodeSpan); ++i)
    if (bfNodeSpanGetNode(nodeSpan, i) == node)
      return true;
  return false;
}

bool bfNodeSpanContainsNodeWithSameRange(BfNodeSpan const *nodeSpan, BfTreeNode const *node) {
  /* TODO: we could speed this up by assuming that the nodes in
   * `nodeSpan` are sorted and then sorting using `node`'s i0 and i1
   * indices. */
  BfSize i0 = bfTreeNodeGetFirstIndex(node);
  BfSize i1 = bfTreeNodeGetLastIndex(node);
  for (BfSize i = 0; i < bfNodeSpanGetSize(nodeSpan); ++i) {
    BfTreeNode const *otherNode = bfNodeSpanGetNode(nodeSpan, i);
    if (i0 == bfTreeNodeGetFirstIndex(otherNode) &&
        i1 == bfTreeNodeGetLastIndex(otherNode))
      return true;
  }
  return false;
}

BfTreeNode *bfNodeSpanGetNode(BfNodeSpan const *nodeSpan, BfSize i) {
  BF_ERROR_BEGIN();

  BfTreeNode *node = NULL;

  if (i >= bfNodeSpanGetSize(nodeSpan))
    RAISE_ERROR(BF_ERROR_OUT_OF_RANGE);

  node = bfPtrArrayGet(nodeSpan->nodes, i);

  BF_ERROR_END() {
    BF_DIE();
  }

  return node;
}

BfTreeNode *bfNodeSpanGetRoot(BfNodeSpan const *nodeSpan) {
  BfTreeNode *root = NULL;
  if (!bfNodeSpanIsEmpty(nodeSpan)) {
    BfSize i0 = bfNodeSpanGetFirstIndex(nodeSpan);
    BfSize i1 = bfNodeSpanGetLastIndex(nodeSpan);
    root = bfNodeSpanGetFirstNode(nodeSpan);
    while (root->parent &&
           i0 <= bfTreeNodeGetFirstIndex(root->parent) &&
           bfTreeNodeGetLastIndex(root->parent) <= i1)
      root = root->parent;
    BF_ASSERT(bfTreeNodeGetFirstIndex(root) == bfNodeSpanGetFirstIndex(nodeSpan));
    BF_ASSERT(bfTreeNodeGetLastIndex(root) == bfNodeSpanGetLastIndex(nodeSpan));
  }
  return root;
}

BfSize bfNodeSpanGetFirstIndex(BfNodeSpan const *nodeSpan) {
  return bfNodeSpanIsEmpty(nodeSpan) ?
    BF_SIZE_BAD_VALUE :
    bfTreeNodeGetFirstIndex(bfNodeSpanGetFirstNode(nodeSpan));
}

BfSize bfNodeSpanGetLastIndex(BfNodeSpan const *nodeSpan) {
  return bfNodeSpanIsEmpty(nodeSpan) ?
    BF_SIZE_BAD_VALUE :
    bfTreeNodeGetLastIndex(bfNodeSpanGetLastNode(nodeSpan));
}

BfTreeNode *bfNodeSpanGetFirstNode(BfNodeSpan const *nodeSpan) {
  return bfNodeSpanIsEmpty(nodeSpan) ?
    NULL :
    bfPtrArrayGet(nodeSpan->nodes, 0);
}

BfTreeNode *bfNodeSpanGetLastNode(BfNodeSpan const *nodeSpan) {
  return bfNodeSpanIsEmpty(nodeSpan) ?
    NULL :
    bfPtrArrayGet(nodeSpan->nodes, bfPtrArraySize(nodeSpan->nodes) - 1);
}

#include <bf/tree_node_span.h>

BfTreeNode const *bfTreeNodeSpanGetByFirstIndex(BfTreeNodeSpan const *span, BfSize i0) {
  BfTreeNode const *node = NULL;
  for (BfSize k = 0; k < bfConstPtrArraySize(nodes); ++k) {
    node = bfConstPtrArrayGet(nodes, k);
    if (bfTreeNodeGetFirstIndex(node) == i0)
      break;
  }
  return node;
}

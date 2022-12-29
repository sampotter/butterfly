#include <bf/interval_tree_node.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

/** Interface: */

static BfTreeNodeVtable TreeNodeVtable = {
  .GetType = (__typeof__(&bfTreeNodeGetType))bfIntervalTreeNodeGetType
};

BfType bfIntervalTreeNodeGetType(BfTreeNode const *treeNode) {
  (void)treeNode;
  return BF_TYPE_QUADTREE_NODE;
}

/** Upcasting: IntervalTreeNode -> TreeNode */

BfTreeNode *bfIntervalTreeNodeToTreeNode(BfIntervalTreeNode *node) {
  return &node->super;
}

BfTreeNode const *bfIntervalTreeNodeConstToTreeNodeConst(BfIntervalTreeNode const *node) {
  return &node->super;
}

/** Downcasting: TreeNode -> IntervalTreeNode */

BfIntervalTreeNode *bfTreeNodeToIntervalTreeNode(BfTreeNode *treeNode) {
  if (!bfTreeNodeInstanceOf(treeNode, BF_TYPE_INTERVAL_TREE_NODE)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfIntervalTreeNode *)treeNode;
  }
}

BfIntervalTreeNode const *bfTreeNodeConstToIntervalTreeNodeConst(BfTreeNode const *treeNode) {
  if (!bfTreeNodeInstanceOf(treeNode, BF_TYPE_INTERVAL_TREE_NODE)) {
    bfSetError(BF_ERROR_RUNTIME_ERROR);
    return NULL;
  } else {
    return (BfIntervalTreeNode const *)treeNode;
  }
}

BfIntervalTreeNode *bfIntervalTreeNodeNew() {
  BEGIN_ERROR_HANDLING();

  BfIntervalTreeNode *node = malloc(sizeof(BfIntervalTreeNode));
  if (node == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    node = NULL;

  return node;
}

void bfIntervalTreeNodeInitEmpty(BfIntervalTreeNode *intervalTreeNode,
                                 BfIntervalTree const *intervalTree,
                                 BfReal a, BfReal b, BfSize k) {
  BEGIN_ERROR_HANDLING();

  (void)a;
  (void)b;

  assert(false);

  bfTreeNodeInit(&intervalTreeNode->super, &TreeNodeVtable,
                 true, (void *)intervalTree, k, BF_SIZE_BAD_VALUE, 0);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

#include <bf/interval_tree_node.h>

#include <assert.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

/** Interface: */

static BfTreeNodeVtable TreeNodeVtable = {
  .GetType = (__typeof__(&bfTreeNodeGetType))bfIntervalTreeNodeGetType
};

BfType bfIntervalTreeNodeGetType(BfTreeNode const *treeNode) {
  (void)treeNode;
  return BF_TYPE_INTERVAL_TREE_NODE;
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
  BF_ERROR_BEGIN();

  BfIntervalTreeNode *node = bfMemAlloc(1, sizeof(BfIntervalTreeNode));
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    node = NULL;

  return node;
}

static void intervalTreeNodeInitEmptyRecursive(BfIntervalTreeNode *intervalTreeNode,
                                               BfSize k, BfSize depth, BfSize maxDepth) {
  BF_ERROR_BEGIN();

  BfReal const a = intervalTreeNode->a;
  BfReal const b = intervalTreeNode->b;
  BfReal const delta = (b - a)/k;

  BfTreeNode *treeNode = bfIntervalTreeNodeToTreeNode(intervalTreeNode);

  if (depth < maxDepth) {
    for (BfSize i = 0; i < k; ++i) {
      BfIntervalTreeNode *child = bfIntervalTreeNodeNew();
      HANDLE_ERROR();

      bfTreeNodeInit(&child->super, &TreeNodeVtable, false, (void *)treeNode, k, i, depth + 1);
      HANDLE_ERROR();

      child->a = a + delta*i;

      /* Avoid any possibility of floating-point error here */
      child->b = i == k - 1 ? b : a + delta*(i + 1);

      child->isLeftmost = intervalTreeNode->isLeftmost && i == 0;
      child->isRightmost = intervalTreeNode->isRightmost && i == k - 1;

      intervalTreeNodeInitEmptyRecursive(child, k, depth + 1, maxDepth);
      HANDLE_ERROR();

      treeNode->child[i] = bfIntervalTreeNodeToTreeNode(child);
    }
  }

  END_ERROR_HANDLING() {}
}

void bfIntervalTreeNodeInitEmptyRoot(BfIntervalTreeNode *intervalTreeNode,
                                     BfIntervalTree const *intervalTree,
                                     BfReal a, BfReal b, BfSize k, BfSize maxDepth) {
  BF_ERROR_BEGIN();

  bfTreeNodeInit(&intervalTreeNode->super, &TreeNodeVtable,
                 true, (void *)intervalTree, k, BF_SIZE_BAD_VALUE, 0);
  HANDLE_ERROR();

  intervalTreeNode->a = a;
  intervalTreeNode->b = b;
  intervalTreeNode->isLeftmost = true;
  intervalTreeNode->isRightmost = true;

  intervalTreeNodeInitEmptyRecursive(intervalTreeNode, k, 0, maxDepth);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    // bfIntervalTreeNodeDeinit(intervalTreeNode);
  }
}

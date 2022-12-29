#include "interval_tree.h"

#include "tree_node.h"

/** Interface: */

BfType bfIntervalTreeNodeGetType(BfTreeNode const *treeNode);

/** Upcasting: IntervalTreeNode -> TreeNode */

BfTreeNode *bfIntervalTreeNodeToTreeNode(BfIntervalTreeNode *node);
BfTreeNode const *bfIntervalTreeNodeConstToTreeNodeConst(BfIntervalTreeNode const *node);

/** Downcasting: TreeNode -> IntervalTreeNode */

BfIntervalTreeNode *bfTreeNodeToIntervalTreeNode(BfTreeNode *treeNode);
BfIntervalTreeNode const *bfTreeNodeConstToIntervalTreeNodeConst(BfTreeNode const *treeNode);

/** Implementation: IntervalTreeNode */

struct BfIntervalTreeNode {
  BfTreeNode super;

  /* The interval for the current node. If `isRightmost` is true, this
   * is interpreted as the interval `[a, b]`, otherwise `[a, b)`. */
  BfReal a, b;

  /* Whether this node corresponds to the rightmost interval at the
   * current level. */
  bool isRightmost;
};

BfIntervalTreeNode *bfIntervalTreeNodeNew();
void bfIntervalTreeNodeInitEmpty(BfIntervalTreeNode *node,
                                 BfIntervalTree const *intervalTree,
                                 BfReal a, BfReal b, BfSize k);

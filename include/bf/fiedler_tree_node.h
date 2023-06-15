#include "fiedler_tree.h"
#include "tree_node.h"

/** Interface: TreeNode */

BfType bfFiedlerTreeNodeGetType(BfFiedlerTreeNode const *node);
void bfFiedlerTreeNodeDelete(BfFiedlerTreeNode **node);

/** FiedlerTreeNode: */

/** Upcasting: FiedlerTreeNode -> TreeNode */

BfTreeNode *bfFiedlerTreeNodeToTreeNode(BfFiedlerTreeNode *node);

/** Downcasting: TreeNode -> FiedlerTreeNode */

/** Implementation: FiedlerTreeNode */

struct BfFiedlerTreeNode {
  BfTreeNode super;
};

BfFiedlerTreeNode *bfFiedlerTreeNodeNew();
void bfFiedlerTreeNodeInitRoot(BfFiedlerTreeNode *node, BfFiedlerTree const *tree);
void bfFiedlerTreeNodeDeinit(BfFiedlerTreeNode *node);
void bfFiedlerTreeNodeDealloc(BfFiedlerTreeNode **node);
void bfFiedlerTreeNodeDeinitAndDealloc(BfFiedlerTreeNode **node);

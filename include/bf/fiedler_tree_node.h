#include "fiedler_tree.h"
#include "tree_node.h"

/** Interface: TreeNode */

BfType bfFiedlerTreeNodeGetType(BfFiedlerTreeNode const *node);
void bfFiedlerTreeNodeDelete(BfFiedlerTreeNode **node);

/** FiedlerTreeNode: */

/** Upcasting: FiedlerTreeNode -> TreeNode */

BfTreeNode *bfFiedlerTreeNodeToTreeNode(BfFiedlerTreeNode *node);

/** Downcasting: TreeNode -> FiedlerTreeNode */

BfFiedlerTreeNode const *bfTreeNodeConstToFiedlerTreeNodeConst(BfTreeNode const *treeNode);

/** Implementation: FiedlerTreeNode */

struct BfFiedlerTreeNode {
  BfTreeNode super;

  /* The triangle mesh corresponding to this node. This is stored
   * optionally (see `keepNodeTrimeshes` flag in the creation
   * functions). */
  BfTrimesh const *trimesh;
};

BfFiedlerTreeNode *bfFiedlerTreeNodeNew(void);
void bfFiedlerTreeNodeInitRoot(BfFiedlerTreeNode *node, BfFiedlerTree const *tree, BfReal tol, bool keepNodeTrimeshes);
void bfFiedlerTreeNodeDeinit(BfFiedlerTreeNode *node);
void bfFiedlerTreeNodeDealloc(BfFiedlerTreeNode **node);
void bfFiedlerTreeNodeDeinitAndDealloc(BfFiedlerTreeNode **node);

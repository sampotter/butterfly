#pragma once

#include "def.h"
#include "types.h"

typedef void (*BfTreeMapFunc)(BfTree *, BfTreeNode *, void *);
typedef void (*BfTreeMapConstFunc)(BfTree const *, BfTreeNode const *, void *);

/** Interface: TreeNode */

BfTreeNode *bfTreeNodeCopy(BfTreeNode const *treeNode);
BfType bfTreeNodeGetType(BfTreeNode const *treeNode);
void bfTreeNodeDelete(BfTreeNode **treeNode);

typedef struct BfTreeNodeVtable {
  __typeof__(&bfTreeNodeCopy) Copy;
  __typeof__(&bfTreeNodeGetType) GetType;
  __typeof__(&bfTreeNodeDelete) Delete;
} BfTreeNodeVtable;

/** Implementation: TreeNode */

struct BfTreeNode {
  BfTreeNodeVtable *vtbl;

  /* Whether or not this node is the root node. */
  // TODO: we can combine `isRoot` and `index` by signaling `isRoot`
  // with `index == BF_SIZE_BAD_VALUE`
  bool isRoot;

  /* The index of this node (`super.child[index]` gives this node). */
  BfSize index;

  /* The parent of this node. If this is the root node, then this is a
   * pointer to the containing `BfTree`. Otherwise, it points to the
   * parent node, which will be a `BfTreeNode`. */
  void *parent;

  BfSize maxNumChildren;

  /* An array of length `maxNumChildren + 1` containing index offsets for
   * this node. These are relative to this node's parent. To find the
   * absolute indices, we follow the parent pointer back to the
   * root. The first and last entries hold sentinel values: offset[0]
   * == 0 and offset[maxNumChildren] == size (the number of points contained in
   * this node). */
  BfSize *offset;

  /* Array of length `maxNumChildren` containing pointers to this
   * node's children. If child[i] == NULL, then there is no child node
   * and offset[i + 1] - offset[i] == 0 holds. */
  BfTreeNode **child;

  /* The depth of this node. */
  BfSize depth;
};

BfTreeNode *bfTreeNodeAlloc(void);
void bfTreeNodeInitRoot(BfTreeNode *treeNode, BfTreeNodeVtable *vtbl, BfTree const *tree, BfSize maxNumChildren);
void bfTreeNodeInit(BfTreeNode *treeNode, BfTreeNodeVtable *vtbl, bool isRoot, void *parent, BfSize maxNumChildren, BfSize index, BfSize depth);
void bfTreeNodeDeinit(BfTreeNode *treeNode);
BfTreeNodeVtable *bfTreeNodeGetVtable(void);
bool bfTreeNodeInstanceOf(BfTreeNode const *treeNode, BfType type);
BfSize bfTreeNodeGetMaxDepthBelow(BfTreeNode const *treeNode);
BfTreeNode *bfTreeNodeGetParent(BfTreeNode *treeNode);
BfTreeNode const *bfTreeNodeGetParentConst(BfTreeNode const *treeNode);
BfSize bfTreeNodeGetNumChildren(BfTreeNode const *treeNode);
BfSize bfTreeNodeGetDepth(BfTreeNode const *treeNode);
BfSize bfTreeNodeGetNumPoints(BfTreeNode const *treeNode);
bool bfTreeNodeIsRoot(BfTreeNode const *treeNode);
bool bfTreeNodeIsLeaf(BfTreeNode const *treeNode);
BfSize bfTreeNodeGetFirstIndex(BfTreeNode const *treeNode);
BfSize bfTreeNodeGetLastIndex(BfTreeNode const *treeNode);
BfSize const *bfTreeNodeGetIndexPtrConst(BfTreeNode const *treeNode, BfTree const *tree);
BfTree *bfTreeNodeGetTree(BfTreeNode *treeNode);
BfTree const *bfTreeNodeGetTreeConst(BfTreeNode const *treeNode);
bool bfTreeNodeIsEmpty(BfTreeNode const *treeNode);
bool bfTreeNodeIsDescendant(BfTreeNode const *treeNode, BfTreeNode const *otherTreeNode);
BfTreeNode *bfTreeNodeGetChild(BfTreeNode *treeNode, BfSize i);
BfSize bfTreeNodeGetMaxNumChildren(BfTreeNode const *treeNode);
bool bfTreeNodeHasChild(BfTreeNode const *treeNode, BfSize i);

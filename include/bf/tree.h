#pragma once

#include "node_span.h"
#include "perm.h"
#include "ptr_array.h"
#include "tree_traversals.h"
#include "types.h"

enum BfTreeNodeFlags {
  BF_TREE_NODE_FLAG_NONE = 0
};

/** Interface: Tree */

typedef void (*BfTreeMapFunc)(BfTree *, BfTreeNode *, void *);
typedef void (*BfTreeMapConstFunc)(BfTree const *, BfTreeNode const *, void *);

void bfTreeCopyInto(BfTree *tree, BfTree *dstTree);
void bfTreeDelete(BfTree **tree);
BfType bfTreeGetType(BfTree const *tree);

typedef struct BfTreeVtable {
  __typeof__(&bfTreeCopyInto) CopyInto;
  __typeof__(&bfTreeDelete) Delete;
  __typeof__(&bfTreeGetType) GetType;
} BfTreeVtable;

/** Implementation: Tree */

struct BfTree {
  BfTreeVtable *vtable;

  /* Root node of the tree. */
  BfTreeNode *root;

  /* Permutation from domain ordering to the tree ordering. To undo,
   * use `bfPermGetReversePerm` to get the inverse permutation. */
  BfPerm *perm;
};

BfTree *bfTreeAlloc(void);
BfTree *bfTreeNewFromNodeType(BfType type);
BfTree *bfTreeNewFromNodeSpan(BfNodeSpan *nodeSpan, BfPerm **perm);
BfTree *bfTreeNewFromNode(BfTreeNode *node, BfPerm **perm);
BfTree *bfTreeNewForMiddleFac(BfTree const *templateTree, BfSize p);
void bfTreeInit(BfTree *tree, BfTreeVtable *vtbl, BfTreeNode *root, BfSize size);
void bfTreeInitFromNodeSpan(BfTree *tree, BfNodeSpan *nodeSpan, BfPerm **perm);
void bfTreeInitFromNode(BfTree *tree, BfTreeNode *node, BfPerm **perm);
void bfTreeInitForMiddleFac(BfTree *tree, BfTree const *templateTree, BfSize p);
void bfTreeDeinit(BfTree *tree);
bool bfTreeInstanceOf(BfTree const *tree, BfType type);
BfTreeNode *bfTreeGetRootNode(BfTree *);
BfTreeNode const *bfTreeGetRootNodeConst(BfTree const *);
BfPerm *bfTreeGetPerm(BfTree *tree);
BfPerm const *bfTreeGetPermConst(BfTree const *tree);
BfSize bfTreeGetMaxDepth(BfTree const *);
void bfTreeMap(BfTree *, BfTreeNode *, BfTreeTraversal, BfTreeMapFunc, void *);
void bfTreeMapConst(BfTree const *, BfTreeNode const *, BfTreeTraversal, BfTreeMapConstFunc, void *);
BfSize bfTreeGetNumPoints(BfTree const *);
BfTreeNode *bfTreeGetNode(BfTree *, BfSize, BfSize);
BfPtrArray *bfTreeGetLevelPtrArray(BfTree *, BfSize);

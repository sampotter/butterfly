#pragma once

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

BfType bfTreeGetType(BfTree const *);

typedef struct BfTreeVtable {
  __typeof__(&bfTreeGetType) GetType;
} BfTreeVtable;

/** Implementation: Tree */

struct BfTree {
  BfTreeVtable *vtable;

  /* Root node of the tree. */
  BfTreeNode *root;

  /* Permutation from domain ordering to the tree ordering. To undo,
   * use `bfPermGetReversePerm` to gete the inverse permutation. */
  BfPerm perm;
};

void bfTreeInit(BfTree *tree, BfTreeVtable *vtbl, BfTreeNode *root, BfSize size);
void bfTreeDeinit(BfTree *tree);
bool bfTreeInstanceOf(BfTree const *tree, BfType type);
BfTreeNode *bfTreeGetRootNode(BfTree *);
BfTreeNode const *bfTreeGetRootNodeConst(BfTree const *);
BfPerm const *bfTreeGetPermConst(BfTree const *tree);
BfSize bfTreeGetMaxDepth(BfTree const *);
void bfTreeMap(BfTree *, BfTreeNode *, BfTreeTraversal, BfTreeMapFunc, void *);
void bfTreeMapConst(BfTree const *, BfTreeNode const *, BfTreeTraversal, BfTreeMapConstFunc, void *);
BfSize bfTreeGetNumPoints(BfTree const *);
BfTreeNode *bfTreeGetNode(BfTree *, BfSize, BfSize);
BfPtrArray bfTreeGetLevelPtrArray(BfTree *, BfSize);
BfConstPtrArray bfTreeGetLevelConstPtrArray(BfTree const *, BfSize);

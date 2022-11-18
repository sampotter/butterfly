#pragma once

#include "interface.h"
#include "ptr_array.h"
#include "tree_traversals.h"
#include "types.h"

typedef void (*BfTreeMapFunc)(BfTree *, BfTreeNode *, void *);

#define BF_INTERFACE_Tree(Type, Subtype, _)                             \
  _(Type, Subtype, BfSize, GetMaxDepth, BfTree const *)                 \
  _(Type, Subtype, void, Map, BfTree *, BfTreeTraversal, BfTreeMapFunc, void *)

#define INTERFACE BF_INTERFACE_Tree
BF_DEFINE_VTABLE_STRUCT(Tree)
BF_DECLARE_INTERFACE(Tree)
#undef INTERFACE

struct BfTree {
  BfTreeVtable *vtbl;
};

#define BF_INTERFACE_TreeNode(Type, Subtype, _) \
  _(Type, Subtype, BfSize, GetNumChildren, BfTreeNode const *)  \
  _(Type, Subtype, bool, IsLeaf, BfTreeNode const *)

#define INTERFACE BF_INTERFACE_TreeNode
BF_DEFINE_VTABLE_STRUCT(TreeNode)
BF_DECLARE_INTERFACE(TreeNode)
#undef INTERFACE

struct BfTreeNode {
  BfTreeNodeVtable *vtbl;
};

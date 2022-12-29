#pragma once

#include "mat.h"
#include "tree.h"
#include "trimesh.h"

struct BfFiedlerTreeNode {
  BfTreeNode super;
  BfFiedlerTreeNode *parent;
};

// #define INTERFACE BF_INTERFACE_TreeNode
// BF_DECLARE_INTERFACE(FiedlerTreeNode)
// #undef INTERFACE

struct BfFiedlerTree {
  BfTree base;
  BfFiedlerTreeNode *root;
};

// #define INTERFACE BF_INTERFACE_Tree
// BF_DECLARE_INTERFACE(FiedlerTree)
// #undef INTERFACE

/* Upcasting: */
BfTree *bfFiedlerTreeToTree(BfFiedlerTree *fiedlerTree);

void bfFiedlerTreeInitFromTrimesh(BfFiedlerTree *fiedlerTree, BfTrimesh const *trimesh);

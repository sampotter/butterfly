#pragma once

#include "mat.h"
#include "tree.h"
#include "trimesh.h"

struct BfFiedlerTreeNode {
  BfTreeNode super;
  BfFiedlerTreeNode *parent;
};

struct BfFiedlerTree {
  BfTree base;
  BfFiedlerTreeNode *root;
};

/* Upcasting: */
BfTree *bfFiedlerTreeToTree(BfFiedlerTree *fiedlerTree);

void bfFiedlerTreeInitFromTrimesh(BfFiedlerTree *fiedlerTree, BfTrimesh const *trimesh);

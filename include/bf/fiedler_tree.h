#pragma once

#include "mat.h"
#include "tree.h"
#include "tree_node.h"
#include "trimesh.h"

/** FiedlerTree: */

struct BfFiedlerTree {
  BfTree super;
  BfFiedlerTreeNode *root;
  BfTrimesh const *trimesh;
};

/** Interface(Tree, FiedlerTree) */

BfType bfFiedlerTreeGetType(BfFiedlerTree const *tree);

/** Upcasting: FiedlerTree -> Tree */

BfTree *bfFiedlerTreeToTree(BfFiedlerTree *fiedlerTree);
BfTree const *bfFiedlerTreeConstToTreeConst(BfFiedlerTree const *fiedlerTree);

/** Downcasting: Tree -> FiedlerTree */

/** Implementation: FiedlerTree */

BfFiedlerTree *bfFiedlerTreeNewFromTrimesh(BfTrimesh const *trimesh);
void bfFiedlerTreeInitFromTrimesh(BfFiedlerTree *fiedlerTree, BfTrimesh const *trimesh);
void bfFiedlerTreeDeinit(BfFiedlerTree *fiedlerTree);
void bfFiedlerTreeDealloc(BfFiedlerTree **fiedlerTree);
void bfFiedlerTreeDeinitAndDealloc(BfFiedlerTree **fiedlerTree);

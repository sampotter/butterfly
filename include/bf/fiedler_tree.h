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

BfFiedlerTree const *bfTreeConstToFiedlerTreeConst(BfTree const *tree);

/** Implementation: FiedlerTree */

BfFiedlerTree *bfFiedlerTreeNewFromTrimesh(BfTrimesh const *trimesh, BfReal tol, bool keepNodeTrimeshes);
void bfFiedlerTreeInitFromTrimesh(BfFiedlerTree *fiedlerTree, BfTrimesh const *trimesh, BfReal tol, bool keepNodeTrimeshes);
void bfFiedlerTreeDeinit(BfFiedlerTree *fiedlerTree);
void bfFiedlerTreeDealloc(BfFiedlerTree **fiedlerTree);
void bfFiedlerTreeDeinitAndDealloc(BfFiedlerTree **fiedlerTree);

#include <bf/fiedler_tree.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/fiedler_tree_node.h>
#include <bf/lbo.h>
#include <bf/linalg.h>
#include <bf/mem.h>
#include <bf/real_array.h>
#include <bf/size_array.h>
#include <bf/vec_real.h>

#include "macros.h"

/** Interface(Tree, FiedlerTree) */

static BfTreeVtable TreeVtable = {
  .GetType = (__typeof__(&bfTreeGetType))bfFiedlerTreeGetType,
};

BfType bfFiedlerTreeGetType(BfFiedlerTree const *tree) {
  (void)tree;
  return BF_TYPE_FIEDLER_TREE;
}

BfFiedlerTree *bfFiedlerTreeNewFromTrimesh(BfTrimesh const *trimesh, BfReal tol, bool keepNodeTrimeshes) {
  BF_ERROR_BEGIN();

  BfFiedlerTree *fiedlerTree = bfMemAlloc(1, sizeof(BfFiedlerTree));
  HANDLE_ERROR();

  bfFiedlerTreeInitFromTrimesh(fiedlerTree, trimesh, tol, keepNodeTrimeshes);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return fiedlerTree;
}

/** Upcasting: FiedlerTree -> Tree */

BfTree *bfFiedlerTreeToTree(BfFiedlerTree *tree) {
  return &tree->super;
}

BfTree const *bfFiedlerTreeConstToTreeConst(BfFiedlerTree const *tree) {
  return &tree->super;
}

/** Downcasting: Tree -> FiedlerTree */

BfFiedlerTree const *bfTreeConstToFiedlerTreeConst(BfTree const *tree) {
  if (!bfTreeInstanceOf(tree, BF_TYPE_FIEDLER_TREE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfFiedlerTree const *)tree;
  }
}

/** Implementation: FiedlerTree */

void bfFiedlerTreeInitFromTrimesh(BfFiedlerTree *tree, BfTrimesh const *trimesh, BfReal tol, bool keepNodeTrimeshes) {
  BF_ERROR_BEGIN();

  BfFiedlerTreeNode *root = bfFiedlerTreeNodeNew();
  HANDLE_ERROR();

  bfTreeInit(
    &tree->super,
    &TreeVtable,
    bfFiedlerTreeNodeToTreeNode(root),
    bfTrimeshGetNumVerts(trimesh));
  HANDLE_ERROR();

  tree->trimesh = trimesh;

  bfFiedlerTreeNodeInitRoot(root, tree, tol, keepNodeTrimeshes);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfFiedlerTreeDeinit(BfFiedlerTree *fiedlerTree) {
  bfTreeDeinit(&fiedlerTree->super);
}

void bfFiedlerTreeDealloc(BfFiedlerTree **fiedlerTree) {
  bfMemFree(*fiedlerTree);
  *fiedlerTree = NULL;
}

void bfFiedlerTreeDeinitAndDealloc(BfFiedlerTree **fiedlerTree) {
  bfFiedlerTreeDeinit(*fiedlerTree);
  bfFiedlerTreeDealloc(fiedlerTree);
}

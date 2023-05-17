#include <bf/tree.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/tree_level_iter.h>
#include <bf/tree_node.h>

/** Interface: Tree */

void bfTreeDelete(BfTree **tree) {
  (*tree)->vtable->Delete(tree);
}

BfType bfTreeGetType(BfTree const *tree) {
  return tree->vtable->GetType(tree);
}

/** Implementation: Tree */

void bfTreeInit(BfTree *tree, BfTreeVtable *vtable, BfTreeNode *root, BfSize size) {
  BF_ERROR_BEGIN();

  tree->vtable = vtable;
  tree->root = root;

  tree->perm = bfPermIdentity(size);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfTreeDeinit(tree);
}

void bfTreeDeinit(BfTree *tree) {
  (void)tree;
  BF_DIE();
}

bool bfTreeInstanceOf(BfTree const *tree, BfType type) {
  return bfTypeDerivedFrom(bfTreeGetType(tree), type);
}

BfTreeNode *bfTreeGetRootNode(BfTree *tree) {
  return tree->root;
}

BfTreeNode const *bfTreeGetRootNodeConst(BfTree const *tree) {
  return tree->root;
}

BfPerm const *bfTreeGetPermConst(BfTree const *tree) {
  return &tree->perm;
}

BfSize bfTreeGetMaxDepth(BfTree const *tree) {
  return bfTreeNodeGetMaxDepthBelow(tree->root);
}

static void
map_levelOrder(BfTree *tree, BfTreeNode *node,
               BfTreeTraversal traversal, BfTreeMapFunc func,
               void *arg) {
  BF_ERROR_BEGIN();

  BfTreeLevelIter iter;
  bfTreeLevelIterInit(&iter, traversal, node);
  HANDLE_ERROR();

  while (!bfTreeLevelIterIsDone(&iter)) {
    for (BfSize i = 0; i < bfPtrArraySize(&iter.levelNodes); ++i) {
      BfTreeNode *currentNode = bfPtrArrayGet(&iter.levelNodes, i);
      func(tree, currentNode, arg);
      HANDLE_ERROR();
    }
    bfTreeLevelIterNext(&iter);
  }

  bfTreeLevelIterDeinit(&iter);

  BF_ERROR_END() {}
}

void bfTreeMap(BfTree *tree, BfTreeNode *node,
                   BfTreeTraversal traversal, BfTreeMapFunc func,
                   void *arg) {
  BF_ERROR_BEGIN();

  if (node == NULL)
    node = tree->root;

  switch (traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER: {
    map_levelOrder(tree, node, traversal, func, arg);
    HANDLE_ERROR();
    break;
  } default:
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  BF_ERROR_END() {}
}

static void
mapConst_levelOrder(BfTree const *tree, BfTreeNode const *node,
                    BfTreeTraversal traversal, BfTreeMapConstFunc func,
                    void *arg) {
  BF_ERROR_BEGIN();

  BfTreeLevelIter iter;
  bfTreeLevelIterInit(&iter, traversal, (BfTreeNode *)node);
  HANDLE_ERROR();

  while (!bfTreeLevelIterIsDone(&iter)) {
    for (BfSize i = 0; i < bfPtrArraySize(&iter.levelNodes); ++i) {
      BfTreeNode const *currentNode = bfPtrArrayGet(&iter.levelNodes, i);
      func(tree, currentNode, arg);
      HANDLE_ERROR();
    }
    bfTreeLevelIterNext(&iter);
  }

  bfTreeLevelIterDeinit(&iter);

  BF_ERROR_END() {}
}

void bfTreeMapConst(BfTree const *tree, BfTreeNode const *node,
                    BfTreeTraversal traversal, BfTreeMapConstFunc func, void *arg) {
  BF_ERROR_BEGIN();

  if (node == NULL)
    node = tree->root;

  switch (traversal) {
  case BF_TREE_TRAVERSAL_LR_LEVEL_ORDER:
  case BF_TREE_TRAVERSAL_LR_REVERSE_LEVEL_ORDER: {
    mapConst_levelOrder(tree, node, traversal, func, arg);
    HANDLE_ERROR();
    break;
  } default:
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
  }

  BF_ERROR_END() {}
}

BfSize bfTreeGetNumPoints(BfTree const *tree) {
  return bfTreeNodeGetNumPoints(tree->root);
}

BfTreeNode *bfTreeGetNode(BfTree *tree, BfSize depth, BfSize nodeIndex) {
  BF_ERROR_BEGIN();

  BfPtrArray levelNodes;
  BfTreeNode *node = NULL;

  levelNodes = bfTreeGetLevelPtrArray(tree, depth);
  HANDLE_ERROR();

  node = bfPtrArrayGet(&levelNodes, nodeIndex);
  HANDLE_ERROR();

  BF_ERROR_END() {
    node = NULL;
  }

  bfPtrArrayDeinit(&levelNodes);

  return node;
}

BfPtrArray bfTreeGetLevelPtrArray(BfTree *tree, BfSize depth) {
  BF_ERROR_BEGIN();

  BfPtrArray levelNodes = bfGetUninitializedPtrArray();

  BfSize maxDepth = bfTreeGetMaxDepth(tree);

  if (depth > maxDepth)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfTreeLevelIter iter;
  bfTreeLevelIterInit(&iter, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, tree->root);
  HANDLE_ERROR();

  for (BfSize _ = 0; _ < depth; ++_)
    bfTreeLevelIterNext(&iter);

  levelNodes = bfPtrArrayCopy(&iter.levelNodes);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfPtrArrayDeinit(&levelNodes);

  bfTreeLevelIterDeinit(&iter);

#if BF_DEBUG
  for (BfSize k = 0; k < bfPtrArraySize(&levelNodes); ++k)
    BF_ASSERT(bfPtrArrayGet(&levelNodes, k) != NULL);
#endif

  return levelNodes;
}

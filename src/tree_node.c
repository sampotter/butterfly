#include <bf/tree_node.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/tree.h>

/** Interface: TreeNode */

BfTreeNode *bfTreeNodeCopy(BfTreeNode const *node) {
  return node->vtbl->Copy(node);
}

BfType bfTreeNodeGetType(BfTreeNode const *node) {
  return node->vtbl->GetType(node);
}

void bfTreeNodeDelete(BfTreeNode **treeNode) {
  (*treeNode)->vtbl->Delete(treeNode);
}

/** Implementation: TreeNode */

BfTreeNode *bfTreeNodeAlloc(void) {
  BF_ERROR_BEGIN();

  BfTreeNode *node = bfMemAlloc(1, sizeof(BfTreeNode));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return node;
}

void bfTreeNodeInitRoot(BfTreeNode *treeNode, BfTreeNodeVtable *vtbl, BfTree const *tree, BfSize maxNumChildren) {
  bfTreeNodeInit(treeNode, vtbl, true, (void *)tree, maxNumChildren, BF_SIZE_BAD_VALUE, 0);
}

void bfTreeNodeInit(BfTreeNode *treeNode, BfTreeNodeVtable *vtbl, bool isRoot,
                    void *parent, BfSize maxNumChildren, BfSize index,
                    BfSize depth) {
  BF_ERROR_BEGIN();

  treeNode->vtbl = vtbl;
  treeNode->isRoot = isRoot;
  treeNode->index = index;
  treeNode->parent = parent;
  treeNode->maxNumChildren = maxNumChildren;
  treeNode->depth = depth;

  /* Allocate space for index offsets */
  treeNode->offset = bfMemAlloc(maxNumChildren + 1, sizeof(BfSize));
  if (treeNode->offset == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  /* Fill index offsets with bad values */
  for (BfSize i = 0; i <= maxNumChildren; ++i)
    treeNode->offset[i] = BF_SIZE_BAD_VALUE;

  /* Make space for child pointers, all `NULL` initially */
  treeNode->child = bfMemAllocAndZero(maxNumChildren, sizeof(BfTreeNode *));
  if (treeNode->child == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END() {
    bfTreeNodeDeinit(treeNode);
  }
}

void bfTreeNodeDeinit(BfTreeNode *treeNode) {
  for (BfSize i = 0; i < treeNode->maxNumChildren; ++i)
    if (treeNode->child[i])
      bfTreeNodeDelete(&treeNode->child[i]);

  treeNode->vtbl = NULL;
  treeNode->isRoot = false;
  treeNode->index = BF_SIZE_BAD_VALUE;
  treeNode->parent = NULL;
  treeNode->maxNumChildren = BF_SIZE_BAD_VALUE;
  treeNode->depth = BF_SIZE_BAD_VALUE;

  bfMemFree(treeNode->offset);
  treeNode->offset = NULL;

  bfMemFree(treeNode->child);
  treeNode->child = NULL;
}

void bfTreeNodeDealloc(BfTreeNode **treeNode) {
  bfMemFree(*treeNode);
  *treeNode = NULL;
}

bool bfTreeNodeInstanceOf(BfTreeNode const *treeNode, BfType type) {
  return bfTypeDerivedFrom(bfTreeNodeGetType(treeNode), type);
}

static void getMaxDepth(BfTree const *tree, BfTreeNode const *treeNode, BfSize *maxDepth) {
  (void)tree;
  BfSize depth = treeNode->depth;
  *maxDepth = depth > *maxDepth ? depth : *maxDepth;
}

BfSize bfTreeNodeGetMaxDepthBelow(BfTreeNode const *node) {
  BfSize maxDepth = bfTreeNodeGetDepth(node);
  bfTreeMapConst(bfTreeNodeGetTreeConst(node), node,
                 BF_TREE_TRAVERSAL_LR_LEVEL_ORDER,
                 (BfTreeMapConstFunc)getMaxDepth, &maxDepth);
  return maxDepth;
}

BfTreeNode *bfTreeNodeGetParent(BfTreeNode *node) {
  return node->parent;
}

BfTreeNode const *bfTreeNodeGetParentConst(BfTreeNode const *node) {
  return node->parent;
}

BfSize bfTreeNodeGetNumChildren(BfTreeNode const *node) {
  BfSize numChildren = 0;
  for (BfSize i = 0; i < node->maxNumChildren; ++i)
    numChildren += node->child[i] != NULL;
  return numChildren;
}

BfSize bfTreeNodeGetDepth(BfTreeNode const *node) {
  return node->depth;
}

BfSize bfTreeNodeGetNumPoints(BfTreeNode const *node) {
  if (bfTreeNodeIsLeaf(node)) {
    BfSize const *parentOffset = bfTreeNodeGetParentConst(node)->offset;
    return parentOffset[node->index + 1] - parentOffset[node->index];
  } else {
    return node->offset[node->maxNumChildren] - node->offset[0];
  }
}

bool bfTreeNodeIsRoot(BfTreeNode const *node) {
  return node->isRoot;
}

bool bfTreeNodeIsLeaf(BfTreeNode const *node) {
  for (BfSize i = 0; i < node->maxNumChildren; ++i)
    if (node->child[i])
      return false;
  return true;
}

BfSize bfTreeNodeGetFirstIndex(BfTreeNode const *treeNode) {
  if (bfTreeNodeIsLeaf(treeNode)) {
    BfTreeNode const *parent = bfTreeNodeGetParentConst(treeNode);
    return parent->offset[treeNode->index];
  } else {
    return treeNode->offset[0];
  }
}

BfSize bfTreeNodeGetLastIndex(BfTreeNode const *treeNode) {
  if (bfTreeNodeIsLeaf(treeNode)) {
    BfTreeNode const *parent = bfTreeNodeGetParentConst(treeNode);
    return parent->offset[treeNode->index + 1];
  } else {
    return treeNode->offset[treeNode->maxNumChildren];
  }
}

BfSize const *bfTreeNodeGetIndexPtrConst(BfTreeNode const *treeNode, BfTree const *tree) {
  if (tree == NULL)
    tree = bfTreeNodeGetTreeConst(treeNode);
  BfSize index = bfTreeNodeGetFirstIndex(treeNode);
  return &tree->perm->index[index];
}

BfTree *bfTreeNodeGetTree(BfTreeNode *node) {
  while (!bfTreeNodeIsRoot(node))
    node = bfTreeNodeGetParent(node);
  return (BfTree *)bfTreeNodeGetParent(node);
}

BfTree const *bfTreeNodeGetTreeConst(BfTreeNode const *node) {
  while (!node->isRoot)
    node = bfTreeNodeGetParentConst(node);
  return (BfTree const *)bfTreeNodeGetParentConst(node);
}

bool bfTreeNodeIsEmpty(BfTreeNode const *treeNode) {
  return bfTreeNodeGetNumPoints(treeNode) == 0;
}

bool bfTreeNodeIsDescendant(BfTreeNode const *treeNode, BfTreeNode const *otherTreeNode) {
  BF_ERROR_BEGIN();

  bool isDescendant = false;

  if (treeNode == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (otherTreeNode == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  if (bfTreeNodeGetType(treeNode) != bfTreeNodeGetType(otherTreeNode))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  isDescendant = treeNode == otherTreeNode;
repeat:
  if (!isDescendant && !treeNode->isRoot) {
    treeNode = treeNode->parent;
    isDescendant = treeNode == otherTreeNode;
    goto repeat;
  }

  BF_ERROR_END() {
    isDescendant = false;
  }

  return isDescendant;
}

BfTreeNode *bfTreeNodeGetChild(BfTreeNode *treeNode, BfSize i) {
  BF_ERROR_BEGIN();

  BfTreeNode *child = NULL;

  if (i >= treeNode->maxNumChildren)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  child = treeNode->child[i];

  BF_ERROR_END() {
    BF_DIE();
  }

  return child;
}

BfSize bfTreeNodeGetMaxNumChildren(BfTreeNode const *treeNode) {
  return treeNode->maxNumChildren;
}

bool bfTreeNodeHasChild(BfTreeNode const *treeNode, BfSize i) {
  return treeNode->child[i] != NULL;
}

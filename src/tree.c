#include <bf/tree.h>

#include <bf/array.h>
#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
#include <bf/node_span.h>
#include <bf/quadtree.h>
#include <bf/tree_iter_post_order.h>
#include <bf/tree_level_iter.h>
#include <bf/tree_node.h>
#include <bf/util.h>

#include "macros.h"

/** Interface: Tree */

BfTreeVtable TreeVtable;

void bfTreeCopyInto(BfTree *tree, BfTree *dstTree) {
  tree->vtable->CopyInto(tree, dstTree);
}

void bfTreeDelete(BfTree **tree) {
  (*tree)->vtable->Delete(tree);
}

BfType bfTreeGetType(BfTree const *tree) {
  return tree->vtable->GetType(tree);
}

/** Implementation: Tree */

BfTree *bfTreeAlloc(void) {
  BF_ERROR_BEGIN();

  BfTree *tree = bfMemAlloc(1, sizeof(BfTree));
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return tree;
}

BfTree *bfTreeNewFromNodeType(BfType type) {
  switch (type) {
  case BF_TYPE_QUADTREE_NODE:
    return bfQuadtreeToTree(bfQuadtreeNew());
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

BfTree *bfTreeNewFromNodeSpan(BfNodeSpan *nodeSpan, BfPerm **perm) {
  BF_ERROR_BEGIN();

  BfType nodeType = bfNodeSpanGetNodeType(nodeSpan);

  BfTree *tree = bfTreeNewFromNodeType(nodeType);
  HANDLE_ERROR();

  bfTreeInitFromNodeSpan(tree, nodeSpan, perm);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return tree;
}

BfTree *bfTreeNewFromNode(BfTreeNode *node, BfPerm **perm) {
  BF_ERROR_BEGIN();

  BfTree *tree = bfTreeNewFromNodeType(bfTreeNodeGetType(node));
  HANDLE_ERROR();

  bfTreeInitFromNode(tree, node, perm);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return tree;
}

BfTree *bfTreeNewForMiddleFac(BfTree const *templateTree, BfSize p) {
  BF_ERROR_BEGIN();

  BfTree *tree = bfTreeAlloc();
  HANDLE_ERROR();

  bfTreeInitForMiddleFac(tree, templateTree, p);
  HANDLE_ERROR();

  printf("max depth: %lu\n", bfTreeGetMaxDepth(tree));

  BF_ERROR_END() {
    BF_DIE();
  }

  return tree;
}

void bfTreeInit(BfTree *tree, BfTreeVtable *vtable, BfTreeNode *root, BfSize size) {
  BF_ERROR_BEGIN();

  tree->vtable = vtable;
  tree->root = root;

  tree->perm = bfPermNewIdentity(size);
  HANDLE_ERROR();

  BF_ERROR_END()
    bfTreeDeinit(tree);
}

static void shiftOffsets(BfTree *tree, BfSize i0) {
  BF_ERROR_BEGIN();

  BfPtrArray *queue = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  BfTreeNode *root = bfTreeGetRootNode(tree);
  if (root == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  bfPtrArrayAppend(queue, root);
  while (!bfPtrArrayIsEmpty(queue)) {
    BfTreeNode *node = bfPtrArrayPopFirst(queue);

    for (BfSize i = 0; i <= node->maxNumChildren; ++i) {
      if (node->offset[i] < i0)
        RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);
      node->offset[i] -= i0;
    }

    for (BfSize i = 0; i < node->maxNumChildren; ++i) {
      if (node->child[i] != NULL) {
        bfPtrArrayAppend(queue, node->child[i]);
        HANDLE_ERROR();
      }
    }
  }

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfTreeInitFromNodeSpan(BfTree *tree, BfNodeSpan *nodeSpan, BfPerm **perm) {
  BF_ERROR_BEGIN();

  if (perm == NULL || *perm != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfTree *nodeSpanTree = bfNodeSpanGetTree(nodeSpan);
  if (nodeSpanTree == NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  /* Initially just set `tree` to be an exact copy of
   * `nodeSpanTree`... We'll modify it as we go. */
  bfTreeCopyInto(nodeSpanTree, tree);

  BfTreeNode *nodeSpanRoot = bfNodeSpanGetRoot(nodeSpan);
  if (nodeSpanRoot == NULL)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfSize i0 = bfTreeNodeGetFirstIndex(nodeSpanRoot);
  BfSize i1 = bfTreeNodeGetLastIndex(nodeSpanRoot);

  *perm = bfPermGetRangeView(tree->perm, i0, i1);
  HANDLE_ERROR();

  tree->perm = bfPermNewIdentity(i1 - i0);
  HANDLE_ERROR();

  tree->root = bfTreeNodeCopy(nodeSpanRoot);
  HANDLE_ERROR();

  /* Mark this as the root of a tree: */
  tree->root->index = BF_SIZE_BAD_VALUE;
  tree->root->isRoot = true;
  tree->root->depth = 0;

  /* We'll use a queue to traverse the subtree rooted at `tree->root`,
   * doing two things:
   *
   * - duplicate the nodes we encounter
   * - stop short at the nodes in `nodeSpan`---we'll terminate the
   *   tree at these nodes, turning them into the new leaves of our
   *   "node span tree" */

  BfPtrArray *queue = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  bfPtrArrayAppend(queue, tree->root);
  HANDLE_ERROR();

  while (!bfPtrArrayIsEmpty(queue)) {
    BfTreeNode *parent = bfPtrArrayPopFirst(queue);

    /* If this node is in `nodeSpan`, erase its children, turning it
     * into a leaf node: */
    if (bfNodeSpanContainsNodeWithSameRange(nodeSpan, parent)) {
      for (BfSize i = 0; i < parent->maxNumChildren; ++i)
        parent->child[i] = NULL;
      continue;
    }

    /* ... otherwise, copy and update the node, then push it onto the
     * queue: */
    for (BfSize i = 0; i < parent->maxNumChildren; ++i) {
      if (parent->child[i] == NULL)
        continue;

      parent->child[i] = bfTreeNodeCopy(parent->child[i]);
      parent->child[i]->parent = parent;
      parent->child[i]->depth = parent->depth + 1;

      bfPtrArrayAppend(queue, parent->child[i]);
    }
  }

  shiftOffsets(tree, i0);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfTreeInitFromNode(BfTree *tree, BfTreeNode *node, BfPerm **perm) {
  BF_ERROR_BEGIN();

  if (perm == NULL || *perm != NULL)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfTree *nodeTree = bfTreeNodeGetTree(node);

  bfTreeCopyInto(nodeTree, tree);

  BfSize i0 = bfTreeNodeGetFirstIndex(node);
  BfSize i1 = bfTreeNodeGetLastIndex(node);

  *perm = bfPermGetRangeView(tree->perm, i0, i1);
  HANDLE_ERROR();

  tree->perm = bfPermNewIdentity(i1 - i0);
  HANDLE_ERROR();

  tree->root = bfTreeNodeCopy(node);
  HANDLE_ERROR();

  /* Mark this as the root of a tree: */
  tree->root->index = BF_SIZE_BAD_VALUE;
  tree->root->isRoot = true;
  tree->root->depth = 0;

  BfPtrArray *queue = bfPtrArrayNewWithDefaultCapacity();
  HANDLE_ERROR();

  bfPtrArrayAppend(queue, tree->root);
  HANDLE_ERROR();

  while (!bfPtrArrayIsEmpty(queue)) {
    BfTreeNode *parent = bfPtrArrayPopFirst(queue);
    for (BfSize i = 0; i < parent->maxNumChildren; ++i) {
      if (parent->child[i] == NULL)
        continue;

      parent->child[i] = bfTreeNodeCopy(parent->child[i]);
      parent->child[i]->parent = parent;
      parent->child[i]->depth = parent->depth + 1;

      bfPtrArrayAppend(queue, parent->child[i]);
    }
  }

  shiftOffsets(tree, i0);

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfTreeInitForMiddleFac(BfTree *tree, BfTree const *templateTree, BfSize p) {
  BF_ERROR_BEGIN();

  tree->vtable = &TreeVtable;

  tree->root = bfTreeNodeAlloc();
  HANDLE_ERROR();

  BfTreeNodeVtable *vtable = bfTreeNodeGetVtable();

  struct {
    BfTreeNode const *templateNode;
    BfTreeNode *node;
    void *parent;
  } nodePair;

  nodePair.templateNode = bfTreeGetRootNodeConst(templateTree);
  nodePair.node = tree->root;
  nodePair.parent = NULL;

  BfArray *queue = bfArrayNewEmpty(sizeof(nodePair));
  HANDLE_ERROR();

  bfArrayAppend(queue, &nodePair);
  HANDLE_ERROR();

  /* We make a first pass over the template tree to instantiate the
   * new tree and to the nodes' parents and children: */

  while (!bfArrayIsEmpty(queue)) {
    bfArrayPopFirst(queue, &nodePair);
    HANDLE_ERROR();

    BfTreeNode const *templateNode = nodePair.templateNode;
    BfTreeNode *node = nodePair.node;
    void *parent = nodePair.parent;

    if (parent == NULL)
      bfTreeNodeInitRoot(node, vtable, tree, templateNode->maxNumChildren);
    else
      bfTreeNodeInit(
        node, vtable, false, parent, templateNode->maxNumChildren,
        templateNode->index, templateNode->depth);
    HANDLE_ERROR();

    for (BfSize i = 0; i < templateNode->maxNumChildren; ++i) {
      if (templateNode->child[i] == NULL)
        continue;

      nodePair.templateNode = templateNode->child[i];
      nodePair.parent = node;
      nodePair.node = node->child[i] = bfTreeNodeAlloc();
      HANDLE_ERROR();

      bfArrayAppend(queue, &nodePair);
    }
  }

  BfTreeIterPostOrder iter;
  bfTreeIterPostOrderInit(&iter, tree);
  HANDLE_ERROR();

  while (!bfTreeIterPostOrderIsDone(&iter)) {
    BfTreeNode *node = bfTreeIterPostOrderGetCurrentNode(&iter);

    bool isRoot = bfTreeNodeIsRoot(node);
    BfSize index = node->index;

    if (!isRoot && index == 0)
      bfTreeNodeGetParent(node)->offset[0] = 0;

    BfSize numPoints = BF_SIZE_BAD_VALUE;

    if (bfTreeNodeIsLeaf(node)) {
      /* TODO: handle this case... hmm, not sure if the current tree
       * structure will even support it! That's not good :-( */
      if (isRoot)
        RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

      numPoints = p;
    } else {
      /* TODO: handle this case... */
      for (BfSize i = 0; i <= node->maxNumChildren; ++i)
        BF_ASSERT(BF_SIZE_OK(node->offset[i]));

      bfSizeRunningSum(node->maxNumChildren + 1, node->offset);

      numPoints = node->offset[node->maxNumChildren];
    }

    if (!isRoot)
      bfTreeNodeGetParent(node)->offset[index + 1] = numPoints;

    bfTreeIterPostOrderNext(&iter);
  }

  bfTreeIterPostOrderDeinit(&iter);

  tree->perm = bfPermNewIdentity(bfTreeNodeGetNumPoints(tree->root));

  BF_ERROR_END() {
    BF_DIE();
  }

  bfArrayDeinitAndDealloc(&queue);
}

void bfTreeDeinit(BfTree *tree) {
  tree->vtable = NULL;
  bfPermDeinitAndDealloc(&tree->perm);
  bfTreeNodeDelete(&tree->root);
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

BfPerm *bfTreeGetPerm(BfTree *tree) {
  return tree->perm;
}

BfPerm const *bfTreeGetPermConst(BfTree const *tree) {
  return tree->perm;
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
    for (BfSize i = 0; i < bfPtrArraySize(iter.levelNodes); ++i) {
      BfTreeNode *currentNode = bfPtrArrayGet(iter.levelNodes, i);
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
    for (BfSize i = 0; i < bfPtrArraySize(iter.levelNodes); ++i) {
      BfTreeNode const *currentNode = bfPtrArrayGet(iter.levelNodes, i);
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

  BfTreeNode *node = NULL;

  BfPtrArray *levelNodes = bfTreeGetLevelPtrArray(tree, depth);
  HANDLE_ERROR();

  node = bfPtrArrayGet(levelNodes, nodeIndex);
  HANDLE_ERROR();

  BF_ERROR_END() {
    node = NULL;
  }

  bfPtrArrayDeinitAndDealloc(&levelNodes);

  return node;
}

BfPtrArray *bfTreeGetLevelPtrArray(BfTree *tree, BfSize depth) {
  BF_ERROR_BEGIN();

  BfSize maxDepth = bfTreeGetMaxDepth(tree);

  if (depth > maxDepth)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfTreeLevelIter iter;
  bfTreeLevelIterInit(&iter, BF_TREE_TRAVERSAL_LR_LEVEL_ORDER, tree->root);
  HANDLE_ERROR();

  for (BfSize _ = 0; _ < depth; ++_)
    bfTreeLevelIterNext(&iter);

  BfPtrArray *levelNodes = bfPtrArrayCopy(iter.levelNodes);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  bfTreeLevelIterDeinit(&iter);

#if BF_DEBUG
  for (BfSize k = 0; k < bfPtrArraySize(levelNodes); ++k)
    BF_ASSERT(bfPtrArrayGet(levelNodes, k) != NULL);
#endif

  return levelNodes;
}

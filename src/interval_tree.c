#include <bf/interval_tree.h>

#include <bf/assert.h>
#include <bf/const.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/interval_tree_node.h>
#include <bf/mem.h>
#include <bf/tree.h>

#include "macros.h"

/** Interface: Tree -> IntervalTree */

static BfTreeVtable TreeVtable = {
  .GetType = (__typeof__(&bfTreeGetType))bfIntervalTreeGetType
};

BfType bfIntervalTreeGetType(BfTree const *tree) {
  (void)tree;
  return BF_TYPE_INTERVAL_TREE;
}

/** Downcasting: Tree -> IntervalTree */

BfIntervalTree *bfTreeToIntervalTree(BfTree *tree) {
  if (!bfTreeInstanceOf(tree, BF_TYPE_INTERVAL_TREE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfIntervalTree *)tree;
  }
}

BfIntervalTree const *bfTreeConstToIntervalTreeConst(BfTree const *tree) {
  if (!bfTreeInstanceOf(tree, BF_TYPE_INTERVAL_TREE)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfIntervalTree const *)tree;
  }
}

/** Upcasting: IntervalTree -> Tree */

BfTree *bfIntervalTreeToTree(BfIntervalTree *intervalTree) {
  return &intervalTree->super;
}

/** Implementation: IntervalTree */

BfIntervalTree *bfIntervalTreeNew() {
  BEGIN_ERROR_HANDLING();

  BfIntervalTree *intervalTree = bfMemAlloc(1, sizeof(BfIntervalTree));
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    intervalTree = NULL;

  return intervalTree;
}

void bfIntervalTreeInitEmpty(BfIntervalTree *intervalTree, BfReal a, BfReal b, BfSize k, BfSize maxDepth) {
  BEGIN_ERROR_HANDLING();

  BfIntervalTreeNode *root = bfIntervalTreeNodeNew();
  HANDLE_ERROR();

  bfTreeInit(&intervalTree->super, &TreeVtable, bfIntervalTreeNodeToTreeNode(root), 0);
  HANDLE_ERROR();

  intervalTree->points = NULL;

  bfIntervalTreeNodeInitEmptyRoot(root, intervalTree, a, b, k, maxDepth);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {
    // bfIntervalTreeDeinit(intervalTree);
  }
}

static void recursivelySiftNodes(BfIntervalTreeNode *intervalTreeNode,
                                 BfPoints1 const *points, BfPerm *perm) {
  BfTreeNode *treeNode = bfIntervalTreeNodeToTreeNode(intervalTreeNode);

  BfSize K = treeNode->maxNumChildren;

  BfTreeNode *parent = bfTreeNodeGetParent(treeNode);

  BfSize i0 = treeNode->isRoot ? 0 : parent->offset[treeNode->index];
  BfSize i1 = treeNode->isRoot ? points->size : parent->offset[treeNode->index + 1];

  BfSize i = i0;

  BfReal const *x = points->data;
  BfSize *P = perm->index;

  treeNode->offset[0] = i0;

  for (BfSize k = 0; k < K; ++k) {
    if (treeNode->child[k] == NULL)
      continue;

    BfIntervalTreeNode *child = bfTreeNodeToIntervalTreeNode(treeNode->child[k]);

    BfReal a = child->isLeftmost ? -BF_INFINITY : child->a;
    BfReal b = child->isRightmost ? BF_INFINITY : child->b;

    /* Skip over leading points which are already in position */
    while (i < i1 && a <= x[P[i]] && x[P[i]] < b)
      ++i;

    BfSize j = i + 1;
    while (i < i1 && j < i1) {
      BF_ASSERT(x[P[i]] < a || b <= x[P[i]]);

      /* Find the next point that's in the current child interval */
      while (j < i1 && (x[P[j]] < a || b <= x[P[j]]))
        ++j;

      if (j == i1)
        break;

      SWAP(P[i], P[j]);
      ++i;
      ++j;

      while (i < i1 && a <= x[P[i]] && x[P[i]] < b)
        ++i;
    }

    treeNode->offset[k + 1] = i;

    recursivelySiftNodes(child, points, perm);
  }
}

void bfIntervalTreeSetPoints(BfIntervalTree *intervalTree, BfPoints1 const *points, bool rebuildTree) {
  BEGIN_ERROR_HANDLING();

  if (rebuildTree)
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  intervalTree->points = points;

  bfPermDeinit(&intervalTree->super.perm);

  intervalTree->super.perm = bfPermIdentity(points->size);
  HANDLE_ERROR();

  BfIntervalTreeNode *root = bfTreeNodeToIntervalTreeNode(intervalTree->super.root);

  recursivelySiftNodes(root, points, &intervalTree->super.perm);
  HANDLE_ERROR();

  END_ERROR_HANDLING() {}
}

#include <bf/tree_iter_post_order.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/tree.h>
#include <bf/tree_node.h>

/* We use this type to wrap a TreeNode and an index to the next of its
 * children that should be visited in the post order traversal. */
typedef struct {
  BfTreeNode *treeNode;
  BfSize nextChild;
} StackNode;

/* Make a new StackNode from a TreeNode and initialize the index to
 * the first of its children which should be visited. */
StackNode *makeStackNode(BfTreeNode *treeNode) {
  BEGIN_ERROR_HANDLING();

  StackNode *stackNode = malloc(sizeof(StackNode));
  if (stackNode == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  stackNode->treeNode = treeNode;

  stackNode->nextChild = BF_SIZE_BAD_VALUE;
  for (BfSize i = 0; i < treeNode->maxNumChildren; ++i) {
    if (treeNode->child[i] != NULL) {
      stackNode->nextChild = i;
      break;
    }
  }

  END_ERROR_HANDLING() {}

  return stackNode;
}

/** Interface: TreeIter */

BfType bfTreeIterPostOrderGetType(BfTreeIterPostOrder const *iter) {
  (void)iter;
  return BF_TYPE_TREE_ITER_POST_ORDER;
}

BfTreeNode *bfTreeIterPostOrderGetCurrentNode(BfTreeIterPostOrder *iter) {
  BEGIN_ERROR_HANDLING();

  StackNode *stackNode = NULL;
  BfTreeNode *treeNode = NULL;

  if (bfPtrArrayIsEmpty(&iter->stack))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  bfPtrArrayGetLast(&iter->stack, (BfPtr *)&stackNode);

  if (stackNode->nextChild != BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  treeNode = stackNode->treeNode;

  END_ERROR_HANDLING() {
    treeNode = NULL;
  }

  return treeNode;
}

bool bfTreeIterPostOrderIsDone(BfTreeIterPostOrder const *iter) {
  return bfPtrArrayIsEmpty(&iter->stack);
}

void bfTreeIterPostOrderNext(BfTreeIterPostOrder *iter) {
  (void)iter;
  assert(false);
}

#define _(FUNC_NAME)                                                    \
  .FUNC_NAME = (__typeof__(&bfTreeIter##FUNC_NAME))bfTreeIterPostOrder##FUNC_NAME
static BfTreeIterVtable TREE_ITER_VTABLE = {
  _(GetType),
  _(GetCurrentNode),
  _(IsDone),
  _(Next)
};
#undef _

/** Upcasting: TreeIterPostOrder -> TreeIter */

BfTreeIter *bfTreeIterPostOrderToTreeIter(BfTreeIterPostOrder *iter) {
  return &iter->super;
}

/** Implementation: TreeIterPostOrder */

BfTreeIterPostOrder *bfTreeIterPostOrderNew() {
  BEGIN_ERROR_HANDLING();

  BfTreeIterPostOrder *iter = malloc(sizeof(BfTreeIterPostOrder));
  if (iter == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING()
    iter = NULL;

  return iter;
}

void bfTreeIterPostOrderInit(BfTreeIterPostOrder *iter, BfTree const *tree) {
  BEGIN_ERROR_HANDLING();

  bfTreeIterInit(&iter->super, &TREE_ITER_VTABLE, tree);

  iter->stack = bfGetUninitializedPtrArray();
  bfInitPtrArray(&iter->stack, BF_ARRAY_DEFAULT_CAPACITY);
  HANDLE_ERROR();

  StackNode *stackNode = makeStackNode(tree->root);
  HANDLE_ERROR();

  bfPtrArrayAppend(&iter->stack, stackNode);
  HANDLE_ERROR();

  while (stackNode->nextChild != BF_SIZE_BAD_VALUE) {
    stackNode = makeStackNode(stackNode->treeNode->child[stackNode->nextChild]);
    HANDLE_ERROR();

    bfPtrArrayAppend(&iter->stack, stackNode);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {}
}

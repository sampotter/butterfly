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

void incStackNodeChild(StackNode *stackNode) {
  assert(stackNode->nextChild != BF_SIZE_BAD_VALUE);

  BfTreeNode *treeNode = stackNode->treeNode;

  /* Seek to the next non-NULL child */
  BfSize i = stackNode->nextChild + 1;
  while (i < treeNode->maxNumChildren) {
    if (treeNode->child[i] != NULL) {
      stackNode->nextChild = i;
      break;
    }
    ++i;
  }

  /* If we make it to the end without finding a non-NULL child, we've
   * visited all of `treeNode`'s children and can mark this stack node
   * as done by setting `nextChild` to `BF_SIZE_BAD_VALUE` */
  if (i == treeNode->maxNumChildren)
    stackNode->nextChild = BF_SIZE_BAD_VALUE;
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

static void pushNewChildNodes(BfPtrArray *stack, StackNode *stackNode) {
  BEGIN_ERROR_HANDLING();

  while (stackNode->nextChild != BF_SIZE_BAD_VALUE) {
    stackNode = makeStackNode(
      stackNode->treeNode->child[stackNode->nextChild]);
    HANDLE_ERROR();

    bfPtrArrayAppend(stack, stackNode);
    HANDLE_ERROR();
  }

  END_ERROR_HANDLING() {}
}

void bfTreeIterPostOrderNext(BfTreeIterPostOrder *iter) {
  BEGIN_ERROR_HANDLING();

  if (bfPtrArrayIsEmpty(&iter->stack))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Start by popping the current stack node off the stack used to
   * drive the iterator */
  bfPtrArrayPopLast(&iter->stack);

  /* If the stack still has nodes on it, we need to inspect the next
   * stack node to decide what to do. */
  if (!bfPtrArrayIsEmpty(&iter->stack)) {
    StackNode *stackNode = NULL;
    bfPtrArrayGetLast(&iter->stack, (BfPtr *)&stackNode);

    /* If this node has more children to visit beneath it, we
     * increment the current stack node's child pointer and
     * recursively push the next one (and all subsequent leftmost
     * children beneath it) onto the stack */
    if (stackNode->nextChild != BF_SIZE_BAD_VALUE) {
      incStackNodeChild(stackNode);

      pushNewChildNodes(&iter->stack, stackNode);
      HANDLE_ERROR();
    }

    /* ... otherwise, `stackNode` is the next node to visit */
  }

  END_ERROR_HANDLING() {
    assert(false);
  }
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

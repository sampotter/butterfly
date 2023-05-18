#include <bf/tree_iter_post_order.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>
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
  BF_ERROR_BEGIN();

  StackNode *stackNode = bfMemAlloc(1, sizeof(StackNode));
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

  BF_ERROR_END() {}

  return stackNode;
}

void deleteStackNode(StackNode *stackNode) {
  bfMemFree(stackNode);
}

void incStackNodeChild(StackNode *stackNode) {
  BF_ASSERT(stackNode->nextChild != BF_SIZE_BAD_VALUE);

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

void bfTreeIterPostOrderDelete(BfTreeIterPostOrder **iter) {
  bfTreeIterPostOrderDeinit(*iter);
  bfTreeIterPostOrderDealloc(iter);
}

BfTreeNode *bfTreeIterPostOrderGetCurrentNode(BfTreeIterPostOrder *iter) {
  BF_ERROR_BEGIN();

  StackNode *stackNode = NULL;
  BfTreeNode *treeNode = NULL;

  if (bfPtrArrayIsEmpty(&iter->stack))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  bfPtrArrayGetLast(&iter->stack, (BfPtr *)&stackNode);

  if (stackNode->nextChild != BF_SIZE_BAD_VALUE)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  treeNode = stackNode->treeNode;

  BF_ERROR_END() {
    treeNode = NULL;
  }

  return treeNode;
}

bool bfTreeIterPostOrderIsDone(BfTreeIterPostOrder const *iter) {
  return bfPtrArrayIsEmpty(&iter->stack);
}

static void pushNewChildNodes(BfPtrArray *stack, StackNode *stackNode) {
  BF_ERROR_BEGIN();

  while (stackNode->nextChild != BF_SIZE_BAD_VALUE) {
    stackNode = makeStackNode(
      stackNode->treeNode->child[stackNode->nextChild]);
    HANDLE_ERROR();

    bfPtrArrayAppend(stack, stackNode);
    HANDLE_ERROR();
  }

  BF_ERROR_END() {}
}

void bfTreeIterPostOrderNext(BfTreeIterPostOrder *iter) {
  BF_ERROR_BEGIN();

  if (bfPtrArrayIsEmpty(&iter->stack))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  /* Start by popping the current stack node off the stack used to
   * drive the iterator */
  StackNode *stackNode = bfPtrArrayPopLast(&iter->stack);
  deleteStackNode(stackNode);

  /* If the stack still has nodes on it, we need to inspect the next
   * stack node to decide what to do. */
  if (!bfPtrArrayIsEmpty(&iter->stack)) {
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

  BF_ERROR_END() {
    BF_DIE();
  }
}

#define _(FUNC_NAME)                                                    \
  .FUNC_NAME = (__typeof__(&bfTreeIter##FUNC_NAME))bfTreeIterPostOrder##FUNC_NAME
static BfTreeIterVtable TREE_ITER_VTABLE = {
  _(GetType),
  _(Delete),
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
  BF_ERROR_BEGIN();

  BfTreeIterPostOrder *iter = bfMemAlloc(1, sizeof(BfTreeIterPostOrder));
  if (iter == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  BF_ERROR_END()
    iter = NULL;

  return iter;
}

void bfTreeIterPostOrderInit(BfTreeIterPostOrder *iter, BfTree const *tree) {
  BF_ERROR_BEGIN();

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

  BF_ERROR_END() {}
}

void bfTreeIterPostOrderDeinit(BfTreeIterPostOrder *iter) {
  for (BfSize i = 0; i < bfPtrArraySize(&iter->stack); ++i) {
    StackNode *stackNode = bfPtrArrayGet(&iter->stack, i);
    deleteStackNode(stackNode);
  }
  bfPtrArrayDeinit(&iter->stack);

  bfTreeIterDeinit(&iter->super);
}

void bfTreeIterPostOrderDealloc(BfTreeIterPostOrder **iter) {
  bfMemFree(*iter);
  *iter = NULL;
}

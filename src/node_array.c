#include <bf/node_array.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mem.h>

BfConstNodeArray *bfConstNodeArrayNewWithDefaultCapacity() {
  BF_ERROR_BEGIN();

  BfConstNodeArray *nodeArray = bfMemAlloc(1, sizeof(BfConstNodeArray));
  HANDLE_ERROR();

  bfConstNodeArrayInitWithDefaultCapacity(nodeArray);

  BF_ERROR_END() {
    BF_DIE();
  }

  return nodeArray;
}

void bfConstNodeArrayInitWithDefaultCapacity(BfConstNodeArray *nodes) {
  bfConstPtrArrayInitWithDefaultCapacity(&nodes->ptrArray);
}

void bfConstNodeArrayDeinit(BfConstNodeArray *nodes) {
  bfConstPtrArrayDeinit(&nodes->ptrArray);
}

void bfConstNodeArrayDealloc(BfConstNodeArray **nodes) {
  bfMemFree(*nodes);
  *nodes = NULL;
}

void bfConstNodeArrayDeinitAndDealloc(BfConstNodeArray **nodes) {
  bfConstNodeArrayDeinit(*nodes);
  bfConstNodeArrayDealloc(nodes);
}

void bfConstNodeArrayAppend(BfConstNodeArray *nodes, BfTreeNode const *node) {
  bfConstPtrArrayAppend(&nodes->ptrArray, node);
}

BfSize bfConstNodeArraySize(BfConstNodeArray const *nodes) {
  return bfConstPtrArraySize(&nodes->ptrArray);
}

bool bfConstNodeArrayIsEmpty(BfConstNodeArray const *nodes) {
  return bfConstPtrArrayIsEmpty(&nodes->ptrArray);
}

BfTreeNode const *bfConstNodeArrayGet(BfConstNodeArray const *nodes, BfSize i) {
  return bfConstPtrArrayGet(&nodes->ptrArray, i);
}

BfTreeNode const *bfConstNodeArrayGetFirst(BfConstNodeArray const *nodes) {
  return bfConstPtrArrayGetFirst(&nodes->ptrArray);
}

BfTreeNode const *bfConstNodeArrayGetLast(BfConstNodeArray const *nodes) {
  return bfConstPtrArrayGetLast(&nodes->ptrArray);
}

void bfConstNodeArrayExtend(BfConstNodeArray *nodes, BfConstNodeArray const *otherNodes) {
  bfConstPtrArrayExtend(&nodes->ptrArray, &otherNodes->ptrArray);
}

bool bfConstNodeArrayIsContiguous(BfConstNodeArray const *nodes) {
  BfSize numNodes = bfConstNodeArraySize(nodes);
  if (numNodes <= 1)
    return true;

  BfTreeNode const *node = bfConstNodeArrayGetFirst(nodes);
  for (BfSize i = 1; i < numNodes; ++i) {
    BfSize i1 = bfTreeNodeGetLastIndex(node);
    node = bfConstNodeArrayGet(nodes, i);
    if (bfTreeNodeGetFirstIndex(node) != i1)
      return false;
  }

  return true;
}

bool bfConstNodeArrayIsSameSpan(BfConstNodeArray const *nodes, BfConstNodeArray const *otherNodes) {
  BfSize numNodes = bfConstNodeArraySize(nodes);
  if (numNodes <= 1)
    return true;

  if (numNodes != bfConstNodeArraySize(otherNodes))
    return false;

  BfTreeNode const *node = bfConstNodeArrayGetFirst(nodes);
  BfTreeNode const *otherNode = bfConstNodeArrayGetFirst(otherNodes);

  if (bfTreeNodeGetFirstIndex(node) != bfTreeNodeGetFirstIndex(otherNode))
    return false;

  for (BfSize i = 1; i < numNodes; ++i) {
    BfSize i1 = bfTreeNodeGetLastIndex(node);
    if (i1 != bfTreeNodeGetLastIndex(otherNode))
      return false;

    node = bfConstNodeArrayGet(nodes, i);
    if (bfTreeNodeGetFirstIndex(node) != i1)
      return false;

    otherNode = bfConstNodeArrayGet(otherNodes, i);
    if (bfTreeNodeGetFirstIndex(otherNode) != i1)
      return false;
  }

  return true;
}

bool nodesHaveSameFirstIndex(BfConstNodeArray const *nodes) {
  BF_ERROR_BEGIN();

  bool sameFirstIndex = true;

  BfSize numNodes = bfConstNodeArraySize(nodes);
  if (numNodes == 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfTreeNode const *node = bfConstNodeArrayGetFirst(nodes);
  BfSize i0 = bfTreeNodeGetFirstIndex(node);
  for (BfSize k = 1; k < numNodes; ++k) {
    node = bfConstNodeArrayGet(nodes, k);
    if (bfTreeNodeGetFirstIndex(node) != i0) {
      sameFirstIndex = false;
      break;
    }
  }

  BF_ERROR_END() {}

  return sameFirstIndex;
}

BfSize getMaxLastIndexForRowNodes(BfConstNodeArray const *nodes,
                                  BfTreeNode const **argmaxNodePtr) {
  BF_ERROR_BEGIN();

  BfSize i1Max = BF_SIZE_BAD_VALUE;
  BfTreeNode const *argmaxNode = NULL;

  BfSize numNodes = bfConstNodeArraySize(nodes);
  if (numNodes == 0)
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  argmaxNode = bfConstNodeArrayGetFirst(nodes);
  i1Max = bfTreeNodeGetLastIndex(argmaxNode);
  for (BfSize k = 1; k < numNodes; ++k) {
    BfTreeNode const *node = bfConstNodeArrayGet(nodes, k);
    BfSize i1 = bfTreeNodeGetLastIndex(node);
    if (i1 > i1Max) {
      i1Max = i1;
      argmaxNode = node;
    }
  }

  BF_ERROR_END() {
    i1Max = BF_SIZE_BAD_VALUE;
    argmaxNode = NULL;
  }

  if (argmaxNodePtr != NULL)
    *argmaxNodePtr = argmaxNode;

  return i1Max;
}

BfTreeNode const *getNodeByFirstIndex(BfConstNodeArray const *nodes, BfSize i0) {
  BfTreeNode const *node = NULL;
  for (BfSize k = 0; k < bfConstNodeArraySize(nodes); ++k) {
    node = bfConstNodeArrayGet(nodes, k);
    if (bfTreeNodeGetFirstIndex(node) == i0)
      break;
  }
  return node;
}

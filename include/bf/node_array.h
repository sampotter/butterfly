#pragma once

#include "ptr_array.h"
#include "tree_node.h"

struct BfConstNodeArray {
  BfConstPtrArray ptrArray;
};

BfConstNodeArray *bfConstNodeArrayNewWithDefaultCapacity();
void bfConstNodeArrayInitWithDefaultCapacity(BfConstNodeArray *nodes);
void bfConstNodeArrayDeinit(BfConstNodeArray *nodes);
void bfConstNodeArrayDealloc(BfConstNodeArray **nodes);
void bfConstNodeArrayDeinitAndDealloc(BfConstNodeArray **nodes);
void bfConstNodeArrayAppend(BfConstNodeArray *nodes, BfTreeNode const *node);
BfSize bfConstNodeArraySize(BfConstNodeArray const *nodes);
bool bfConstNodeArrayIsEmpty(BfConstNodeArray const *nodes);
BfTreeNode const *bfConstNodeArrayGet(BfConstNodeArray const *nodes, BfSize i);
BfTreeNode const *bfConstNodeArrayGetFirst(BfConstNodeArray const *nodes);
BfTreeNode const *bfConstNodeArrayGetLast(BfConstNodeArray const *nodes);
void bfConstNodeArrayExtend(BfConstNodeArray *nodes, BfConstNodeArray const *otherNodes);

bool nodesHaveSameFirstIndex(BfConstNodeArray const *nodes);
BfSize getMaxLastIndexForRowNodes(BfConstNodeArray const *nodes, BfTreeNode const **argmaxNodePtr);
BfTreeNode const *getNodeByFirstIndex(BfConstNodeArray const *nodes, BfSize i0);

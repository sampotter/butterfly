#pragma once

#include "ptr_array.h"

struct BfNodeSpan {
  BfPtrArray *nodes;
};

BfNodeSpan *bfNodeSpanNewFromPtrArray(BfPtrArray *ptrArray, BfPolicy policy);
BfSize bfNodeSpanGetSize(BfNodeSpan const *nodeSpan);
bool bfNodeSpanIsEmpty(BfNodeSpan const *nodeSpan);
BfType bfNodeSpanGetNodeType(BfNodeSpan const *nodeSpan);
BfTree *bfNodeSpanGetTree(BfNodeSpan const *nodeSpan);
bool bfNodeSpanContains(BfNodeSpan const *nodeSpan, BfTreeNode const *node);
bool bfNodeSpanContainsNodeWithSameRange(BfNodeSpan const *nodeSpan, BfTreeNode const *node);
BfTreeNode *bfNodeSpanGetNode(BfNodeSpan const *nodeSpan, BfSize i);
BfTreeNode *bfNodeSpanGetRoot(BfNodeSpan const *nodeSpan);
BfSize bfNodeSpanGetFirstIndex(BfNodeSpan const *nodeSpan);
BfSize bfNodeSpanGetLastIndex(BfNodeSpan const *nodeSpan);
BfTreeNode *bfNodeSpanGetFirstNode(BfNodeSpan const *nodeSpan);
BfTreeNode *bfNodeSpanGetLastNode(BfNodeSpan const *nodeSpan);

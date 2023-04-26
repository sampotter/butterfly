#pragma once

#include "types.h"

struct BfTreeNodeSpan {
  BfNodeArray *nodeArray;
};

BfTreeNode const *bfTreeNodeSpanGetByFirstIndex(BfTreeNodeSpan const *span, BfSize i0);

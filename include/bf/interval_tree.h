#pragma once

#include "tree.h"

/** IntervalTree: */

struct BfIntervalTree {
  BfTree super;
};

BfType bfIntervalTreeGetType(BfTree const *tree);

/** Interface: Tree -> IntervalTree */

/** Upcasting: IntervalTree -> Tree */

BfTree *bfIntervalTreeToTree(BfIntervalTree *intervalTree);

BfIntervalTree *bfIntervalTreeNew();
void bfIntervalTreeInitEmpty(BfIntervalTree *intervalTree, BfReal a, BfReal b, BfSize k, BfSize depth);

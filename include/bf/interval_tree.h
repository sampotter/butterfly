#pragma once

#include "points.h"
#include "tree.h"

/** IntervalTree: */

struct BfIntervalTree {
  BfTree super;
  BfPoints1 const *points;
};

BfType bfIntervalTreeGetType(BfTree const *tree);

/** Interface: Tree -> IntervalTree */

BfIntervalTree *bfTreeToIntervalTree(BfTree *tree);
BfIntervalTree const *bfTreeConstToIntervalTreeConst(BfTree const *tree);

/** Upcasting: IntervalTree -> Tree */

BfTree *bfIntervalTreeToTree(BfIntervalTree *intervalTree);

BfIntervalTree *bfIntervalTreeNew();
void bfIntervalTreeInitEmpty(BfIntervalTree *intervalTree, BfReal a, BfReal b, BfSize k, BfSize depth);
void bfIntervalTreeSetPoints(BfIntervalTree *intervalTree, BfPoints1 const *points, bool rebuildTree);

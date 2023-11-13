#pragma once

#include "interval.h"
#include "points.h"
#include "tree.h"

/** IntervalTree: */

struct BfIntervalTree {
  BfTree super;
  BfPoints1 const *points;
};

BfType bfIntervalTreeGetType(BfTree const *tree);
void bfIntervalTreeDelete(BfIntervalTree **intervalTree);

/** Interface: Tree -> IntervalTree */

BfIntervalTree *bfTreeToIntervalTree(BfTree *tree);
BfIntervalTree const *bfTreeConstToIntervalTreeConst(BfTree const *tree);

/** Upcasting: IntervalTree -> Tree */

BfTree *bfIntervalTreeToTree(BfIntervalTree *intervalTree);
BfTree const *bfIntervalTreeConstToTreeConst(BfIntervalTree const *intervalTree);

BfIntervalTree *bfIntervalTreeNew();
void bfIntervalTreeInitEmpty(BfIntervalTree *intervalTree, BfReal a, BfReal b, BfSize k, BfSize depth);
void bfIntervalTreeDeinit(BfIntervalTree *intervalTree);
void bfIntervalTreeDealloc(BfIntervalTree **intervalTree);
void bfIntervalTreeSetPoints(BfIntervalTree *intervalTree, BfPoints1 const *points, bool rebuildTree);
BfInterval bfIntervalTreeGetInterval(BfIntervalTree const *intervalTree);

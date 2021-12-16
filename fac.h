#pragma once

#include "quadtree.h"

/* This is a temporary type that we'll use until we sort out how we
 * want to organize our lightweight matrix library
 *
 * A `BfFactor` is a block sparse matrix in DOK (dictionary of keys)
 * format. */
typedef struct BfFactor {
  BfSize numBlockRows, numBlockCols;
  BfSize numBlocks;
  BfSize *rowInd, *colInd;
  BfSize *numRows, *numCols;
  BfMat *block;
} BfFactor;

enum BfError
bfMakeFac(BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode,
          BfReal K, BfSize *numFactors, BfFactor **factors);

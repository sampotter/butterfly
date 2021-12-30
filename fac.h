#pragma once

#include "quadtree.h"

/* This is a temporary type that we'll use until we sort out how we
 * want to organize our lightweight matrix library
 *
 * A `BfFactor` is a block sparse matrix in DOK (dictionary of keys)
 * format.
 *
 * The shapes of each block in this structure must be compatible:
 *
 * 1) If `block[i]` and `block[j]` are on the same block row, then
 *    they must have the same number of rows.
 *
 * 2) If `block[i]` and `block[j]` are on the same block column, then
 *    they must have the sae number of columns. */
typedef struct BfFactor {
  /* The number of blocks in each row of the block matrix. */
  BfSize numBlockRows;

  /* The number of blocks in each column of the block matrix. */
  BfSize numBlockCols;

  /* The total number of blocks. */
  BfSize numBlocks;

  /* The row which the block resides in (i.e., `block[i]` is in the
   * `rowInd[i]`th block row). */
  BfSize *rowInd;

  /* The column which the block resides in (i.e., `block[i]` in in the
   * `colInd[i]`th block column). */
  BfSize *colInd;

  /* Array of `numBlockRows + 1` entries containing the offset in rows
   * of each block row in the matrix (i.e., `block[i]` starts on row
   * `rowOffset[rowInd[i]]`). The final entry,
   * `rowOffset[numBlockRows]`, contains a sentinel value indicating
   * the total number of rows in the block matrix. */
  BfSize *rowOffset;

  /* Array of `numBlockCols + 1` entries containing the offset in
   * columns of each block column in the matrix (i.e., `block[i]`
   * starts on column `colOffset[colInd[i]]`). The final entry,
   * `colOffset[numBlockCols]`, contains a sentinel value indicating
   * the total number of columns in the block matrix. */
  BfSize *colOffset;

  /* The blocks in the factor (`numBlocks` total). */
  BfMat *block;

#if BF_DEBUG
  BfPoints2 *srcPtsOrig, *srcPtsEquiv, *tgtPts;
#endif
} BfFactor;

void bfFreeFactor(BfFactor *factor);

void bfMakeFac(BfQuadtree const *tree,
               BfQuadtreeNode const *srcNode, BfQuadtreeNode const *tgtNode,
               BfReal K, BfSize *numFactors, BfFactor **factors);

void bfFreeFac(BfSize numFactors, BfFactor **factorPtr);

void bfMulFac(BfFactor const *factor, BfMat const *X, BfMat *Y);

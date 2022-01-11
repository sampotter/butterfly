#pragma once

#include "mat.h"

struct BfMatBlockCoo {
  BfMat super;

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

  /* The blocks themselves (`numBlocks` total). */
  BfMat **block;
};

BfMatBlockCoo *bfMatBlockCooNew();
void bfMatBlockCooDeinit(BfMatBlockCoo *mat);
void bfMatBlockCooDelete(BfMatBlockCoo **mat);
BfMat *bfMatBlockCooGetMatPtr(BfMatBlockCoo *mat);
BfSize bfMatBlockCooGetNumBlockRows(BfMatBlockCoo const *mat);

#pragma once

#include "mat_block.h"

/** Interface: Mat */

BfMat *bfMatBlockCooCopy(BfMat const *mat);
BfVec *bfMatBlockCooGetRowCopy(BfMat const *mat, BfSize i);
void bfMatBlockCooDelete(BfMat **mat);
BfType bfMatBlockCooGetType(BfMat const *mat);
BfSize bfMatBlockCooNumBytes(BfMat const *mat);
BfSize bfMatBlockCooGetNumRows(BfMat const *mat);
BfSize bfMatBlockCooGetNumCols(BfMat const *mat);
BfMat *bfMatBlockCooGetRowRangeCopy(BfMatBlockCoo const *matBlockCoo, BfSize i0, BfSize i1);
BfMat *bfMatBlockCooMul(BfMat const *op1, BfMat const *op2);
void bfMatBlockCooNegate(BfMat *mat);

/** Interface: MatBlock */

BfSize bfMatBlockCooNumBlocks(BfMatBlock const *mat);

struct BfMatBlockCoo {
  BfMatBlock super;

  BfSize numBlocks;

  /* The row which the block resides in (i.e., `block[i]` is in the
   * `rowInd[i]`th block row). */
  BfSize *rowInd;

  /* The column which the block resides in (i.e., `block[i]` in in the
   * `colInd[i]`th block column). */
  BfSize *colInd;
};

/* Upcasting: */
BfMat *bfMatBlockCooToMat(BfMatBlockCoo *matBlockCoo);
BfMat const *bfMatBlockCooConstToMatConst(BfMatBlockCoo const *matBlockCoo);

/* Downcasting: */
BfMatBlockCoo *bfMatToMatBlockCoo(BfMat *mat);
BfMatBlockCoo const *bfMatConstToMatBlockCooConst(BfMat const *mat);
BfMatBlockCoo const *bfMatBlockConstToMatBlockCooConst(BfMatBlock const *matBlock);

BfMatBlockCoo *bfMatBlockCooNew();
void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows,
                       BfSize numBlockCols, BfSize numBlocks);
void bfMatBlockCooDeinit(BfMatBlockCoo *mat);
void bfMatBlockCooDealloc(BfMatBlockCoo **mat);
void bfMatBlockCooDeinitAndDealloc(BfMatBlockCoo **mat);

#pragma once

#include "mat_block.h"

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

BfMatBlockCoo *bfMatBlockCooNew();
void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows,
                       BfSize numBlockCols, BfSize numBlocks);
BfMat *bfMatBlockCooGetMatPtr(BfMatBlockCoo *mat);

/* BfMat interface */
void bfMatBlockCooDeinit(BfMatBlockCoo *mat);
void bfMatBlockCooDelete(BfMatBlockCoo **mat);
void bfMatBlockCooDeinitAndDelete(BfMatBlockCoo **mat);
BfMatType bfMatBlockCooGetType(BfMatBlockCoo *mat);
BfSize bfMatBlockCooNumBytes(BfMatBlockCoo *mat);
void bfMatBlockCooSave(BfMatBlockCoo const *mat, char const *path);
BfMat *bfMatBlockCooMul(BfMatBlockCoo const *op1, BfMat const *op2);
BfMat *bfMatBlockCooLstSq(BfMatBlockCoo const *lhs, BfMat const *rhs);

/* BfMatBlock interface */
BfSize bfMatBlockCooNumBlocks(BfMatBlockCoo const *mat);

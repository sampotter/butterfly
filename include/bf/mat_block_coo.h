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

BF_DECLARE_INTERFACE_MAT(MatBlockCoo);

BfMatBlockCoo *bfMatBlockCooNew();
void bfMatBlockCooInit(BfMatBlockCoo *mat, BfSize numBlockRows,
                       BfSize numBlockCols, BfSize numBlocks);
BfMat *bfMatBlockCooGetMatPtr(BfMatBlockCoo *mat);

/* BfMatBlock interface */
BfSize bfMatBlockCooNumBlocks(BfMatBlockCoo const *mat);

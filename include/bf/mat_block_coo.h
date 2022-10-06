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

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatBlockCoo)
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DECLARE_INTERFACE(MatBlockCoo)
#undef INTERFACE

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

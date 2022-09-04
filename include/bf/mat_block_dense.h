#pragma once

#include "mat_block.h"

struct BfMatBlockDense {
  BfMatBlock super;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatBlockDense)
#undef INTERFACE

#define INTERFACE BF_INTERFACE_MatBlock
BF_DECLARE_INTERFACE(MatBlockDense)
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatBlockDenseToMat(BfMatBlockDense *matBlockDense);
BfMat const *bfMatBlockDenseConstToMatConst(BfMatBlockDense const *matBlockDense);

/* Downcasting: */
BfMatBlockDense const *bfMatConstToMatBlockDenseConst(BfMat const *mat);

BfMatBlockDense *bfMatBlockDenseNew();
void bfMatBlockDenseInit(BfMatBlockDense *mat, BfSize numBlockRows,
                         BfSize numBlockCols);
void bfMatBlockDenseDeinit(BfMatBlockDense *mat);
void bfMatBlockDenseDealloc(BfMatBlockDense **mat);
void bfMatBlockDenseDeinitAndDealloc(BfMatBlockDense **mat);
BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *mat, BfSize i, BfSize j);
void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j,
                             BfMat *block);

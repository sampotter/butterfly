#pragma once

#include "mat_block.h"

struct BfMatBlockDense {
  BfMatBlock super;
};

BfMatBlockDense *bfMatBlockDenseNew();
void bfMatBlockDenseInit(BfMatBlockDense *mat, BfSize numBlockRows,
                         BfSize numBlockCols);
BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat);

/* BfMat interface */
void bfMatBlockDenseDeinit(BfMatBlockDense *mat);
void bfMatBlockDenseDelete(BfMatBlockDense **mat);
void bfMatBlockDenseDeinitAndDelete(BfMatBlockDense **mat);
BfMatType bfMatBlockDenseGetType(BfMatBlockDense *mat);
BfSize bfMatBlockDenseNumBytes(BfMatBlockDense *mat);
void bfMatBlockDenseSave(BfMatBlockDense const *mat, char const *path);
BfMat *bfMatBlockDenseMul(BfMatBlockDense const *op1, BfMat const *op2);
BfMat *bfMatBlockDenseLstSq(BfMatBlockDense const *lhs, BfMat const *rhs);

/* BfMatBlock interface */
BfSize bfMatBlockDenseNumBlocks(BfMatBlockDense const *mat);

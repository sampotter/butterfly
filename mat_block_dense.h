#pragma once

#include "mat_block.h"

struct BfMatBlockDense {
  BfMatBlock super;
};

BfMatBlockDense *bfMatBlockDenseNew();
void bfMatBlockDenseInit(BfMatBlockDense *mat, BfSize numBlockRows,
                         BfSize numBlockCols);
BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat);
BfMat const *bfMatBlockDenseGetMatConstPtr(BfMatBlockDense const *mat);
BfMat *bfMatBlockDenseGetBlock(BfMatBlockDense *mat, BfSize i, BfSize j);
void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j,
                             BfMat *block);

/* BfMat interface */
void bfMatBlockDenseDeinit(BfMatBlockDense *mat);
void bfMatBlockDenseDelete(BfMatBlockDense **mat);
void bfMatBlockDenseDeinitAndDelete(BfMatBlockDense **mat);
BfMat *bfMatBlockDenseZerosLike(BfMatBlockDense *mat, BfSize numRows, BfSize numCols);
BfMatType bfMatBlockDenseGetType(BfMatBlockDense *mat);
BfSize bfMatBlockDenseNumBytes(BfMatBlockDense *mat);
void bfMatBlockDenseSave(BfMatBlockDense const *mat, char const *path);
BfSize bfMatBlockDenseGetNumRows(BfMatBlockDense const *mat);
BfSize bfMatBlockDenseGetNumCols(BfMatBlockDense const *mat);
BfMatBlockDense *bfMatBlockDenseGetRowRange(BfMatBlockDense *mat, BfSize i0, BfSize i1);
BfMatBlockDense *bfMatBlockDenseGetColRange(BfMatBlockDense *mat, BfSize j0, BfSize j1);
void *bfMatBlockDenseAddInplace(BfMatBlockDense *op1, BfMat *op2);
BfMat *bfMatBlockDenseMul(BfMatBlockDense const *op1, BfMat const *op2);
BfMat *bfMatBlockDenseLstSq(BfMatBlockDense const *lhs, BfMat const *rhs);

/* BfMatBlock interface */
BfSize bfMatBlockDenseNumBlocks(BfMatBlockDense const *mat);

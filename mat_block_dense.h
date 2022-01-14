#pragma once

#include "mat.h"

typedef struct BfMatBlockDense {
  BfMat super;
  BfMat **data;
} BfMatBlockDense;

BfMatBlockDense *bfMatBlockDenseNew();
void bfMatBlockInit(BfMatBlockDense *mat, BfSize numRows, BfSize numCols);
void bfMatBlockDenseDeinit(BfMatBlockDense *mat);
void bfMatBlockDenseDelete(BfMatBlockDense **mat);
void bfMatBlockDenseDeinitAndDelete(BfMatBlockDense **mat);
BfMat *bfMatBlockDenseGetMatPtr(BfMatBlockDense *mat);
BfMatType bfMatBlockDenseGetType(BfMatBlockDense const *mat);
BfSize bfMatBlockDenseNumBytes(BfMatBlockDense const *mat);
void bfMatBlockDenseSave(BfMatBlockDense const *mat, char const *path);
BfMat *bfMatBlockDenseMul(BfMatBlockDense const *op1, BfMat const *op2);
BfMat *bfMatBlockDenseLstSq(BfMatBlockDense const *op1, BfMat const *op2);
void bfMatBlockDenseSetBlock(BfMatBlockDense *mat, BfSize i, BfSize j, BfMat *block);

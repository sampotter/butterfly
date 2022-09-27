#pragma once

#include "mat.h"

struct BfMatDenseReal {
  BfMat super;
  BfSize rowStride;
  BfSize colStride;
  BfReal *data;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatDenseReal)
#undef INTERFACE

BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal);
BfMat const *bfMatDenseRealConstToMatConst(BfMatDenseReal const *matDenseReal);

BfMatDenseReal *bfMatToMatDenseReal(BfMat *mat);
BfMatDenseReal const *bfMatConstToMatDenseRealConst(BfMat const *mat);

BfMatDenseReal *bfMatDenseRealNew();
BfMatDenseReal *bfMatDenseRealFromFile(char const *path, BfSize numRows, BfSize numCols);
void bfMatDenseRealInit(BfMatDenseReal *mat, BfSize numRows, BfSize numCols);
void bfMatDenseRealInitWithValue(BfMatDenseReal *mat, BfSize numRows,
                                 BfSize numCols, BfReal fillValue);
void bfMatDenseRealDeinit(BfMatDenseReal *mat);
void bfMatDenseRealDealloc(BfMatDenseReal **mat);
void bfMatDenseRealDeinitAndDealloc(BfMatDenseReal **mat);

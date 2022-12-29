#pragma once

#include "mat.h"

/** Interface: MatDense */

void bfMatDenseSvd(BfMatDense const *mat, BfMatDense *U, BfMatDiagReal *S, BfMatDense *VH);

typedef struct BfMatDenseVtable {
  __typeof__(&bfMatDenseSvd) Svd;
} BfMatDenseVtable;

/** Implementation: MatDense */

struct BfMatDense {
  BfMat super;
  BfMatDenseVtable *vtable;
  BfSize rowStride;
  BfSize colStride;
};

BfMat *bfMatDenseToMat(BfMatDense *matDense);
BfMat const *bfMatDenseConstToMatConst(BfMatDense const *matDense);

BfMatDense *bfMatToMatDense(BfMat *mat);
BfMatDense const *bfMatConstToMatDenseConst(BfMat const *mat);

BfMatDense *bfMatDenseNew();
void bfMatDenseInit(BfMatDense *matDense, BfMatVtable *matVtable, BfMatDenseVtable *matDenseVtable, BfSize numRows, BfSize numCols, BfSize rowStride, BfSize colStride);
void bfMatDenseDeinit(BfMatDense *matDense);
BfSize bfMatDenseGetRowStride(BfMatDense const *matDense);
BfSize bfMatDenseGetColStride(BfMatDense const *matDense);

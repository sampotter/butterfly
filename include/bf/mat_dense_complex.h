#pragma once

#include "mat.h"

struct BfMatDenseComplex {
  BfMat super;
  BfSize rowStride;
  BfSize colStride;
  BfComplex *data;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatDenseComplex)
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatDenseComplexToMat(BfMatDenseComplex *matDenseComplex);
BfMat const *bfMatDenseComplexConstToMatConst(BfMatDenseComplex const *matDenseComplex);

/* Downcasting: */
BfMatDenseComplex *bfMatToMatDenseComplex(BfMat *mat);
BfMatDenseComplex const *bfMatConstToMatDenseComplexConst(BfMat const *mat);

BfMatDenseComplex *bfMatDenseComplexNew();
BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols);
void bfMatDenseComplexInit(BfMatDenseComplex *mat, BfSize numRows, BfSize numCols);
void bfMatDenseComplexDeinit(BfMatDenseComplex *mat);
void bfMatDenseComplexDealloc(BfMatDenseComplex **mat);
void bfMatDenseComplexDeinitAndDealloc(BfMatDenseComplex **mat);
BfSize bfMatDenseComplexGetRowStride(BfMatDenseComplex const *mat);
BfSize bfMatDenseComplexGetColStride(BfMatDenseComplex const *mat);
void bfMatDenseComplexCopy(BfMatDenseComplex *dst, BfMatDenseComplex const *src);
void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexLstSq(BfMatDenseComplex const *lhs,
                                   BfMatDenseComplex const *rhs);
bool bfMatDenseComplexIsFinite(BfMatDenseComplex const *mat);

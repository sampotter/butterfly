#pragma once

#include "mat.h"

struct BfMatDenseComplex {
  BfMat super;
  BfSize rowStride;
  BfSize colStride;
  BfComplex *data;
};

BF_DECLARE_INTERFACE_MAT(MatDenseComplex);

BfMatDenseComplex *bfMatDenseComplexNew();
BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols);
void bfMatDenseComplexInit(BfMatDenseComplex *mat, BfSize numRows, BfSize numCols);
void bfMatDenseComplexDeinit(BfMatDenseComplex *mat);
void bfMatDenseComplexDealloc(BfMatDenseComplex **mat);
void bfMatDenseComplexDeinitAndDealloc(BfMatDenseComplex **mat);
BfMatDenseComplex *bfMatDenseComplexFromMatPtr(BfMat *mat);
BfMatDenseComplex const *bfMatDenseComplexFromMatConstPtr(BfMat const *mat);
BfMat *bfMatDenseComplexGetMatPtr(BfMatDenseComplex *mat);
BfSize bfMatDenseComplexGetRowStride(BfMatDenseComplex const *mat);
BfSize bfMatDenseComplexGetColStride(BfMatDenseComplex const *mat);
void bfMatDenseComplexCopy(BfMatDenseComplex *dst, BfMatDenseComplex const *src);
BfMat const *bfMatDenseComplexGetMatConstPtr(BfMatDenseComplex const *mat);
void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexLstSq(BfMatDenseComplex const *lhs,
                                   BfMatDenseComplex const *rhs);
bool bfMatDenseComplexIsFinite(BfMatDenseComplex const *mat);

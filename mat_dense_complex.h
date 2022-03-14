#pragma once

#include "mat.h"

struct BfMatDenseComplex {
  BfMat super;
  BfSize rowStride;
  BfSize colStride;
  BfComplex *data;
};

BfMatDenseComplex *bfMatDenseComplexNew();
BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols);
void bfMatDenseComplexInit(BfMatDenseComplex *mat, BfSize numRows, BfSize numCols);
BfMatDenseComplex *bfMatDenseComplexFromMatPtr(BfMat *mat);
BfMatDenseComplex const *bfMatDenseComplexFromMatConstPtr(BfMat const *mat);
BfMat *bfMatDenseComplexGetMatPtr(BfMatDenseComplex *mat);
BfMat const *bfMatDenseComplexGetMatConstPtr(BfMatDenseComplex const *mat);
void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                 BfMatDenseComplex const *op2);
BfMatDenseComplex *
bfMatDenseComplexDenseComplexLstSq(BfMatDenseComplex const *lhs,
                                   BfMatDenseComplex const *rhs);

/* BfMat interface: */
BfMat *bfMatDenseComplexZerosLike(BfMatDenseComplex const *mat, BfSize numRows, BfSize numCols);
void bfMatDenseComplexDeinit(BfMatDenseComplex *mat);
void bfMatDenseComplexDelete(BfMatDenseComplex **mat);
void bfMatDenseComplexDeinitAndDelete(BfMatDenseComplex **mat);
BfMatType bfMatDenseComplexGetType(BfMatDenseComplex const *mat);
BfSize bfMatDenseComplexNumBytes(BfMatDenseComplex const *mat);
void bfMatDenseComplexSave(BfMatDenseComplex const *mat, char const *path);
BfSize bfMatDenseComplexGetNumRows(BfMatDenseComplex const *mat);
BfSize bfMatDenseComplexGetNumCols(BfMatDenseComplex const *mat);
BfMatDenseComplex *bfMatDenseComplexGetRowRange(BfMatDenseComplex *mat, BfSize i0, BfSize i1);
BfMatDenseComplex *bfMatDenseComplexGetColRange(BfMatDenseComplex *mat, BfSize j0, BfSize j1);
void bfMatDenseComplexAddInplace(BfMatDenseComplex *op1, BfMat const *op2);
BfMat *bfMatDenseComplexMul(BfMatDenseComplex const *op1, BfMat const *op2);
BfMat *bfMatDenseComplexLstSq(BfMatDenseComplex const *lhs, BfMat const *rhs);

#pragma once

#include "mat.h"

struct BfMatDenseComplex {
  BfMat super;
  BfSize rowStride;
  BfSize colStride;
  BfComplex *data;
};

BfMatDenseComplex *bfMatDenseComplexNew();
void bfMatDenseComplexInit(BfMatDenseComplex *mat, BfSize numRows, BfSize numCols);
void bfMatDenseComplexDeinit(BfMatDenseComplex *mat);
void bfMatDenseComplexDelete(BfMatDenseComplex **mat);
void bfMatDenseComplexDeinitAndDelete(BfMatDenseComplex **mat);
BfMatDenseComplex *bfMatDenseComplexFromMatPtr(BfMat *mat);
BfMatDenseComplex const *bfMatDenseComplexFromMatConstPtr(BfMat const *mat);
BfMatType bfMatDenseComplexGetType(BfMatDenseComplex const *mat);
BfSize bfMatDenseComplexNumBytes(BfMatDenseComplex const *mat);
void bfMatDenseComplexSave(BfMatDenseComplex const *mat, char const *path);
BfMat *bfMatDenseComplexGetMatPtr(BfMatDenseComplex *mat);
BfMat const *bfMatDenseComplexGetMatConstPtr(BfMatDenseComplex const *mat);
BfMatDenseComplex *bfMatDenseComplexGetColRange(BfMatDenseComplex *mat, BfSize j0, BfSize j1);
BfMat *bfMatDenseComplexMul(BfMatDenseComplex const *op1, BfMat const *op2);
BfMatDenseComplex *bfMatDenseComplexDenseComplexMul(BfMatDenseComplex const *op1,
                                                    BfMatDenseComplex const *op2);
BfMat *bfMatDenseComplexLstSq(BfMatDenseComplex const *lhs, BfMat const *rhs);
BfMatDenseComplex *bfMatDenseComplexDenseComplexLstSq(BfMatDenseComplex const *lhs,
                                                      BfMatDenseComplex const *rhs);
void bfMatDenseComplexSvd(BfMatDenseComplex const *mat, BfMatDenseComplex *U,
                          BfMatDiagReal *S, BfMatDenseComplex *VH);

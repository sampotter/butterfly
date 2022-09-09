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
void bfMatDenseComplexInit(BfMatDenseComplex *mat, BfSize numRows, BfSize numCols);
void bfMatDenseComplexDeinit(BfMatDenseComplex *mat);
void bfMatDenseComplexDealloc(BfMatDenseComplex **mat);
void bfMatDenseComplexDeinitAndDealloc(BfMatDenseComplex **mat);
BfMatDenseComplex *bfMatDenseComplexZeros(BfSize numRows, BfSize numCols);
BfMatDenseComplex *bfMatDenseComplexFromFile(char const *path, BfSize numRows, BfSize numCols);
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

void bf_zmat_add_diag(BfMatDenseComplex *mat, BfReal value);
BfMatDenseComplex *bf_zmat_solve(BfMatDenseComplex *A, BfMatDenseComplex *B);

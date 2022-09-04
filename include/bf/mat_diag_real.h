#pragma once

#include "mat_dense_complex.h"

struct BfMatDiagReal {
  BfMat super;
  BfSize numElts;
  BfReal *data;
};

BfMatDiagReal *bfMatDiagRealNew();
BfMatDiagReal *bfMatDiagRealNewView(BfMatDiagReal *mat);
void bfMatDiagRealInit(BfMatDiagReal *mat, BfSize numRows, BfSize numCols);
void bfMatDiagRealDeinit(BfMatDiagReal *mat);
void bfMatDiagRealDealloc(BfMatDiagReal **mat);
void bfMatDiagRealDeinitAndDealloc(BfMatDiagReal **mat);
BfMat *bfMatDiagRealGetMatPtr(BfMatDiagReal *mat);
BfMat const *bfMatDiagRealGetMatConstPtr(BfMatDiagReal const *mat);
BfMatType bfMatDiagRealGetType(BfMatDiagReal const *mat);
BfSize bfMatDiagRealNumBytes(BfMatDiagReal const *mat);
void bfMatDiagRealSave(BfMatDiagReal const *mat, char const *path);
BfMatDiagReal *bfMatDiagRealGetDiagBlock(BfMatDiagReal *mat, BfSize i0, BfSize i1);
// BfMat *bfMatDiagRealMul(BfMatDiagReal const *op1, BfMat const *op2);
// BfMat *bfMatDiagRealLstSq(BfMatDiagReal const *lhs, BfMat const *rhs);
// BfMat *bfMatDiagRealSolve(BfMatDiagReal const *lhs, BfMat const *rhs);
BfMatDenseComplex *bfMatDiagRealDenseComplexSolve(BfMatDiagReal const *lhs,
                                                  BfMatDenseComplex const *rhs);

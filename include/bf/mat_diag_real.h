#pragma once

#include "mat_dense_complex.h"

BfMat *bfMatDiagRealGetView(BfMat *mat);
BfVec *bfMatDiagRealGetRowCopy(BfMat const *mat, BfSize i);
void bfMatDiagRealDelete(BfMat **mat);
BfType bfMatDiagRealGetType(BfMat const *mat);
bool bfMatDiagRealInstanceOf(BfMat const *mat, BfType type);
BfSize bfMatDiagRealGetNumRows(BfMatDiagReal const *matDiagReal);
BfSize bfMatDiagRealGetNumCols(BfMatDiagReal const *matDiagReal);
BfMat *bfMatDiagRealGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1);
BfVec *bfMatDiagRealMulVec(BfMatDiagReal const *matDiagReal, BfVec const *vec);

struct BfMatDiagReal {
  BfMat super;
  BfSize numElts;
  BfReal *data;
};

/* Upcasting: */
BfMat *bfMatDiagRealToMat(BfMatDiagReal *mat);
BfMat const *bfMatDiagRealConstToMatConst(BfMatDiagReal const *mat);

/* Downcasting: */
BfMatDiagReal *bfMatToMatDiagReal(BfMat *mat);
BfMatDiagReal const *bfMatConstToMatDiagRealConst(BfMat const *mat);

BfMatDiagReal *bfMatDiagRealNew();
BfMatDiagReal *bfMatDiagRealEye(BfSize numRows, BfSize numCols);
BfMatDiagReal *bfMatDiagRealNewConstant(BfSize numRows, BfSize numCols, BfReal diagValue);
BfMatDiagReal *bfMatDiagRealNewFromData(BfSize numRows, BfSize numCols, BfReal const *data);
void bfMatDiagRealInit(BfMatDiagReal *mat, BfSize numRows, BfSize numCols);
void bfMatDiagRealInitView(BfMatDiagReal *mat, BfSize numRows, BfSize numCols, BfReal *data);
void bfMatDiagRealDeinit(BfMatDiagReal *mat);
void bfMatDiagRealDealloc(BfMatDiagReal **mat);
void bfMatDiagRealDeinitAndDealloc(BfMatDiagReal **mat);
void bfMatDiagRealSetConstant(BfMatDiagReal *mat, BfReal value);
BfMatDiagReal *bfMatDiagRealGetDiagBlock(BfMatDiagReal *mat, BfSize i0, BfSize i1);
BfMatDenseComplex *bfMatDiagRealDenseComplexSolve(BfMatDiagReal const *lhs, BfMatDenseComplex const *rhs);
BfVec *bfMatDiagRealGetVecView(BfMatDiagReal *mat);

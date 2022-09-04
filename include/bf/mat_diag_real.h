#pragma once

#include "mat_dense_complex.h"

struct BfMatDiagReal {
  BfMat super;
  BfSize numElts;
  BfReal *data;
};

#define INTERFACE BF_INTERFACE_Mat
BF_DECLARE_INTERFACE(MatDiagReal)
#undef INTERFACE

/* Upcasting: */
BfMat *bfMatDiagRealToMat(BfMatDiagReal *mat);
BfMat const *bfMatDiagRealConstToMatConst(BfMatDiagReal const *mat);

BfMatDiagReal *bfMatDiagRealNew();
BfMatDiagReal *bfMatDiagRealNewView(BfMatDiagReal *mat);
void bfMatDiagRealInit(BfMatDiagReal *mat, BfSize numRows, BfSize numCols);
void bfMatDiagRealDeinit(BfMatDiagReal *mat);
void bfMatDiagRealDealloc(BfMatDiagReal **mat);
void bfMatDiagRealDeinitAndDealloc(BfMatDiagReal **mat);
BfMatDiagReal *bfMatDiagRealGetDiagBlock(BfMatDiagReal *mat, BfSize i0, BfSize i1);
BfMatDenseComplex *bfMatDiagRealDenseComplexSolve(BfMatDiagReal const *lhs,
                                                  BfMatDenseComplex const *rhs);

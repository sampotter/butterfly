#pragma once

BfMat *bfMatDenseRealToMat(BfMatDenseReal *matDenseReal);
BfMat const *bfMatDenseRealConstToMatConst(BfMatDenseReal const *matDenseReal);

#define BF_TO_MAT(_) _Generic((_),                               \
    BfMatDenseComplex *: bfMatDenseComplexToMat,                 \
    BfMatDenseComplex const *: bfMatDenseComplexConstToMatConst, \
    BfMatDenseReal *: bfMatDenseRealToMat,                       \
    BfMatDenseReal const *: bfMatDenseRealConstToMatConst        \
  )(_)

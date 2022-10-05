#include <bf/mat_givens.h>

#include <math.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatGivensComplex)
#undef INTERFACE

/* Interface: Mat */

BF_STUB(BfMat *, MatGivensComplexCopy, BfMat const *)
BF_STUB(BfMat *, MatGivensComplexGetView, BfMat *)
BF_STUB(BfVec *, MatGivensComplexGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatGivensComplexGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatGivensComplexGetColRangeView, BfMat *, BfSize, BfSize, BfSize)
BF_STUB(void, MatGivensComplexDelete, BfMat **)
BF_STUB(BfMat *, MatGivensComplexEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatGivensComplexZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatGivensComplexGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_GIVENS_COMPLEX;
}

BF_STUB(BfSize, MatGivensComplexNumBytes, BfMat const *)
BF_STUB(void, MatGivensComplexSave, BfMat const *, char const *)
BF_STUB(void, MatGivensComplexPrint, BfMat const *, FILE *)

BfSize bfMatGivensComplexGetNumRows(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_GIVENS_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numRows;
  }
}

BfSize bfMatGivensComplexGetNumCols(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_GIVENS_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return BF_SIZE_BAD_VALUE;
  } else {
    return mat->numCols;
  }
}

BF_STUB(void, MatGivensComplexSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatGivensComplexSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatGivensComplexSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatGivensComplexGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatGivensComplexGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatGivensComplexGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatGivensComplexGetColRangeCopy , BfMat const *, BfSize, BfSize)
BF_STUB(void, MatGivensComplexSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatGivensComplexPermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatGivensComplexPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatGivensComplexRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatGivensComplexColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatGivensComplexColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatGivensComplexColNorms, BfMat const *)
BF_STUB(void, MatGivensComplexScaleCols, BfMat *, BfVec const *)
BF_STUB(BfVec *, MatGivensComplexSumCols, BfMat const *)
BF_STUB(void, MatGivensComplexAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatGivensComplexAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatGivensComplexSub, BfMat const *, BfMat const *)
BF_STUB(void, MatGivensComplexSubInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatGivensComplexMul, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatGivensComplexMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatGivensComplexMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatGivensComplexSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatGivensComplexLstSq, BfMat const *, BfMat const *)
BF_STUB(bool, MatGivensComplexIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatGivensComplexBackwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(bool, MatGivensComplexIsZero, BfMat const *)
BF_STUB(void, MatGivensComplexNegate, BfMat *)

/** Upcasting: */

BfMat *bfMatGivensComplexToMat(BfMatGivensComplex *mat) {
  return &mat->super;
}

/** Downcasting: */

BfMatGivensComplex const *bfMatConstToMatGivensComplexConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_GIVENS_COMPLEX)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatGivensComplex const *)mat;
  }
}

/** Implementation: MatGivensComplex */

BfMatGivensComplex *bfMatGivensComplexNew() {
  BEGIN_ERROR_HANDLING();

  BfMatGivensComplex *mat = malloc(sizeof(BfMatGivensComplex));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatGivensComplexInit(BfMatGivensComplex *mat, BfSize n,
                            BfSize srcInd, BfSize elimInd,
                            BfReal c, BfComplex s) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MatVtbl, n, n);
  HANDLE_ERROR();

  mat->srcInd = srcInd;
  mat->elimInd = elimInd;
  mat->c = c;
  mat->s = s;

  END_ERROR_HANDLING()
    bfMatGivensComplexDeinit(mat);
}

void bfMatGivensComplexDeinit(BfMatGivensComplex *mat) {
  mat->srcInd = BF_SIZE_BAD_VALUE;
  mat->elimInd = BF_SIZE_BAD_VALUE;
  mat->c = NAN + I*NAN;
  mat->s = NAN + I*NAN;
}

void bfMatGivensComplexDealloc(BfMatGivensComplex **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatGivensComplexDeinitAndDealloc(BfMatGivensComplex **mat) {
  bfMatGivensComplexDeinit(*mat);
  bfMatGivensComplexDealloc(mat);
}

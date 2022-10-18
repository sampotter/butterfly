#include <bf/mat_zero.h>

#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/vec_zero.h>

/** Interface: Mat */

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatZero)
#undef INTERFACE

BF_STUB(BfMat *, MatZeroCopy, BfMat const *)
BF_STUB(BfMat *, MatZeroGetView, BfMat *)

BfVec *bfMatZeroGetRowCopy(BfMat const *mat, BfSize i) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = bfMatGetNumRows(mat);
  if (i >= numRows)
    RAISE_ERROR(BF_ERROR_INVALID_ARGUMENTS);

  BfSize numCols = bfMatGetNumCols(mat);

  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_ZERO))
    RAISE_ERROR(BF_ERROR_RUNTIME_ERROR);

  BfVecZero *rowCopy = bfVecZeroNew();
  HANDLE_ERROR();

  bfVecZeroInit(rowCopy, numCols);

  END_ERROR_HANDLING() {}

  return bfVecZeroToVec(rowCopy);
}

BF_STUB(BfVec *, MatZeroGetRowView, BfMat *, BfSize)
BF_STUB(BfVec *, MatZeroGetColView, BfMat *, BfSize)
BF_STUB(BfVec *, MatZeroGetColRangeView, BfMat *, BfSize, BfSize, BfSize)
BF_STUB(void, MatZeroDelete, BfMat **)
BF_STUB(BfMat *, MatZeroEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatZeroZerosLike, BfMat const *, BfSize, BfSize)

BfType bfMatZeroGetType(BfMat const *mat) {
  (void)mat;
  return BF_TYPE_MAT_ZERO;
}

BF_STUB(BfSize, MatZeroNumBytes, BfMat const *)
BF_STUB(void, MatZeroSave, BfMat const *, char const *)
BF_STUB(void, MatZeroPrint, BfMat const *, FILE *)

BfSize bfMatZeroGetNumRows(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfSize numRows = BF_SIZE_BAD_VALUE;

  BfMatZero const *matZero = bfMatConstToMatZeroConst(mat);
  HANDLE_ERROR();

  numRows = matZero->super.numRows;

  END_ERROR_HANDLING() {}

  return numRows;
}

BfSize bfMatZeroGetNumCols(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfSize numCols = BF_SIZE_BAD_VALUE;

  BfMatZero const *matZero = bfMatConstToMatZeroConst(mat);
  HANDLE_ERROR();

  numCols = matZero->super.numCols;

  END_ERROR_HANDLING() {}

  return numCols;
}

BF_STUB(void, MatZeroSetRow, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatZeroSetCol, BfMat *, BfSize, BfVec const *)
BF_STUB(void, MatZeroSetColRange, BfMat *, BfSize, BfSize, BfSize, BfVec const *)
BF_STUB(BfMat *, MatZeroGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatZeroGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatZeroGetRowRangeCopy, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatZeroGetColRangeCopy , BfMat const *, BfSize, BfSize)
BF_STUB(void, MatZeroSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(void, MatZeroPermuteRows, BfMat *, BfPerm const *)
BF_STUB(void, MatZeroPermuteCols, BfMat *, BfPerm const *)
BF_STUB(BfVec *, MatZeroRowDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatZeroColDists, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatZeroColDots, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatZeroColNorms, BfMat const *)
BF_STUB(void, MatZeroScaleCols, BfMat *, BfVec const *)
BF_STUB(BfVec *, MatZeroSumCols, BfMat const *)
BF_STUB(void, MatZeroAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatZeroAddDiag, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatZeroSub, BfMat const *, BfMat const *)
BF_STUB(void, MatZeroSubInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatZeroMul, BfMat const *, BfMat const *)
BF_STUB(BfVec *, MatZeroMulVec, BfMat const *, BfVec const *)
BF_STUB(void, MatZeroMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatZeroSolveLU, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatZeroLstSq, BfMat const *, BfMat const *)
BF_STUB(bool, MatZeroIsUpperTri, BfMat const *)
BF_STUB(BfVec *, MatZeroBackwardSolveVec, BfMat const *, BfVec const *)
BF_STUB(bool, MatZeroIsZero, BfMat const *)
BF_STUB(void, MatZeroNegate, BfMat *)
BF_STUB(BfMat *, MatZeroToType, BfMat const *, BfType)

/** Upcasting: */

BfMat *bfMatZeroToMat(BfMatZero *mat) {
  return &mat->super;
}

BfMat const *bfMatZeroConstToMatConst(BfMatZero const *mat) {
  return &mat->super;
}

/** Downcasting: */

BfMatZero const *bfMatConstToMatZeroConst(BfMat const *mat) {
  if (!bfMatInstanceOf(mat, BF_TYPE_MAT_ZERO)) {
    bfSetError(BF_ERROR_TYPE_ERROR);
    return NULL;
  } else {
    return (BfMatZero const *)mat;
  }
}

/** Implementation: MatZero */

BfMatZero *bfMatZeroNew() {
  BEGIN_ERROR_HANDLING();

  BfMatZero *mat = malloc(sizeof(BfMatZero));
  if (mat == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return mat;
}

void bfMatZeroInit(BfMatZero *mat, BfSize numRows, BfSize numCols) {
  BEGIN_ERROR_HANDLING();

  bfMatInit(&mat->super, &MatVtbl, numRows, numCols);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfMatDeinit(&mat->super);
}

void bfMatZeroDeinit(BfMatZero *mat) {
  (void)mat;
}

void bfMatZeroDealloc(BfMatZero **mat) {
  free(*mat);
  *mat = NULL;
}

void bfMatZeroDeinitAndDealloc(BfMatZero **mat) {
  bfMatZeroDeinit(*mat);
  bfMatZeroDealloc(mat);
}

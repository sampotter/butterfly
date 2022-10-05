#include <bf/mat.h>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/mat_dense_real.h>

/** Interface: Mat */

BfMat *bfMatCopy(BfMat const *mat) {
  return mat->vtbl->Copy(mat);
}

BfMat *bfMatGetView(BfMat *mat) {
  return mat->vtbl->GetView(mat);
}

BfVec *bfMatGetRowView(BfMat *mat, BfSize i) {
  return mat->vtbl->GetRowView(mat, i);
}

BfVec *bfMatGetColView(BfMat *mat, BfSize j) {
  return mat->vtbl->GetColView(mat, j);
}

BfVec *bfMatGetColRangeView(BfMat *mat, BfSize i0, BfSize i1, BfSize j) {
  return mat->vtbl->GetColRangeView(mat, i0, i1, j);
}

void bfMatDelete(BfMat **mat) {
  (*mat)->vtbl->Delete(mat);
}

BfMat *bfMatEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  return mat->vtbl->EmptyLike(mat, numRows, numCols);
}

BfMat *bfMatZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  return mat->vtbl->ZerosLike(mat, numRows, numCols);
}

BfType bfMatGetType(BfMat const *mat) {
  return mat->vtbl->GetType(mat);
}

BfSize bfMatNumBytes(BfMat const *mat) {
  return mat->vtbl->NumBytes(mat);
}

void bfMatSave(BfMat const *mat, char const *path) {
  mat->vtbl->Save(mat, path);
}

void bfMatPrint(BfMat const *mat, FILE *fp) {
  mat->vtbl->Print(mat, fp);
}

BfSize bfMatGetNumRows(BfMat const *mat) {
  return mat->vtbl->GetNumRows(mat);
}

BfSize bfMatGetNumCols(BfMat const *mat) {
  return mat->vtbl->GetNumCols(mat);
}

void bfMatSetRow(BfMat *mat, BfSize i, BfVec const *rowVec) {
  mat->vtbl->SetRow(mat, i, rowVec);
}

void bfMatSetCol(BfMat *mat, BfSize j, BfVec const *colVec) {
  mat->vtbl->SetCol(mat, j, colVec);
}

void bfMatSetColRange(BfMat *mat, BfSize j, BfSize i0, BfSize i1, BfVec const *colVec) {
  mat->vtbl->SetColRange(mat, j, i0, i1, colVec);
}

BfMat *bfMatGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  return mat->vtbl->GetRowRange(mat, i0, i1);
}

BfMat *bfMatGetColRange(BfMat *mat, BfSize i0, BfSize i1) {
  return mat->vtbl->GetColRange(mat, i0, i1);
}

BfMat *bfMatGetRowRangeCopy(BfMat const *mat, BfSize i0, BfSize i1) {
  return mat->vtbl->GetRowRangeCopy(mat, i0, i1);
}

BfMat *bfMatGetColRangeCopy(BfMat const *mat, BfSize j0, BfSize j1) {
  return mat->vtbl->GetColRangeCopy(mat, j0, j1);
}

void bfMatSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows) {
  mat->vtbl->SetRowRange(mat, i0, i1, rows);
}

void bfMatPermuteRows(BfMat *mat, BfPerm const *perm) {
  mat->vtbl->PermuteRows(mat, perm);
}

void bfMatPermuteCols(BfMat *mat, BfPerm const *perm) {
  mat->vtbl->PermuteCols(mat, perm);
}

BfVec *bfMatRowDists(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->RowDists(mat, otherMat);
}

BfVec *bfMatColDists(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->ColDists(mat, otherMat);
}

BfVec *bfMatColDots(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->ColDots(mat, otherMat);
}

BfVec *bfMatColNorms(BfMat const *mat) {
  return mat->vtbl->ColNorms(mat);
}

void bfMatScaleCols(BfMat *mat, BfVec const *vec) {
  mat->vtbl->ScaleCols(mat, vec);
}

BfVec *bfMatSumCols(BfMat const *mat) {
  return mat->vtbl->SumCols(mat);
}

void bfMatAddInplace(BfMat *lhs, BfMat const *rhs) {
  lhs->vtbl->AddInplace(lhs, rhs);
}

void bfMatAddDiag(BfMat *mat, BfMat const *diagMat) {
  mat->vtbl->AddDiag(mat, diagMat);
}

BfMat *bfMatSub(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->Sub(mat, otherMat);
}

void bfMatSubInplace(BfMat *mat, BfMat const *otherMat) {
  mat->vtbl->SubInplace(mat, otherMat);
}

BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->Mul(lhs, rhs);
}

BfVec *bfMatMulVec(BfMat const *mat, BfVec const *vec) {
  return mat->vtbl->MulVec(mat, vec);
}

void bfMatMulInplace(BfMat *mat, BfMat const *otherMat) {
  mat->vtbl->MulInplace(mat, otherMat);
}

BfMat *bfMatSolveLU(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->SolveLU(mat, otherMat);
}

BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->LstSq(lhs, rhs);
}

bool bfMatIsUpperTri(BfMat const *mat) {
  return mat->vtbl->IsUpperTri(mat);
}

BfVec *bfMatBackwardSolveVec(BfMat const *mat, BfVec const *vec) {
  return mat->vtbl->BackwardSolveVec(mat, vec);
}

bool bfMatIsZero(BfMat const *mat) {
  return mat->vtbl->IsZero(mat);
}

void bfMatNegate(BfMat *mat) {
  mat->vtbl->Negate(mat);
}

/** Implementation: Mat */

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols) {
  mat->vtbl = vtbl;
  mat->props = BF_MAT_PROPS_NONE;
  mat->numRows = numRows;
  mat->numCols = numCols;
#if BF_DEBUG
  mat->aux = NULL;
#endif
}

void bfMatDeinit(BfMat *mat) {
#if BF_DEBUG
  mat->vtbl = NULL;
  mat->props = BF_MAT_PROPS_NONE;
  mat->numRows = BF_SIZE_BAD_VALUE;
  mat->numCols = BF_SIZE_BAD_VALUE;
  mat->aux = NULL;
#endif
}

BfMat *bfMatFromFile(char const *path, BfSize numRows, BfSize numCols, BfDtype dtype) {
  switch (dtype) {
  case BF_DTYPE_COMPLEX:
    return (BfMat *)bfMatDenseComplexFromFile(path, numRows, numCols);
  case BF_DTYPE_REAL:
    return (BfMat *)bfMatDenseRealFromFile(path, numRows, numCols);
  default:
    bfSetError(BF_ERROR_NOT_IMPLEMENTED);
    return NULL;
  }
}

bool bfMatInstanceOf(BfMat const *mat, BfType type) {
  return bfTypeDerivedFrom(bfMatGetType(mat), type);
}

bool bfMatIsTransposed(BfMat const *mat) {
  return mat->props & BF_MAT_PROPS_TRANS;
}

BfMat *bfMatConjTrans(BfMat *mat) {
  mat->props ^= (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ);
  return mat;
}

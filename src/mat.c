#include <bf/mat.h>

#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat_dense_complex.h>
#include <bf/mat_dense_real.h>

BfMat *bfMatGetView(BfMat *mat) {
  return mat->vtbl->GetView(mat);
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

enum BfMatTypes bfMatGetType(BfMat const *mat) {
  return mat->vtbl->GetType(mat);
}

bool bfMatInstanceOf(BfMat const *mat, BfMatType matType) {
  return bfMatTypeDerivedFrom(bfMatGetType(mat), matType);
}

BfSize bfMatNumBytes(BfMat const *mat) {
  return mat->vtbl->NumBytes(mat);
}

void bfMatSave(BfMat const *mat, char const *path) {
  mat->vtbl->Save(mat, path);
}

void bfMatPrint(FILE *fp, BfMat const *mat) {
  mat->vtbl->Print(fp, mat);
}

BfSize bfMatGetNumRows(BfMat const *mat) {
  return mat->vtbl->GetNumRows(mat);
}

BfSize bfMatGetNumCols(BfMat const *mat) {
  return mat->vtbl->GetNumCols(mat);
}

BfMat *bfMatGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  return mat->vtbl->GetRowRange(mat, i0, i1);
}

BfMat *bfMatGetColRange(BfMat *mat, BfSize i0, BfSize i1) {
  return mat->vtbl->GetColRange(mat, i0, i1);
}

void bfMatSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows) {
  mat->vtbl->SetRowRange(mat, i0, i1, rows);
}

BfMat *bfMatRowDists(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->RowDists(mat, otherMat);
}

BfMat *bfMatColDists(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->ColDists(mat, otherMat);
}

void bfMatScaleCols(BfMat *mat, BfMat const *otherMat) {
  mat->vtbl->ScaleCols(mat, otherMat);
}

BfMat *bfMatSumCols(BfMat const *mat) {
  return mat->vtbl->SumCols(mat);
}

void bfMatAddInplace(BfMat *lhs, BfMat const *rhs) {
  lhs->vtbl->AddInplace(lhs, rhs);
}

void bfMatAddDiag(BfMat *mat, BfMat const *diagMat) {
  mat->vtbl->AddDiag(mat, diagMat);
}

BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->Mul(lhs, rhs);
}

void bfMatMulInplace(BfMat *mat, BfMat const *otherMat) {
  mat->vtbl->MulInplace(mat, otherMat);
}

BfMat *bfMatSolve(BfMat const *mat, BfMat const *otherMat) {
  return mat->vtbl->Solve(mat, otherMat);
}

BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->LstSq(lhs, rhs);
}

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

bool bfMatIsTransposed(BfMat const *mat) {
  return mat->props & BF_MAT_PROPS_TRANS;
}

BfMat *bfMatConjTrans(BfMat *mat) {
  mat->props ^= (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ);
  return mat;
}

#include <bf/mat.h>

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols) {
  mat->vtbl = vtbl;
  mat->props = BF_MAT_PROPS_NONE;
  mat->numRows = numRows;
  mat->numCols = numCols;
#if BF_DEBUG
  mat->aux = NULL;
#endif
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

void bfMatDeinit(BfMat *mat) {
#if BF_DEBUG
  mat->vtbl = NULL;
  mat->props = BF_MAT_PROPS_NONE;
  mat->numRows = BF_SIZE_BAD_VALUE;
  mat->numCols = BF_SIZE_BAD_VALUE;
  mat->aux = NULL;
#endif
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

bool bfMatIsTransposed(BfMat const *mat) {
  return mat->props & BF_MAT_PROPS_TRANS;
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

BfMat *bfMatConjTrans(BfMat *mat) {
  mat->props ^= (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ);
  return mat;
}

void bfMatAddInplace(BfMat *lhs, BfMat const *rhs) {
  lhs->vtbl->AddInplace(lhs, rhs);
}

BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->Mul(lhs, rhs);
}

BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->LstSq(lhs, rhs);
}

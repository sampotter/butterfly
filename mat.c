#include "mat.h"

void bfMatInit(BfMat *mat, BfMatVtable *vtbl, BfSize numRows, BfSize numCols) {
  mat->vtbl = vtbl;
  mat->props = BF_MAT_PROPS_NONE;
  mat->numRows = numRows;
  mat->numCols = numCols;
#if BF_DEBUG
  mat->aux = NULL;
#endif
}

BfMat *bfMatEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  return mat->vtbl->emptyLike(mat, numRows, numCols);
}

BfMat *bfMatZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  return mat->vtbl->zerosLike(mat, numRows, numCols);
}

void bfMatDeinit(BfMat *mat) {
  mat->vtbl->deinit(mat);

#if BF_DEBUG
  mat->vtbl = NULL;
  mat->props = BF_MAT_PROPS_NONE;
  mat->numRows = BF_SIZE_BAD_VALUE;
  mat->numCols = BF_SIZE_BAD_VALUE;
  mat->aux = NULL;
#endif
}

void bfMatDelete(BfMat **mat) {
  (*mat)->vtbl->delete(mat);
}

void bfMatDeinitAndDelete(BfMat **mat) {
  (*mat)->vtbl->deinitAndDelete(mat);
}

enum BfMatTypes bfMatGetType(BfMat const *mat) {
  return mat->vtbl->getType(mat);
}

BfSize bfMatNumBytes(BfMat const *mat) {
  return mat->vtbl->numBytes(mat);
}

void bfMatSave(BfMat const *mat, char const *path) {
  mat->vtbl->save(mat, path);
}

bool bfMatIsTransposed(BfMat const *mat) {
  return mat->props & BF_MAT_PROPS_TRANS;
}

BfSize bfMatGetNumRows(BfMat const *mat) {
  return mat->vtbl->getNumRows(mat);
}

BfSize bfMatGetNumCols(BfMat const *mat) {
  return mat->vtbl->getNumCols(mat);
}

BfMat *bfMatGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  return mat->vtbl->getRowRange(mat, i0, i1);
}

BfMat *bfMatGetColRange(BfMat *mat, BfSize i0, BfSize i1) {
  return mat->vtbl->getColRange(mat, i0, i1);
}

void bfMatSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *rows) {
  mat->vtbl->setRowRange(mat, i0, i1, rows);
}

BfMat *bfMatConjTrans(BfMat *mat) {
  mat->props ^= (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ);
  return mat;
}

void bfMatAddInplace(BfMat *lhs, BfMat const *rhs) {
  lhs->vtbl->addInplace(lhs, rhs);
}

BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->mul(lhs, rhs);
}

BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->lstSq(lhs, rhs);
}

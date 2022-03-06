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
  return bfMatIsTransposed(mat) ? mat->numCols : mat->numRows;
}

BfSize bfMatGetNumCols(BfMat const *mat) {
  return bfMatIsTransposed(mat) ? mat->numRows : mat->numCols;
}

BfMat *bfMatConjTrans(BfMat *mat) {
  mat->props ^= (BF_MAT_PROPS_TRANS | BF_MAT_PROPS_CONJ);
  return mat;
}

BfMat *bfMatMul(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->mul(lhs, rhs);
}

BfMat *bfMatLstSq(BfMat const *lhs, BfMat const *rhs) {
  return lhs->vtbl->lstSq(lhs, rhs);
}

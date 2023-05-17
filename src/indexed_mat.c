#include <bf/indexed_mat.h>

#include <bf/assert.h>
#include <bf/error.h>
#include <bf/error_macros.h>
#include <bf/mat.h>
#include <bf/mem.h>

BfIndexedMat *bfIndexedMatNewFromMat(BfSize i0, BfSize j0, BfMat *mat, BfPolicy policy) {
  BF_ERROR_BEGIN();

  BfIndexedMat *indexedMat = bfMemAlloc(1, sizeof(BfIndexedMat));
  HANDLE_ERROR();

  bfIndexedMatInitFromMat(indexedMat, i0, j0, mat, policy);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }

  return indexedMat;
}

void bfIndexedMatInitFromMat(BfIndexedMat *indexedMat, BfSize i0, BfSize j0, BfMat *mat, BfPolicy policy) {
  BF_ERROR_BEGIN();

  indexedMat->i0 = i0;
  indexedMat->j0 = j0;

  indexedMat->mat = bfMatGet(mat, policy);
  HANDLE_ERROR();

  BF_ERROR_END() {
    BF_DIE();
  }
}

void bfIndexedMatDeinit(BfIndexedMat *indexedMat) {
  indexedMat->i0 = BF_SIZE_BAD_VALUE;
  indexedMat->j0 = BF_SIZE_BAD_VALUE;
  bfMatDelete(&indexedMat->mat);
}

void bfIndexedMatDealloc(BfIndexedMat **indexedMat) {
  bfMemFree(*indexedMat);
  *indexedMat = NULL;
}

void bfIndexedMatDelete(BfIndexedMat **indexedMat) {
  bfIndexedMatDeinit(*indexedMat);
  bfIndexedMatDealloc(indexedMat);
}

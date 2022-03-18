#include "mat_product.h"

#include <assert.h>
#include <stdlib.h>

#include "error.h"
#include "error_macros.h"

BF_DEFINE_MAT_VTABLE(MatProduct);

BfMatProduct *bfMatProductNew() {
  BEGIN_ERROR_HANDLING();

  BfMatProduct *prod = malloc(sizeof(BfMatProduct));
  if (prod == NULL)
    RAISE_ERROR(BF_ERROR_MEMORY_ERROR);

  END_ERROR_HANDLING() {}

  return prod;
}

void bfMatProductInit(BfMatProduct *mat) {
  BEGIN_ERROR_HANDLING();

  /* We don't store the number of rows or columns in `mat->super`
   * since we always look up the number of rows and columns from the
   * leftmost and rightmost factors at runtime. */
  bfMatInit(&mat->super, &matVtbl, BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE);

  mat->factorArr = bfGetUninitializedPtrArray();

  bfInitPtrArray(&mat->factorArr, /* capacity: */ 4);
  HANDLE_ERROR();

  END_ERROR_HANDLING()
    bfPtrArrayDeinit(&mat->factorArr);
}

BfMat *bfMatProductGetMatPtr(BfMatProduct *mat) {
  return &mat->super;
}

BfMat const *bfMatProductGetMatConstPtr(BfMatProduct const *mat) {
  return &mat->super;
}

BfMatProduct *bfMatProductEmptyLike(BfMatProduct const *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

BfMatProduct *bfMatProductZerosLike(BfMatProduct const *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

void bfMatProductDeinit(BfMatProduct *prod) {
  for (BfSize i = 0; i < bfPtrArraySize(&prod->factorArr); ++i) {
    BfMat *mat = bfPtrArrayGet(&prod->factorArr, i);
    bfMatDeinitAndDelete(&mat);
  }

  bfPtrArrayDeinit(&prod->factorArr);
}

void bfMatProductDelete(BfMatProduct **prod) {
  free(*prod);
  *prod = NULL;
}

void bfMatProductDeinitAndDelete(BfMatProduct **prod) {
  bfMatProductDeinit(*prod);
  bfMatProductDelete(prod);
}

BfMatType bfMatProductGetType(BfMatProduct const *mat) {
  return BF_MAT_TYPE_PRODUCT;
}

BfSize bfMatProductNumBytes(BfMatProduct const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatProductSave(BfMatProduct const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfSize bfMatProductGetNumRows(BfMatProduct const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMat const *factor = NULL;
  BfSize numRows = BF_SIZE_BAD_VALUE;

  if (bfMatIsTransposed(bfMatProductGetMatConstPtr(mat)))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  factor = bfMatProductGetFactor((BfMatProduct *)mat, 0);
  HANDLE_ERROR();

  numRows = bfMatGetNumRows(factor);

  END_ERROR_HANDLING() {}

  return numRows;
}

BfSize bfMatProductGetNumCols(BfMatProduct const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMat const *factor = NULL;
  BfSize numCols = BF_SIZE_BAD_VALUE;

  if (bfMatIsTransposed(bfMatProductGetMatConstPtr(mat)))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numFactors = bfMatProductNumFactors((BfMatProduct *)mat);

  factor = bfMatProductGetFactor((BfMatProduct *)mat, numFactors - 1);
  HANDLE_ERROR();

  numCols = bfMatGetNumCols(factor);

  END_ERROR_HANDLING() {}

  return numCols;
}

BfMatProduct *bfMatProductGetRowRange(BfMatProduct *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

BfMatProduct *bfMatProductGetColRange(BfMatProduct *, BfSize, BfSize) {
  assert(false);
  return NULL;
}

void bfMatProductSetRowRange(BfMatProduct *, BfSize, BfSize, BfMat const *) {
  assert(false);
}

void bfMatProductAddInplace(BfMatProduct *, BfMat const *) {
  assert(false);
}

BfMat *bfMatProductMul(BfMatProduct const *op1, BfMat const *op2) {
  /* TODO: add error handling... */
  BfSize numFactors = bfMatProductNumFactors((BfMatProduct *)op1);
  BfMat *factor = NULL, *result = (BfMat *)op2, *prev = NULL;
  BfSize i = numFactors - 1;
  factor = bfMatProductGetFactor((BfMatProduct *)op1, i);
  prev = bfMatMul(factor, op2);
  while (i > 0) {
    factor = bfMatProductGetFactor((BfMatProduct *)op1, --i);
    result = bfMatMul(factor, prev);
    bfMatDeinitAndDelete(&prev);
    prev = result;
  }
  return result;
}

BfMat *bfMatProductLstSq(BfMatProduct const *, BfMat const *) {
  assert(false);
  return NULL;
}

BfSize bfMatProductNumFactors(BfMatProduct *prod) {
  return bfPtrArraySize(&prod->factorArr);
}

BfMat *bfMatProductGetFactor(BfMatProduct *prod, BfSize i) {
  return bfPtrArrayGet(&prod->factorArr, i);
}

void bfMatProductPostMultiply(BfMatProduct *prod, BfMat *mat) {
  // TODO: verify that shapes are compatible
  bfPtrArrayAppend(&prod->factorArr, mat);
}

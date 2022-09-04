#include <bf/mat_product.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatProduct)
#undef INTERFACE

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
  bfMatInit(&mat->super, &MatVtbl, BF_SIZE_BAD_VALUE, BF_SIZE_BAD_VALUE);

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

void bfMatProductDelete(BfMat **mat) {
  bfMatProductDeinitAndDealloc((BfMatProduct **)mat);
}

BfMat *bfMatProductEmptyLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  (void)mat;
  (void)numRows;
  (void)numCols;
  assert(false);
  return NULL;
}

BfMat *bfMatProductZerosLike(BfMat const *mat, BfSize numRows, BfSize numCols) {
  (void)mat;
  (void)numRows;
  (void)numCols;
  assert(false);
  return NULL;
}

void bfMatProductDeinit(BfMatProduct *prod) {
  for (BfSize i = 0; i < bfPtrArraySize(&prod->factorArr); ++i) {
    BfMat *mat = bfPtrArrayGet(&prod->factorArr, i);
    bfMatDelete(&mat);
  }

  bfPtrArrayDeinit(&prod->factorArr);
}

void bfMatProductDealloc(BfMatProduct **prod) {
  free(*prod);
  *prod = NULL;
}

void bfMatProductDeinitAndDealloc(BfMatProduct **prod) {
  bfMatProductDeinit(*prod);
  bfMatProductDealloc(prod);
}

BfMatType bfMatProductGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_PRODUCT;
}

bool bfMatProductInstanceOf(BfMat const *mat, BfMatType matType) {
  BfMat const *parent = bfMatProductGetMatConstPtr((BfMatProduct const *)mat);
  return bfMatTypeDerivedFrom(bfMatGetType(parent), matType);
}

BfSize bfMatProductNumBytes(BfMat const *mat) {
  (void)mat;
  assert(false);
  return BF_SIZE_BAD_VALUE;
}

void bfMatProductSave(BfMat const *mat, char const *path) {
  (void)mat;
  (void)path;
  assert(false);
}

BfSize bfMatProductGetNumRows(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMat const *factor = NULL;
  BfSize numRows = BF_SIZE_BAD_VALUE;

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  factor = bfMatProductGetFactor((BfMatProduct *)mat, 0);
  HANDLE_ERROR();

  numRows = bfMatGetNumRows(factor);

  END_ERROR_HANDLING() {}

  return numRows;
}

BfSize bfMatProductGetNumCols(BfMat const *mat) {
  BEGIN_ERROR_HANDLING();

  BfMat const *factor = NULL;
  BfSize numCols = BF_SIZE_BAD_VALUE;

  if (bfMatIsTransposed(mat))
    RAISE_ERROR(BF_ERROR_NOT_IMPLEMENTED);

  BfSize numFactors = bfMatProductNumFactors((BfMatProduct *)mat);

  factor = bfMatProductGetFactor((BfMatProduct *)mat, numFactors - 1);
  HANDLE_ERROR();

  numCols = bfMatGetNumCols(factor);

  END_ERROR_HANDLING() {}

  return numCols;
}

BfMat *bfMatProductGetRowRange(BfMat *mat, BfSize i0, BfSize i1) {
  (void)mat; (void)i0; (void)i1;
  assert(false);
  return NULL;
}

BfMat *bfMatProductGetColRange(BfMat *mat, BfSize j0, BfSize j1) {
  (void)mat; (void)j0; (void)j1;
  assert(false);
  return NULL;
}

void bfMatProductSetRowRange(BfMat *mat, BfSize i0, BfSize i1, BfMat const *otherMat) {
  (void)mat; (void)i0; (void)i1; (void)otherMat;
  assert(false);
}

void bfMatProductAddInplace(BfMat *mat, BfMat const *otherMat) {
  (void)mat; (void)otherMat;
  assert(false);
}

BfMat *bfMatProductMul(BfMat const *op1, BfMat const *op2) {
  /* TODO: add error handling... */
  BfSize numFactors = bfMatProductNumFactors((BfMatProduct *)op1);
  BfMat *factor = NULL, *result = (BfMat *)op2, *prev = NULL;
  BfSize i = numFactors - 1;
  factor = bfMatProductGetFactor((BfMatProduct *)op1, i);
  prev = bfMatMul(factor, op2);
  while (i > 0) {
    factor = bfMatProductGetFactor((BfMatProduct *)op1, --i);
    result = bfMatMul(factor, prev);
    bfMatDelete(&prev);
    prev = result;
  }
  return result;
}

BfMat *bfMatProductLstSq(BfMat const *mat, BfMat const *otherMat) {
  (void)mat; (void)otherMat;
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

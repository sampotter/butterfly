#include <bf/mat_product.h>

#include <assert.h>
#include <stdlib.h>

#include <bf/error.h>
#include <bf/error_macros.h>

#define INTERFACE BF_INTERFACE_Mat
BF_DEFINE_VTABLE(Mat, MatProduct)
#undef INTERFACE

BF_STUB(BfMat *, MatProductGetView, BfMat *)

void bfMatProductDelete(BfMat **mat) {
  bfMatProductDeinitAndDealloc((BfMatProduct **)mat);
}

BF_STUB(BfMat *, MatProductEmptyLike, BfMat const *, BfSize, BfSize)
BF_STUB(BfMat *, MatProductZerosLike, BfMat const *, BfSize, BfSize)

BfMatType bfMatProductGetType(BfMat const *mat) {
  (void)mat;
  return BF_MAT_TYPE_PRODUCT;
}

bool bfMatProductInstanceOf(BfMat const *mat, BfMatType matType) {
  return bfMatTypeDerivedFrom(bfMatGetType(mat), matType);
}

BF_STUB(BfSize, MatProductNumBytes, BfMat const *)
BF_STUB(void, MatProductSave, BfMat const *, char const *)
BF_STUB(void, MatProductPrint, FILE *, BfMat const *)

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

BF_STUB(BfMat *, MatProductGetRowRange, BfMat *, BfSize, BfSize)
BF_STUB(BfMat *, MatProductGetColRange, BfMat *, BfSize, BfSize)
BF_STUB(void, MatProductSetRowRange, BfMat *, BfSize, BfSize, BfMat const *)
BF_STUB(BfMat *, MatProductRowDists, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatProductColDists, BfMat const *, BfMat const *)
BF_STUB(void, MatProductScaleCols, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatProductSumCols, BfMat const *)
BF_STUB(void, MatProductAddInplace, BfMat *, BfMat const *)
BF_STUB(void, MatProductAddDiag, BfMat *, BfMat const *)

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

BF_STUB(void, MatProductMulInplace, BfMat *, BfMat const *)
BF_STUB(BfMat *, MatProductSolve, BfMat const *, BfMat const *)
BF_STUB(BfMat *, MatProductLstSq, BfMat const *, BfMat const *)

/* Upcasting: */

BfMat *bfMatProductToMat(BfMatProduct *matProduct) {
  return &matProduct->super;
}

/* Implementation: MatProduct */

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
